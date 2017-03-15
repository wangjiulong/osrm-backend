#include "partition/partitioner.hpp"
#include "partition/bisection_graph.hpp"
#include "partition/bisection_to_partition.hpp"
#include "partition/cell_storage.hpp"
#include "partition/compressed_node_based_graph_reader.hpp"
#include "partition/edge_based_graph_reader.hpp"
#include "partition/io.hpp"
#include "partition/multi_level_partition.hpp"
#include "partition/recursive_bisection.hpp"

#include "extractor/io.hpp"

#include "util/coordinate.hpp"
#include "util/geojson_debug_logger.hpp"
#include "util/geojson_debug_policies.hpp"
#include "util/integer_range.hpp"
#include "util/json_container.hpp"
#include "util/log.hpp"

#include <algorithm>
#include <iterator>
#include <vector>

#include <boost/assert.hpp>

#include "util/geojson_debug_logger.hpp"
#include "util/geojson_debug_policies.hpp"
#include "util/json_container.hpp"
#include "util/timing_util.hpp"

namespace osrm
{
namespace partition
{

void LogGeojson(const std::string &filename, const std::vector<std::uint32_t> &bisection_ids)
{
    // reload graph, since we destroyed the old one
    auto compressed_node_based_graph = LoadCompressedNodeBasedGraph(filename);

    util::Log() << "Loaded compressed node based graph: "
                << compressed_node_based_graph.edges.size() << " edges, "
                << compressed_node_based_graph.coordinates.size() << " nodes";

    groupEdgesBySource(begin(compressed_node_based_graph.edges),
                       end(compressed_node_based_graph.edges));

    auto graph =
        makeBisectionGraph(compressed_node_based_graph.coordinates,
                           adaptToBisectionEdge(std::move(compressed_node_based_graph.edges)));

    const auto get_level = [](const std::uint32_t lhs, const std::uint32_t rhs) {
        auto xored = lhs ^ rhs;
        std::uint32_t level = log(xored) / log(2.0);
        return level;
    };

    std::vector<std::vector<util::Coordinate>> border_vertices(33);

    for (NodeID nid = 0; nid < graph.NumberOfNodes(); ++nid)
    {
        const auto source_id = bisection_ids[nid];
        for (const auto &edge : graph.Edges(nid))
        {
            const auto target_id = bisection_ids[edge.target];
            if (source_id != target_id)
            {
                auto level = get_level(source_id, target_id);
                border_vertices[level].push_back(graph.Node(nid).coordinate);
                border_vertices[level].push_back(graph.Node(edge.target).coordinate);
            }
        }
    }

    util::ScopedGeojsonLoggerGuard<util::CoordinateVectorToMultiPoint> guard(
        "border_vertices.geojson");
    std::size_t level = 0;
    for (auto &bv : border_vertices)
    {
        if (!bv.empty())
        {
            std::sort(bv.begin(), bv.end(), [](const auto lhs, const auto rhs) {
                return std::tie(lhs.lon, lhs.lat) < std::tie(rhs.lon, rhs.lat);
            });
            bv.erase(std::unique(bv.begin(), bv.end()), bv.end());

            util::json::Object jslevel;
            jslevel.values["level"] = util::json::Number(level++);
            guard.Write(bv, jslevel);
        }
    }
}

int Partitioner::Run(const PartitionConfig &config)
{
    auto compressed_node_based_graph =
        LoadCompressedNodeBasedGraph(config.compressed_node_based_graph_path.string());

    util::Log() << "Loaded compressed node based graph: "
                << compressed_node_based_graph.edges.size() << " edges, "
                << compressed_node_based_graph.coordinates.size() << " nodes";

    groupEdgesBySource(begin(compressed_node_based_graph.edges),
                       end(compressed_node_based_graph.edges));

    auto graph =
        makeBisectionGraph(compressed_node_based_graph.coordinates,
                           adaptToBisectionEdge(std::move(compressed_node_based_graph.edges)));

    util::Log() << " running partition: " << config.minimum_cell_size << " " << config.balance
                << " " << config.boundary_factor << " " << config.num_optimizing_cuts << " "
                << config.small_component_size
                << " # max_cell_size balance boundary cuts small_component_size";
    RecursiveBisection recursive_bisection(graph,
                                           config.minimum_cell_size,
                                           config.balance,
                                           config.boundary_factor,
                                           config.num_optimizing_cuts,
                                           config.small_component_size);

    // Up until now we worked on the compressed node based graph.
    // But what we actually need is a partition for the edge based graph to work on.
    // The following loads a mapping from node based graph to edge based graph.
    // Then loads the edge based graph tanslates the partition and modifies it.
    // For details see #3205

    std::vector<extractor::NBGToEBG> mapping;
    extractor::io::read(config.cnbg_ebg_mapping_path.string(), mapping);
    util::Log() << "Loaded node based graph to edge based graph mapping";

    auto edge_based_graph = LoadEdgeBasedGraph(config.edge_based_graph_path.string());
    util::Log() << "Loaded edge based graph for mapping partition ids: "
                << edge_based_graph->GetNumberOfEdges() << " edges, "
                << edge_based_graph->GetNumberOfNodes() << " nodes";

    // TODO: node based graph to edge based graph partition id mapping should be done split off.

    // Partition ids, keyed by node based graph nodes
    const auto &node_based_partition_ids = recursive_bisection.BisectionIDs();

    // Partition ids, keyed by edge based graph nodes
    std::vector<NodeID> edge_based_partition_ids(edge_based_graph->GetNumberOfNodes(),
                                                 SPECIAL_NODEID);

    // Only resolve all easy cases in the first pass
    for (const auto &entry : mapping)
    {
        const auto u = entry.u;
        const auto v = entry.v;
        const auto forward_node = entry.forward_ebg_node;
        const auto backward_node = entry.backward_ebg_node;

        if (node_based_partition_ids[u] == node_based_partition_ids[v])
        {
            // Can use partition_ids[u/v] as partition for edge based graph `node_id`
            edge_based_partition_ids[forward_node] = node_based_partition_ids[u];
            if (backward_node != SPECIAL_NODEID)
                edge_based_partition_ids[backward_node] = node_based_partition_ids[u];
        }
    }

    // Heuristic: Pick the bisection ID of u or v depending on how many border vertexes
    // it would induce, given the current (partital) assignment
    for (const auto &entry : mapping)
    {
        const auto u = entry.u;
        const auto v = entry.v;
        const auto forward_node = entry.forward_ebg_node;
        const auto backward_node = entry.backward_ebg_node;

        if (edge_based_partition_ids[forward_node] == SPECIAL_NODEID)
        {
            BOOST_ASSERT(backward_node == SPECIAL_NODEID ||
                         edge_based_partition_ids[backward_node] == SPECIAL_NODEID);
            // Border nodes u,v - need to be resolved.
            // FIXME: just pick one side for now. See #3205.

            std::size_t u_border_edges = 0;
            std::size_t v_border_edges = 0;

            const auto count_border_egdes = [&](NodeID ebg_node) {
                for (auto edge : edge_based_graph->GetAdjacentEdgeRange(ebg_node))
                {
                    auto target = edge_based_graph->GetTarget(edge);
                    if (edge_based_partition_ids[target] != node_based_partition_ids[u])
                        u_border_edges++;
                    if (edge_based_partition_ids[target] != node_based_partition_ids[v])
                        v_border_edges++;
                    // note: the target node can be neither in u or v's partition
                }
            };

            count_border_egdes(forward_node);
            if (backward_node != SPECIAL_NODEID)
                count_border_egdes(backward_node);

            bool use_u = u_border_edges < v_border_edges;

            // Use partition that introduce less cross cell connections
            edge_based_partition_ids[forward_node] = node_based_partition_ids[use_u ? u : v];
            if (backward_node != SPECIAL_NODEID)
                edge_based_partition_ids[backward_node] = node_based_partition_ids[use_u ? u : v];
        }

        // after this we have resolved all nodes
        BOOST_ASSERT(edge_based_partition_ids[forward_node] != SPECIAL_NODEID);
        BOOST_ASSERT(backward_node == SPECIAL_NODEID ||
                     edge_based_partition_ids[backward_node] != SPECIAL_NODEID);
    }

    std::vector<Partition> partitions;
    std::vector<std::uint32_t> level_to_num_cells;
    std::tie(partitions, level_to_num_cells) =
        bisectionToPartition(edge_based_partition_ids,
                             {config.minimum_cell_size,
                              config.minimum_cell_size * 32,
                              config.minimum_cell_size * 32 * 16,
                              config.minimum_cell_size * 32 * 16 * 32});

    auto num_fixed_unconnected = 0;
    auto num_unconnected = 0;
    for (int level_index = partitions.size()-1; level_index >= 0; level_index--)
    {
        std::vector<std::tuple<CellID, NodeID>> forward_witnesses;
        std::vector<std::tuple<CellID, NodeID>> backward_witnesses;
        for (const auto &entry : mapping)
        {
            forward_witnesses.clear();
            backward_witnesses.clear();

            bool forward_is_source = false;
            bool forward_is_target = false;
            bool backward_is_source = false;
            bool backward_is_target = false;

            const auto find_witnesses =
                [&](const NodeID node, bool &is_source, bool &is_target, auto &witnesses) {
                    const auto cell_id = partitions[level_index][node];
                    for (auto edge : edge_based_graph->GetAdjacentEdgeRange(node))
                    {
                        const auto data = edge_based_graph->GetEdgeData(edge);
                        const auto target = edge_based_graph->GetTarget(edge);
                        const auto target_cell_id = partitions[level_index][target];
                        if (target_cell_id == cell_id)
                        {
                            is_source |= data.forward;
                            is_target |= data.backward;
                        }
                        else
                        {
                            witnesses.push_back(std::make_tuple(target_cell_id, target));
                        }
                    }
                };

            find_witnesses(
                entry.forward_ebg_node, forward_is_source, forward_is_target, forward_witnesses);
            const auto forward_unconnected =
                forward_witnesses.size() > 0 && !forward_is_source && !forward_is_target;

            if (entry.backward_ebg_node != SPECIAL_NODEID)
            {
                find_witnesses(entry.backward_ebg_node,
                               backward_is_source,
                               backward_is_target,
                               backward_witnesses);
            }
            const auto backward_unconnected =
                backward_witnesses.size() > 0 && !backward_is_source && !backward_is_target;

            if (forward_unconnected)
                num_unconnected++;

            if (backward_unconnected)
                num_unconnected++;

            if (forward_unconnected && entry.backward_ebg_node == SPECIAL_NODEID)
            {
                num_fixed_unconnected++;

                // FIXME actually fix
            }
            else if (forward_unconnected && backward_unconnected && entry.backward_ebg_node != SPECIAL_NODEID)
            {

                // both nodes are unconnected to the cell in which they are contained
                if (forward_unconnected && backward_unconnected)
                {
                    // we need to find a witness that puts both forward and
                    // backward node in the same cell
                    std::sort(forward_witnesses.begin(), forward_witnesses.end());
                    std::sort(backward_witnesses.begin(), backward_witnesses.end());

                    decltype(forward_witnesses) merged_witnesses;
                    std::set_intersection(forward_witnesses.begin(),
                                          forward_witnesses.end(),
                                          backward_witnesses.begin(),
                                          backward_witnesses.end(),
                                          std::back_inserter(merged_witnesses),
                                          [](const auto &lhs, const auto &rhs) {
                                              return std::get<0>(lhs) < std::get<0>(rhs);
                                          });

                    if (merged_witnesses.size() > 0)
                    {
                        num_fixed_unconnected += 2;
                    }
                }
                // FIXME actually fix
            }
        }
    }

    util::Log() << "Fixed " << num_fixed_unconnected << " out of " << num_unconnected
                << " unconnected nodes";

    util::Log() << "Edge-based-graph annotation:";
    for (std::size_t level = 0; level < level_to_num_cells.size(); ++level)
    {
        util::Log() << "  level " << level + 1 << " #cells " << level_to_num_cells[level]
                    << " bit size " << std::ceil(std::log2(level_to_num_cells[level] + 1));
    }

    TIMER_START(packed_mlp);
    MultiLevelPartition mlp{partitions, level_to_num_cells};
    TIMER_STOP(packed_mlp);
    util::Log() << "MultiLevelPartition constructed in " << TIMER_SEC(packed_mlp) << " seconds";

    TIMER_START(cell_storage);
    CellStorage storage(mlp, *edge_based_graph);
    TIMER_STOP(cell_storage);
    util::Log() << "CellStorage constructed in " << TIMER_SEC(cell_storage) << " seconds";

    TIMER_START(writing_mld_data);
    io::write(config.mld_partition_path, mlp);
    io::write(config.mld_storage_path, storage);
    TIMER_STOP(writing_mld_data);
    util::Log() << "MLD data writing took " << TIMER_SEC(writing_mld_data) << " seconds";

    return 0;
}

} // namespace partition
} // namespace osrm
