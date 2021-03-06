#include "partition/partitioner.hpp"
#include "partition/bisection_graph.hpp"
#include "partition/bisection_to_partition.hpp"
#include "partition/cell_storage.hpp"
#include "partition/compressed_node_based_graph_reader.hpp"
#include "partition/edge_based_graph_reader.hpp"
#include "partition/io.hpp"
#include "partition/multi_level_partition.hpp"
#include "partition/node_based_graph_to_edge_based_graph_mapping_reader.hpp"
#include "partition/recursive_bisection.hpp"

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

    auto mapping = LoadNodeBasedGraphToEdgeBasedGraphMapping(config.cnbg_ebg_mapping_path.string());
    util::Log() << "Loaded node based graph to edge based graph mapping";

    auto edge_based_graph = LoadEdgeBasedGraph(config.edge_based_graph_path.string());
    util::Log() << "Loaded edge based graph for mapping partition ids: "
                << edge_based_graph->GetNumberOfEdges() << " edges, "
                << edge_based_graph->GetNumberOfNodes() << " nodes";

    // TODO: node based graph to edge based graph partition id mapping should be done split off.

    // Partition ids, keyed by node based graph nodes
    const auto &node_based_partition_ids = recursive_bisection.BisectionIDs();

    // Partition ids, keyed by edge based graph nodes
    std::vector<NodeID> edge_based_partition_ids(edge_based_graph->GetNumberOfNodes());

    // Extract edge based border nodes, based on node based partition and mapping.
    for (const auto node : util::irange(0u, edge_based_graph->GetNumberOfNodes()))
    {
        const auto node_based_nodes = mapping.Lookup(node);

        const auto u = node_based_nodes.u;
        const auto v = node_based_nodes.v;

        if (node_based_partition_ids[u] == node_based_partition_ids[v])
        {
            // Can use partition_ids[u/v] as partition for edge based graph `node_id`
            edge_based_partition_ids[node] = node_based_partition_ids[u];

            auto edges = edge_based_graph->GetAdjacentEdgeRange(node);
            if (edges.size() == 1)
            { // Check the edge case with one adjacent edge-based backward edge
                auto edge = edges.front();
                auto other = edge_based_graph->GetTarget(edge);
                auto &data = edge_based_graph->GetEdgeData(edge);
                auto other_node_based_nodes = mapping.Lookup(other);
                if (data.backward &&
                    node_based_partition_ids[other_node_based_nodes.u] !=
                        node_based_partition_ids[u])
                { // use id of other node if the edge [other_u, other_v] -> [u,v] is a single edge
                    // and nodes other_[u,v] are in  different node-based partitions
                    edge_based_partition_ids[node] =
                        node_based_partition_ids[other_node_based_nodes.u];
                }
            }
        }
        else
        {
            // Border nodes u,v - need to be resolved.
            // FIXME: just pick one side for now. See #3205.

            bool use_u = false;
            for (auto edge : edge_based_graph->GetAdjacentEdgeRange(node))
            {
                auto other = edge_based_graph->GetTarget(edge);
                auto &data = edge_based_graph->GetEdgeData(edge);
                auto other_node_based_nodes = mapping.Lookup(other);

                if (data.backward)
                { // can use id of u if [other_u, other_v] -> [u,v] is in the same partition as u
                    BOOST_ASSERT(u == other_node_based_nodes.v);
                    use_u |= node_based_partition_ids[u] ==
                             node_based_partition_ids[other_node_based_nodes.u];
                }
            }

            // Use partition that introduce less cross cell connections
            edge_based_partition_ids[node] = node_based_partition_ids[use_u ? u : v];
        }
    }

    std::vector<Partition> partitions;
    std::vector<std::uint32_t> level_to_num_cells;
    std::tie(partitions, level_to_num_cells) =
        bisectionToPartition(edge_based_partition_ids,
                             {config.minimum_cell_size,
                              config.minimum_cell_size * 32,
                              config.minimum_cell_size * 32 * 16,
                              config.minimum_cell_size * 32 * 16 * 32});

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
