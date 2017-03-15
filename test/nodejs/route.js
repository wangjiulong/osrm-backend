var OSRM = require('../../');
var test = require('tape');
var berlin_path = require('./osrm-data-path').data_path;

test('route: routes Berlin', function(assert) {
    assert.plan(5);
    var osrm = new OSRM(berlin_path);
    osrm.route({coordinates: [[13.43864,52.51993],[13.415852,52.513191]]}, function(err, route) {
        assert.ifError(err);
        assert.ok(route.waypoints);
        assert.ok(route.routes);
        assert.ok(route.routes.length);
        assert.ok(route.routes[0].geometry);
    });
});

test('route: throws with too few or invalid args', function(assert) {
    assert.plan(3);
    var osrm = new OSRM(berlin_path);
    assert.throws(function() { osrm.route({coordinates: [[13.43864,52.51993],[13.415852,52.513191]]}) },
        /Two arguments required/);
    assert.throws(function() { osrm.route(null, function(err, route) {}) },
        /First arg must be an object/);
    assert.throws(function() { osrm.route({coordinates: [[13.43864,52.51993],[13.415852,52.513191]]}, true)},
        /last argument must be a callback function/);
});

test('route: provides no alternatives by default, but when requested', function(assert) {
    assert.plan(6);
    var osrm = new OSRM(berlin_path);
    var options = {coordinates: [[13.302383,52.490516], [13.418427,52.522070]]};

    osrm.route(options, function(err, route) {
        assert.ifError(err);
        assert.ok(route.routes);
        assert.equal(route.routes.length, 1);
    });
    options.alternatives = true;
    osrm.route(options, function(err, route) {
        assert.ifError(err);
        assert.ok(route.routes);
        assert.equal(route.routes.length, 2);
    });
});

test('route: throws with bad params', function(assert) {
    assert.plan(11);
    var osrm = new OSRM(berlin_path);
    assert.throws(function () { osrm.route({coordinates: []}, function(err) {}) });
    assert.throws(function() { osrm.route({}, function(err, route) {}) },
        /Must provide a coordinates property/);
    assert.throws(function() { osrm.route({coordinates: null}, function(err, route) {}) },
        /Coordinates must be an array of \(lon\/lat\) pairs/);
    assert.throws(function() { osrm.route({coordinates: [13.438640, 52.519930]}, function(err, route) {}) },
        /Coordinates must be an array of \(lon\/lat\) pairs/);
    assert.throws(function() { osrm.route({coordinates: [[true, '52.519930'], [13.438640, 52.519930]]}, function(err, route) {}) },
        /Each member of a coordinate pair must be a number/);
    assert.throws(function() { osrm.route({coordinates: [[213.43864,252.51993],[413.415852,552.513191]]}, function(err, route) {}) },
        /Lng\/Lat coordinates must be within world bounds \(-180 < lng < 180, -90 < lat < 90\)/);
    assert.throws(function() { osrm.route({coordinates: [[13.438640], [52.519930]]}, function(err, route) {}) },
        /Coordinates must be an array of \(lon\/lat\) pairs/);
    assert.throws(function() { osrm.route({coordinates: [[13.43864,52.51993],[13.415852,52.513191]], hints: null}, function(err, route) {}) },
        /Hints must be an array of strings\/null/);
    assert.throws(function() { osrm.route({coordinates: [[13.43864,52.51993],[13.415852,52.513191]], steps: null}, function(err, route) {}) });
    assert.throws(function() { osrm.route({coordinates: [[13.43864,52.51993],[13.415852,52.513191]], annotations: null}, function(err, route) {}) });
    var options = {
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        alternateRoute: false,
        hints: [13.438640, 52.519930]
    };
    assert.throws(function() { osrm.route(options, function(err, route) {}) },
        /Hint must be null or string/);
});

test('route: routes Berlin using shared memory', function(assert) {
    assert.plan(2);
    var osrm = new OSRM();
    osrm.route({coordinates: [[13.43864,52.51993],[13.415852,52.513191]]}, function(err, route) {
        assert.ifError(err);
        assert.ok(Array.isArray(route.routes));
    });
});

test('route: routes Berlin with geometry compression', function(assert) {
    assert.plan(2);
    var osrm = new OSRM(berlin_path);
    var options = {
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
    };
    osrm.route(options, function(err, route) {
        assert.ifError(err);
        assert.equal('string', typeof route.routes[0].geometry);
    });
});

test('route: routes Berlin without geometry compression', function(assert) {
    assert.plan(4);
    var osrm = new OSRM(berlin_path);
    var options = {
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        geometries: 'geojson'
    };
    osrm.route(options, function(err, route) {
        assert.ifError(err);
        assert.ok(Array.isArray(route.routes));
        assert.ok(Array.isArray(route.routes[0].geometry.coordinates));
        assert.equal(route.routes[0].geometry.type, 'LineString');
    });
});

test('Test polyline6 geometries option', function(assert) {
    assert.plan(6);
    var osrm = new OSRM(berlin_path);
    var options = {
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        continue_straight: false,
        overview: 'false',
        geometries: 'polyline6',
        steps: true
    };
    osrm.route(options, function(err, first) {
        assert.ifError(err);
        assert.ok(first.routes);
        assert.equal(first.routes.length, 1);
        assert.notOk(first.routes[0].geometry);
        assert.ok(first.routes[0].legs[0]);
        assert.equal(typeof first.routes[0].legs[0].steps[0].geometry, 'string');
    });
});

test('route: routes Berlin with speed annotations options', function(assert) {
    assert.plan(17);
    var osrm = new OSRM(berlin_path);
    var options = {
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        continue_straight: false,
        overview: 'false',
        geometries: 'polyline',
        steps: true,
        annotations: ['speed']
    };
    osrm.route(options, function(err, first) {
        assert.ifError(err);
        assert.ok(first.routes);
        assert.ok(first.routes[0].legs.every(function(l) { return Array.isArray(l.steps) && l.steps.length > 0; }));
        assert.equal(first.routes.length, 1);
        assert.notOk(first.routes[0].geometry);
        assert.ok(first.routes[0].legs[0]);
        assert.ok(first.routes[0].legs.every(l => { return l.steps.length > 0; }), 'every leg has steps');
        assert.ok(first.routes[0].legs.every(l => { return l.annotation;}), 'every leg has annotations');
        assert.ok(first.routes[0].legs.every(l => { return l.annotation.speed;}), 'every leg has annotations for speed');
        assert.notOk(first.routes[0].legs.every(l => { return l.annotation.weight; }), 'has no annotations for weight')
        assert.notOk(first.routes[0].legs.every(l => { return l.annotation.datasources; }), 'has no annotations for datasources')
        assert.notOk(first.routes[0].legs.every(l => { return l.annotation.duration; }), 'has no annotations for duration')
        assert.notOk(first.routes[0].legs.every(l => { return l.annotation.distance; }), 'has no annotations for distance')
        assert.notOk(first.routes[0].legs.every(l => { return l.annotation.nodes; }), 'has no annotations for nodes')

        options.overview = 'full';
        osrm.route(options, function(err, full) {
            assert.ifError(err);
            options.overview = 'simplified';
            osrm.route(options, function(err, simplified) {
                assert.ifError(err);
                assert.notEqual(full.routes[0].geometry, simplified.routes[0].geometry);
            });
        });
    });
});

test('route: routes Berlin with several (duration, distance, nodes) annotations options', function(assert) {
    assert.plan(17);
    var osrm = new OSRM(berlin_path);
    var options = {
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        continue_straight: false,
        overview: 'false',
        geometries: 'polyline',
        steps: true,
        annotations: ['duration','distance','nodes']
    };
    osrm.route(options, function(err, first) {
        assert.ifError(err);
        assert.ok(first.routes);
        assert.ok(first.routes[0].legs.every(function(l) { return Array.isArray(l.steps) && l.steps.length > 0; }));
        assert.equal(first.routes.length, 1);
        assert.notOk(first.routes[0].geometry);
        assert.ok(first.routes[0].legs[0]);
        assert.ok(first.routes[0].legs.every(l => { return l.steps.length > 0; }), 'every leg has steps');
        assert.ok(first.routes[0].legs.every(l => { return l.annotation;}), 'every leg has annotations');
        assert.ok(first.routes[0].legs.every(l => { return l.annotation.distance;}), 'every leg has annotations for distance');
        assert.ok(first.routes[0].legs.every(l => { return l.annotation.duration;}), 'every leg has annotations for durations');
        assert.ok(first.routes[0].legs.every(l => { return l.annotation.nodes;}), 'every leg has annotations for nodes');
        assert.notOk(first.routes[0].legs.every(l => { return l.annotation.weight; }), 'has no annotations for weight')
        assert.notOk(first.routes[0].legs.every(l => { return l.annotation.datasources; }), 'has no annotations for datasources')
        assert.notOk(first.routes[0].legs.every(l => { return l.annotation.speed; }), 'has no annotations for speed')

        options.overview = 'full';
        osrm.route(options, function(err, full) {
            assert.ifError(err);
            options.overview = 'simplified';
            osrm.route(options, function(err, simplified) {
                assert.ifError(err);
                assert.notEqual(full.routes[0].geometry, simplified.routes[0].geometry);
            });
        });
    });
});

test('route: routes Berlin with options', function(assert) {
    assert.plan(11);
    var osrm = new OSRM(berlin_path);
    var options = {
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        continue_straight: false,
        overview: 'false',
        geometries: 'polyline',
        steps: true,
        annotations: true
    };
    osrm.route(options, function(err, first) {
        assert.ifError(err);
        assert.ok(first.routes);
        assert.ok(first.routes[0].legs.every(function(l) { return Array.isArray(l.steps) && l.steps.length > 0; }));
        assert.equal(first.routes.length, 1);
        assert.notOk(first.routes[0].geometry);
        assert.ok(first.routes[0].legs[0]);
        assert.ok(first.routes[0].legs.every(l => { return l.steps.length > 0; }), 'every leg has steps');
        assert.ok(first.routes[0].legs.every(l => { return l.annotation;}), 'every leg has annotations');

        options.overview = 'full';
        osrm.route(options, function(err, full) {
            assert.ifError(err);
            options.overview = 'simplified';
            osrm.route(options, function(err, simplified) {
                assert.ifError(err);
                assert.notEqual(full.routes[0].geometry, simplified.routes[0].geometry);
            });
        });
    });
});

test('route: routes Berlin with options', function(assert) {
    assert.plan(11);
    var osrm = new OSRM(berlin_path);
    var options = {
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        continue_straight: false,
        overview: 'false',
        geometries: 'polyline',
        steps: true,
        annotations: true
    };
    osrm.route(options, function(err, first) {
        assert.ifError(err);
        assert.ok(first.routes);
        assert.ok(first.routes[0].legs.every(function(l) { return Array.isArray(l.steps) && l.steps.length > 0; }));
        assert.equal(first.routes.length, 1);
        assert.notOk(first.routes[0].geometry);
        assert.ok(first.routes[0].legs[0]);
        assert.ok(first.routes[0].legs.every(l => { return l.steps.length > 0; }), 'every leg has steps');
        assert.ok(first.routes[0].legs.every(l => { return l.annotation;}), 'every leg has annotations');

        options.overview = 'full';
        osrm.route(options, function(err, full) {
            assert.ifError(err);
            options.overview = 'simplified';
            osrm.route(options, function(err, simplified) {
                assert.ifError(err);
                assert.notEqual(full.routes[0].geometry, simplified.routes[0].geometry);
            });
        });
    });
});

test('route: invalid route options', function(assert) {
    assert.plan(8);
    var osrm = new OSRM(berlin_path);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        continue_straight: []
    }, function(err, route) {}); },
        /must be boolean/);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        alternatives: []
    }, function(err, route) {}); },
        /must be boolean/);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        geometries: true
    }, function(err, route) {}); },
        /Geometries must be a string: \[polyline, polyline6, geojson\]/);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        overview: false
    }, function(err, route) {}); },
        /Overview must be a string: \[simplified, full, false\]/);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        overview: false
    }, function(err, route) {}); },
        /Overview must be a string: \[simplified, full, false\]/);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        overview: 'maybe'
    }, function(err, route) {}); },
        /'overview' param must be one of \[simplified, full, false\]/);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        geometries: 'maybe'
    }, function(err, route) {}); },
        /'geometries' param must be one of \[polyline, polyline6, geojson\]/);
    assert.throws(function() { osrm.route({
        coordinates: [[NaN, -NaN],[Infinity, -Infinity]]
    }, function(err, route) {}); },
        /Lng\/Lat coordinates must be valid numbers/);
});

test('route: integer bearing values no longer supported', function(assert) {
    assert.plan(1);
    var osrm = new OSRM(berlin_path);
    var options = {
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        bearings: [200, 250],
    };
    assert.throws(function() { osrm.route(options, function(err, route) {}); },
        /Bearing must be an array of \[bearing, range\] or null/);
});

test('route: valid bearing values', function(assert) {
    assert.plan(4);
    var osrm = new OSRM(berlin_path);
    var options = {
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        bearings: [[200, 180], [250, 180]],
    };
    osrm.route(options, function(err, route) {
        assert.ifError(err);
        assert.ok(route.routes[0]);
    });
    options.bearings = [null, [360, 180]];
    osrm.route(options, function(err, route) {
        assert.ifError(err);
        assert.ok(route.routes[0]);
    });
});

test('route: invalid bearing values', function(assert) {
    assert.plan(6);
    var osrm = new OSRM(berlin_path);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        bearings: [[400, 180], [-250, 180]],
    }, function(err, route) {}) },
        /Bearing values need to be in range 0..360, 0..180/);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        bearings: [[200], [250, 180]],
    }, function(err, route) {}) },
        /Bearing must be an array of/);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        bearings: [[400, 109], [100, 720]],
    }, function(err, route) {}) },
        /Bearing values need to be in range 0..360, 0..180/);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        bearings: 400,
    }, function(err, route) {}) },
        /Bearings must be an array of arrays of numbers/);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        bearings: [[100, 100]],
    }, function(err, route) {}) },
        /Bearings array must have the same length as coordinates array/);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        bearings: [Infinity, Infinity],
    }, function(err, route) {}) },
        /Bearing must be an array of \[bearing, range\] or null/);
});

test('route: routes Berlin with hints', function(assert) {
    assert.plan(5);
    var osrm = new OSRM(berlin_path);
    var options = {
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]]
    };
    osrm.route(options, function(err, first) {
        assert.ifError(err);
        assert.ok(first.waypoints);
        var hints = first.waypoints.map(function(wp) { return wp.hint; });
        assert.ok(hints.every(function(h) { return typeof h === 'string'; }));

        options.hints = hints;

        osrm.route(options, function(err, second) {
            assert.ifError(err);
            assert.deepEqual(first, second);
        });
    });
});

test('route: routes Berlin with null hints', function(assert) {
    assert.plan(1);
    var osrm = new OSRM(berlin_path);
    var options = {
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        hints: [null, null]
    };
    osrm.route(options, function(err, route) {
        assert.ifError(err);
    });
});

test('route: throws on bad hints', function(assert) {
    assert.plan(2);
    var osrm = new OSRM(berlin_path);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        hints: ['', '']
    }, function(err, route) {})}, /Hint cannot be an empty string/);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        hints: [null]
    }, function(err, route) {})}, /Hints array must have the same length as coordinates array/);
});

test('route: routes Berlin with valid radius values', function(assert) {
    assert.plan(3);
    var osrm = new OSRM(berlin_path);
    var options = {
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        radiuses: [100, 100]
    };
    osrm.route(options, function(err, route) {
        assert.ifError(err);
    });
    options.radiuses = [null, null];
    osrm.route(options, function(err, route) {
        assert.ifError(err);
    });
    options.radiuses = [100, null];
    osrm.route(options, function(err, route) {
        assert.ifError(err);
    });
});

test('route: throws on bad radiuses', function(assert) {
    assert.plan(3);
    var osrm = new OSRM(berlin_path);
    var options = {
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        radiuses: [10, 10]
    };
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        radiuses: 10
    }, function(err, route) {}) },
        /Radiuses must be an array of non-negative doubles or null/);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        radiuses: ['magic', 'numbers']
    }, function(err, route) {}) },
        /Radius must be non-negative double or null/);
    assert.throws(function() { osrm.route({
        coordinates: [[13.43864,52.51993],[13.415852,52.513191]],
        radiuses: [10]
    }, function(err, route) {}) },
        /Radiuses array must have the same length as coordinates array/);
});
