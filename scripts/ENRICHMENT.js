/************************************************************
  ENRICHMENT.js — Sobel-Based Ore Deposit Enrichment Analysis
  Author: Dilyara Nazyrova, 2025
  Project: Sobel-Morphostructure-KZ
  Purpose: Quantitative evaluation of mineral deposit enrichment
           within Sobel-derived high-gradient morphostructural zones.
************************************************************/

// ====================== 0. INPUT DATA ======================

// AOI and DEM
var AOI = ee.FeatureCollection('projects/ee-euniversity2024/assets/AOI');
var DEM = ee.Image('MERIT/DEM/v1_0_3').select('dem').clip(AOI);

// Deposit dataset
var DEPOSITS = ee.FeatureCollection('projects/ee-euniversity2024/assets/deposits_720');

// ====================== 1. PARAMETERS ======================
var SCALE = 90;
var MAXPX = 1e13;
var TOP_PCTL = 85; // top 15%
var REGION = AOI.geometry();

// ====================== 2. SOBEL FILTER ======================

// Compute gradients
var kx = ee.Kernel.fixed(3,3,[[-1,0,1],[-2,0,2],[-1,0,1]]);
var ky = ee.Kernel.fixed(3,3,[[-1,-2,-1],[0,0,0],[1,2,1]]);
var gx = DEM.convolve(kx);
var gy = DEM.convolve(ky);
var sobelMag = gx.pow(2).add(gy.pow(2)).sqrt().rename('sobelMag');

// Threshold and mask
var sobelThr = ee.Number(
  sobelMag.reduceRegion({
    reducer: ee.Reducer.percentile([TOP_PCTL]),
    geometry: REGION,
    scale: SCALE, bestEffort:true, tileScale:4, maxPixels:MAXPX
  }).get('sobelMag')
);
print('Sobel threshold (top 15%):', sobelThr);

var sobelMask = sobelMag.gte(sobelThr).selfMask();
Map.centerObject(AOI, 8);
Map.addLayer(sobelMask, {palette:['#ffff00']}, 'Sobel top mask');

// ====================== 3. TOPOGRAPHIC ZONE SHARE ======================

var areaImg = ee.Image.pixelArea();
var A_total = ee.Number(areaImg.reduceRegion({
  reducer: ee.Reducer.sum(), geometry: REGION, scale:SCALE, bestEffort:true, tileScale:4, maxPixels:MAXPX
}).get('area')).divide(1e6);

var A_top = ee.Number(areaImg.updateMask(sobelMask).reduceRegion({
  reducer: ee.Reducer.sum(), geometry: REGION, scale:SCALE, bestEffort:true, tileScale:4, maxPixels:MAXPX
}).get('area')).divide(1e6);

var top_share_pct = A_top.divide(A_total).multiply(100);
print('Top-zone area share (%):', top_share_pct);

// ====================== 4. DEPOSIT FLAGS ======================

var depositsFlagged = DEPOSITS.map(function(f){
  var v = sobelMask.sample({
    region: f.geometry(), scale: SCALE, numPixels:1, geometries:false
  }).first();
  var inside = ee.Algorithms.If(v, 1, 0);
  return f.set('inside_top', inside);
});

var dep_total = depositsFlagged.size();
var dep_inside = depositsFlagged.aggregate_sum('inside_top');
var dep_inside_pct = ee.Number(dep_inside).divide(dep_total).multiply(100);

print('Total deposits:', dep_total);
print('Deposits inside top zones (%):', dep_inside_pct);

// ====================== 5. ENRICHMENT INDEX ======================

var enrichment = ee.Number(dep_inside_pct).divide(top_share_pct);
print('Enrichment index (E):', enrichment);

// ====================== 6. DISTANCE TO SOBEL ZONES ======================

var background = ee.Image(1).subtract(sobelMask.unmask(0)).rename('bg');
var dist_m = background.fastDistanceTransform(30).sqrt().multiply(SCALE).rename('dist_m');

var depositsWithDist = dist_m.reduceRegions({
  collection: depositsFlagged,
  reducer: ee.Reducer.first(),
  scale: SCALE, tileScale: 4
}).map(function(f){
  var d = ee.Number(f.get('first'));
  return ee.Feature(f.geometry(), {
    inside_top: f.get('inside_top'),
    dist_m: d
  });
});

var distStats = depositsWithDist.reduceColumns({
  reducer: ee.Reducer.median().combine({
    reducer2: ee.Reducer.percentile([25, 75]), sharedInputs: true
  }),
  selectors: ['dist_m']
});
print('Distance to Sobel-top zone (m): median / p25 / p75 =', distStats);

// ====================== 7. EXPORTS ======================

// Summary table
Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, {
      top_share_pct: top_share_pct,
      dep_total: dep_total,
      dep_inside_pct: dep_inside_pct,
      enrichment: enrichment
    })
  ]),
  description: 'Sobel_enrichment_summary',
  fileNamePrefix: 'Sobel_enrichment_summary',
  fileFormat: 'CSV'
});

// Individual deposit distances
Export.table.toDrive({
  collection: depositsWithDist,
  description: 'Sobel_deposits_distance',
  fileNamePrefix: 'Sobel_deposits_distance',
  fileFormat: 'CSV'
});

// ====================== 8. DISTANCE TO SOBEL ZONE EDGES ======================

var scl = SCALE;

// Distance from outside pixels to nearest inside (top zone)
var outside = ee.Image(1).subtract(sobelMask.unmask(0));
var distOutside = outside.fastDistanceTransform(1).sqrt().multiply(scl);

// Distance from inside pixels to nearest outside
var inside = sobelMask.unmask(0);
var distInside = inside.fastDistanceTransform(1).sqrt().multiply(scl);

// Distance to the nearest edge = max(distInside, distOutside)
var distEdge_m = ee.Image.cat(distInside, distOutside).reduce(ee.Reducer.max())
  .rename('dist_edge_m');

// Assign distances to deposit points
var depositsWithEdgeDist = distEdge_m.reduceRegions({
  collection: depositsFlagged,
  reducer: ee.Reducer.first(),
  scale: scl, tileScale: 4
}).map(function(f){
  return ee.Feature(f.geometry(), {
    inside_top: f.get('inside_top'),
    dist_edge_m: ee.Number(f.get('first'))
  });
});

// Edge distance statistics
var edgeStats_all = depositsWithEdgeDist.reduceColumns({
  reducer: ee.Reducer.median().combine({
    reducer2: ee.Reducer.percentile([25, 75]), sharedInputs: true
  }),
  selectors: ['dist_edge_m']
});
print('Distance to Sobel edge (ALL): median / p25 / p75 (m) =', edgeStats_all);

// Filter out zero distances (points exactly on the edge)
var depositsNonZero = depositsWithEdgeDist.filter(ee.Filter.gt('dist_edge_m', 0));
var edgeStats_nz = depositsNonZero.reduceColumns({
  reducer: ee.Reducer.median().combine({
    reducer2: ee.Reducer.percentile([25, 75]), sharedInputs: true
  }),
  selectors: ['dist_edge_m']
});
print('Distance to Sobel edge (>0 only): median / p25 / p75 (m) =', edgeStats_nz);

// Export edge distances
Export.table.toDrive({
  collection: depositsWithEdgeDist,
  description: 'Sobel_deposits_edge_distance',
  fileNamePrefix: 'Sobel_deposits_edge_distance',
  fileFormat: 'CSV'
});

// ====================== 9. ENRICHMENT TRADE-OFF ANALYSIS ======================

var PCTL_LIST = [70, 75, 80, 85, 90, 95];

var results = ee.FeatureCollection(PCTL_LIST.map(function(pct){
  var thr = ee.Number(
    sobelMag.reduceRegion({
      reducer: ee.Reducer.percentile([pct]),
      geometry: REGION,
      scale: SCALE, bestEffort: true, tileScale: 4, maxPixels: MAXPX
    }).get('sobelMag')
  );
  var mask = sobelMag.gte(thr).selfMask();

  var A_top = ee.Number(areaImg.updateMask(mask).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: REGION, scale:SCALE, bestEffort:true, tileScale:4, maxPixels:MAXPX
  }).get('area')).divide(1e6);
  var A_share = A_top.divide(A_total).multiply(100);

  var depIn = DEPOSITS.map(function(f){
    var v = mask.sample({
      region: f.geometry(), scale: SCALE, numPixels:1, geometries:false
    }).first();
    return f.set('inside', ee.Algorithms.If(v, 1, 0));
  });
  var depInPct = ee.Number(depIn.aggregate_sum('inside')).divide(depIn.size()).multiply(100);
  var enr = depInPct.divide(A_share);

  return ee.Feature(null, {
    percentile: pct,
    sobel_thr: thr,
    coverage_pct: A_share,
    inside_pct: depInPct,
    enrichment: enr
  });
}));

print('SOBEL Enrichment–Coverage trade-off:', results);

print(ui.Chart.feature.byFeature(results, 'coverage_pct', 'enrichment')
  .setOptions({
    title: 'Sobel Enrichment vs Coverage',
    hAxis: {title: 'Top-zone coverage (%)'},
    vAxis: {title: 'Enrichment (Deposits / Area)'},
    pointSize: 6,
    lineWidth: 2,
    colors: ['#d95f02']
  })
);

Export.table.toDrive({
  collection: results,
  description: 'Sobel_Enrichment_Tradeoff',
  fileNamePrefix: 'Sobel_Enrichment_Tradeoff',
  fileFormat: 'CSV'
});

// ====================== 10. QUALITY CONTROL ======================

print('Sample of first 10 edge-distance records:');
print(depositsWithEdgeDist.limit(10));

var distSummary = depositsWithEdgeDist.reduceColumns({
  reducer: ee.Reducer.minMax()
    .combine({reducer2: ee.Reducer.count(), sharedInputs: true})
    .combine({reducer2: ee.Reducer.percentile([25, 50, 75]), sharedInputs: true}),
  selectors: ['dist_edge_m']
});
print('Statistics of dist_edge_m (min/max/count/p25/p50/p75):', distSummary);

var zeroCount = depositsWithEdgeDist.filter(ee.Filter.eq('dist_edge_m', 0)).size();
print('Number of deposits with dist_edge_m = 0:', zeroCount);

var uniqueVals = depositsWithEdgeDist.aggregate_array('dist_edge_m').distinct().size();
print('Number of unique dist_edge_m values:', uniqueVals);

print('Sample of first 10 Sobel-zone distances (dist_m):');
print(depositsWithDist.limit(10));

var distSummary2 = depositsWithDist.reduceColumns({
  reducer: ee.Reducer.minMax()
    .combine({reducer2: ee.Reducer.count(), sharedInputs: true})
    .combine({reducer2: ee.Reducer.percentile([25, 50, 75]), sharedInputs: true}),
  selectors: ['dist_m']
});
print('Statistics of dist_m (min/max/count/p25/p50/p75):', distSummary2);
