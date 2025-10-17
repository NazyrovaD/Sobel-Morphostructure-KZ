/********* SOBEL — Final version (top-15%, consistent with MDPI article) *********/

// 0) Inputs
var AOI    = ee.FeatureCollection('projects/ee-euniversity2024/assets/AOI');
var DEM    = ee.Image('MERIT/DEM/v1_0_3').select('dem').clip(AOI);
var FAULTS = ee.FeatureCollection('projects/ee-euniversity2024/assets/Razlom');
Map.centerObject(AOI, 8);

// 1) Parameters
var DEM_PROJ = DEM.projection();
var SCALE    = DEM_PROJ.nominalScale();   // native MERIT DEM ≈90 m
var CRS      = DEM_PROJ.crs();
var MAXPX    = 1e13;

// threshold = top-15 % (≈85th percentile)
var TOP_SHARE = 15;                        // % of strongest gradients
var TOP_PCTL  = ee.Number(100).subtract(TOP_SHARE);  // 85th percentile
var BUF_LIST_M = [1000, 3000, 6000];
var BIN_SIZE  = 15;                         // histogram step (°)
var BINS      = ee.List.sequence(0, 180 - BIN_SIZE, BIN_SIZE);

// 2) Sobel filter (X, Y derivatives and magnitude)
var kx = ee.Kernel.fixed(3,3,[[-1,0,1],[-2,0,2],[-1,0,1]]);
var ky = ee.Kernel.fixed(3,3,[[-1,-2,-1],[0,0,0],[1,2,1]]);
var sobelX = DEM.convolve(kx);
var sobelY = DEM.convolve(ky);
var sobelMag = sobelX.pow(2).add(sobelY.pow(2)).sqrt().rename('sobelMag');

Map.addLayer(sobelMag, {min:0, max:100, palette:['black','white']}, 'Sobel magnitude');
Map.addLayer(FAULTS.style({color:'red', width:2}), {}, 'Faults');

// 3) Top-15 % mask
var sobelThr = ee.Number(
  sobelMag.reduceRegion({
    reducer: ee.Reducer.percentile([TOP_PCTL]),
    geometry: AOI.geometry(),
    scale: SCALE, bestEffort:true, tileScale:4, maxPixels: MAXPX
  }).get('sobelMag')
);
print('Sobel threshold (P'+TOP_PCTL.format()+' ~ top-'+TOP_SHARE+'%):', sobelThr);

var sobelMask = sobelMag.gte(sobelThr).selfMask();
Map.addLayer(sobelMask, {palette:['#ffff00']}, 'Sobel top-'+TOP_SHARE+'% mask');

// 4) Orientation 0–180° (gradient + 90° → lineament strike)
var theta = sobelY.atan2(sobelX).multiply(180/Math.PI).add(360).mod(180);
var orientation = theta.add(90).mod(180).updateMask(sobelMask).rename('orientation');

// 5) Manual orientation histogram (15° bins)
var histFC = ee.FeatureCollection(
  BINS.map(function(start){
    start = ee.Number(start);
    var end = start.add(BIN_SIZE);
    var inBin = orientation.gte(start).and(orientation.lt(end));
    var cnt = orientation.updateMask(inBin).reduceRegion({
      reducer: ee.Reducer.count(),
      geometry: AOI.geometry(),
      scale: SCALE, bestEffort:true, tileScale:4, maxPixels: MAXPX
    }).get('orientation');
    cnt = ee.Number(ee.Algorithms.If(cnt, cnt, 0));
    return ee.Feature(null, {binStart: start, binEnd: end, count: cnt});
  })
);
print('Orientation histogram (manual bins):', histFC);

print(ui.Chart.feature.byFeature(histFC, 'binStart', 'count')
  .setOptions({
    title:'Sobel lineation orientation (0–180°)',
    hAxis:{title:'Azimuth bin start (°)'},
    vAxis:{title:'Pixel count'},
    legend:'none',
    bars:'vertical'
  })
);

// 6) Validation vs. fault buffers (1–3–6 km)
var areaImg = ee.Image.pixelArea();

function overlapMetrics(bufferMeters){
  var bufGeom = FAULTS.geometry().buffer(bufferMeters).intersection(AOI.geometry(), 1);

  var A_sobel = ee.Number(
    areaImg.updateMask(sobelMask).reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: AOI.geometry(),
      scale: SCALE, bestEffort:true, tileScale:4, maxPixels: MAXPX
    }).get('area')
  ).divide(1e6); // km²

  var A_fault = ee.Number(
    areaImg.updateMask(ee.Image.constant(1)).reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: bufGeom,
      scale: SCALE, bestEffort:true, tileScale:4, maxPixels: MAXPX
    }).get('area')
  ).divide(1e6); // km²

  var A_intersect = ee.Number(
    areaImg.updateMask(sobelMask).reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: bufGeom,
      scale: SCALE, bestEffort:true, tileScale:4, maxPixels: MAXPX
    }).get('area')
  ).divide(1e6); // km²

  var precision = A_intersect.divide(A_sobel).multiply(100);
  var recall    = A_intersect.divide(A_fault).multiply(100);
  var f1        = precision.multiply(recall).multiply(2).divide(precision.add(recall));

  return ee.Feature(null, {
    buffer_km: bufferMeters/1000,
    area_sobel_km2: A_sobel,
    area_fault_km2: A_fault,
    area_intersect_km2: A_intersect,
    precision_pct: precision,
    recall_pct: recall,
    f1_pct: f1
  });
}

var metricsFC = ee.FeatureCollection(BUF_LIST_M.map(overlapMetrics));
print('Precision / Recall / F1 table:', metricsFC);

// 7) Area share of top-zone (≈13 %)
var areaTop = areaImg.updateMask(sobelMask).reduceRegion({
  reducer: ee.Reducer.sum(), geometry: AOI, scale: SCALE,
  bestEffort:true, tileScale:4, maxPixels: MAXPX
}).getNumber('area');
var areaAOI = areaImg.reduceRegion({
  reducer: ee.Reducer.sum(), geometry: AOI, scale: SCALE,
  bestEffort:true, tileScale:4, maxPixels: MAXPX
}).getNumber('area');
var areaSharePct = areaTop.divide(areaAOI).multiply(100);
print('Top-zone area share in AOI (%):', areaSharePct);

// 8) Exports
var REGION = AOI.geometry();

// RAW gradient
Export.image.toDrive({
  image: sobelMag,
  description:'SobelMag_RAW',
  fileNamePrefix:'SobelMag_RAW',
  region: REGION, scale: SCALE, crs: CRS,
  fileFormat:'GeoTIFF', maxPixels: MAXPX
});

// VIS gradient map (electric palette)
Export.image.toDrive({
  image: sobelMag.visualize({
    min:0, max:100,
    palette:['#000000','#440154','#3b528b','#21918c','#5ec962','#fde725']
  }),
  description:'SobelMag_VIS',
  fileNamePrefix:'SobelMag_VIS',
  region: REGION, scale: SCALE, crs: CRS,
  fileFormat:'GeoTIFF', maxPixels: MAXPX
});

// Binary top-15 % mask
Export.image.toDrive({
  image: sobelMask,
  description:'Sobel_topMask',
  fileNamePrefix:'Sobel_topMask',
  region: REGION, scale: SCALE, crs: CRS,
  fileFormat:'GeoTIFF', maxPixels: MAXPX
});

// Orientation map 0–180°
Export.image.toDrive({
  image: orientation.visualize({
    min:0, max:180,
    palette:['#00204D','#005792','#00B2CA','#7DCFB6','#FBD1A2','#F79256','#F05D5E']
  }),
  description:'Sobel_orientation_VIS',
  fileNamePrefix:'Sobel_orientation_VIS',
  region: REGION, scale: SCALE, crs: CRS,
  fileFormat:'GeoTIFF', maxPixels: MAXPX
});

// Orientation histogram and validation tables
Export.table.toDrive({
  collection: histFC,
  description:'Sobel_orientation_table',
  fileNamePrefix:'Sobel_orientation_table',
  fileFormat:'CSV'
});

Export.table.toDrive({
  collection: metricsFC,
  description:'Sobel_overlap_metrics',
  fileNamePrefix:'Sobel_overlap_metrics',
  fileFormat:'CSV'
});
