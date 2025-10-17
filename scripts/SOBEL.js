/***************************************************************
* SOBEL.js — Robust Morphostructural Gradient Analysis
* 
* This script computes the Sobel gradient magnitude and orientation
* from the MERIT DEM, performs percentile thresholding, 
* generates azimuth histograms, and evaluates spatial overlap 
* with known fault buffers.
*
* Author: Dilyara Nazyrova, 2025
* Project: Sobel-Morphostructure-KZ
***************************************************************/

// 0) INPUTS
var AOI    = ee.FeatureCollection('projects/ee-euniversity2024/assets/AOI');
var DEM    = ee.Image('MERIT/DEM/v1_0_3').select('dem').clip(AOI);
var FAULTS = ee.FeatureCollection('projects/ee-euniversity2024/assets/Razlom');
Map.centerObject(AOI, 8);

// 1) PARAMETERS
var SCALE       = 90;
var CRS         = 'EPSG:4326';
var MAXPX       = 1e13;
var TOP_PCTL    = 70;                 // top 30% (for visualization stability)
var BUF_LIST_M  = [1000, 5000, 10000];
var BIN_SIZE    = 15;                 // histogram bin width (°)
var BINS        = ee.List.sequence(0, 180 - BIN_SIZE, BIN_SIZE);

// 2) SOBEL FILTER (X, Y, Magnitude)
var kx = ee.Kernel.fixed(3,3,[[-1,0,1],[-2,0,2],[-1,0,1]]);
var ky = ee.Kernel.fixed(3,3,[[-1,-2,-1],[0,0,0],[1,2,1]]);
var sobelX = DEM.convolve(kx);
var sobelY = DEM.convolve(ky);
var sobelMag = sobelX.pow(2).add(sobelY.pow(2)).sqrt().rename('sobelMag');

Map.addLayer(sobelMag, {min:0, max:100, palette:['black','white']}, 'Sobel magnitude');
Map.addLayer(FAULTS.style({color:'red', width:2}), {}, 'Faults');

// 3) PERCENTILE MASK
var sobelThr = ee.Number(
  sobelMag.reduceRegion({
    reducer: ee.Reducer.percentile([TOP_PCTL]),
    geometry: AOI.geometry(),
    scale: SCALE, bestEffort:true, tileScale:4, maxPixels: MAXPX
  }).get('sobelMag')
);
var sobelMask = sobelMag.gte(sobelThr).selfMask();
Map.addLayer(sobelMask, {palette:['#ffff00']}, 'Sobel top30% mask');

// 4) ORIENTATION (0–180°)
var theta = sobelY.atan2(sobelX).multiply(180/Math.PI).add(360).mod(180);
var orientation = theta.add(90).mod(180).updateMask(sobelMask).rename('orientation');

// 5) MANUAL HISTOGRAM (no fixedHistogram)
var histFC = ee.FeatureCollection(
  BINS.map(function(start){
    start = ee.Number(start);
    var end = start.add(BIN_SIZE);
    var inBin = orientation.gte(start).and(orientation.lt(end));
    var cnt = orientation.updateMask(inBin).reduceRegion({
      reducer: ee.Reducer.count(),
      geometry: AOI.geometry(),
      scale: SCALE, bestEffort:true, tileScale: 4, maxPixels: MAXPX
    }).get('orientation');
    cnt = ee.Number(ee.Algorithms.If(cnt, cnt, 0)); // handle empty bins
    return ee.Feature(null, {binStart: start, binEnd: end, count: cnt});
  })
);

print('Sobel orientation histogram (manual bins):', histFC);
print(ui.Chart.feature.byFeature(histFC, 'binStart', 'count')
  .setOptions({
    title:'Sobel lineation orientation (0–180°)',
    hAxis:{title:'Azimuth bin start (°)'},
    vAxis:{title:'Pixel count'},
    legend:'none',
    bars:'vertical'
  })
);

// 6) OVERLAP ANALYSIS (1, 5, 10 km)
var areaImg = ee.Image.pixelArea();

function overlapStats(bufferMeters){
  var bufGeom = FAULTS.geometry().buffer(bufferMeters);

  var A_total = ee.Number(
    areaImg.updateMask(sobelMask).reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: AOI.geometry(),
      scale: SCALE, bestEffort:true, tileScale:4, maxPixels: MAXPX
    }).get('area')
  ).divide(1e6); // km²

  var A_inbuf = ee.Number(
    areaImg.updateMask(sobelMask).reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: bufGeom,
      scale: SCALE, bestEffort:true, tileScale:4, maxPixels: MAXPX
    }).get('area')
  ).divide(1e6); // km²

  var pct = ee.Algorithms.If(A_total.gt(0), A_inbuf.divide(A_total).multiply(100), 0);

  return ee.Feature(null, {
    buffer_km: bufferMeters/1000,
    area_total_km2: A_total,
    area_inbuf_km2: A_inbuf,
    overlap_pct: ee.Number(pct)
  });
}

var overlapFC = ee.FeatureCollection(BUF_LIST_M.map(overlapStats));
print('Overlap table (Sobel mask vs fault buffers):', overlapFC);

// 7) EXPORTS
var REGION = AOI.geometry();

Export.image.toDrive({
  image: sobelMag,
  description:'SobelMag_RAW',
  fileNamePrefix:'SobelMag_RAW',
  region: REGION, scale: SCALE, crs: CRS, fileFormat:'GeoTIFF', maxPixels: MAXPX
});
Export.image.toDrive({
  image: sobelMag.visualize({min:0, max:100, palette:['#000000','#440154','#3b528b','#21918c','#5ec962','#fde725']}),
  description:'SobelMag_VIS',
  fileNamePrefix:'SobelMag_VIS',
  region: REGION, scale: SCALE, crs: CRS, fileFormat:'GeoTIFF', maxPixels: MAXPX
});
Export.image.toDrive({
  image: sobelMask,
  description:'Sobel_topMask',
  fileNamePrefix:'Sobel_topMask',
  region: REGION, scale: SCALE, crs: CRS, fileFormat:'GeoTIFF', maxPixels: MAXPX
});
Export.image.toDrive({
  image: orientation.visualize({min:0, max:180, palette:['#00204D','#005792','#00B2CA','#7DCFB6','#FBD1A2','#F79256','#F05D5E']}),
  description:'Sobel_orientation_VIS',
  fileNamePrefix:'Sobel_orientation_VIS',
  region: REGION, scale: SCALE, crs: CRS, fileFormat:'GeoTIFF', maxPixels: MAXPX
});

// EXPORT TABLES
Export.table.toDrive({
  collection: histFC,
  description: 'Sobel_orientation_table',
  fileNamePrefix: 'Sobel_orientation_table',
  fileFormat: 'CSV'
});
Export.table.toDrive({
  collection: overlapFC,
  description: 'Sobel_overlap_table',
  fileNamePrefix: 'Sobel_overlap_table',
  fileFormat: 'CSV'
});
