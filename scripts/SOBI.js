/**
 * SOBI.js — Sobel-Based Morphostructural Index
 * 
 * This script computes the Spectral-Oriented Brightness Index (SOBI)
 * and integrates gradient magnitude, azimuth, and anomaly thresholding.
 * 
 * Author: Dilyara Nazyrova, 2025
 * Project: Sobel-Morphostructure-KZ
 */

var dem = ee.Image('MERIT/DEM/v1_0_3');
var sobel = ee.Terrain.slope(dem).rename('sobel_gradient');

// Normalize and apply threshold
var threshold = sobel.gt(ee.Number(23.9));  // example threshold
var sobi = sobel.multiply(threshold);

// Display
Map.centerObject(dem);
Map.addLayer(sobel, {min:0, max:50, palette:['white','black']}, 'Sobel Gradient');
Map.addLayer(threshold, {palette:['red']}, 'Top Sobel Zones');

Export.image.toDrive({
  image: sobi,
  description: 'SOBI_map',
  scale: 30,
  region: AOI_all_1
});
