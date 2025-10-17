/**
 * BridgeMap.js — Combined Structural and Metallogenic Visualization
 * 
 * Integrates Sobel high-gradient zones with ore deposit locations
 * and computes enrichment ratios for visualization.
 * 
 * Author: Dilyara Nazyrova, 2025
 * Project: Sobel-Morphostructure-KZ
 */

var sobelZones = ee.Image('users/dilyaranazyrova/Sobel_Zones');
var deposits = ee.FeatureCollection('users/dilyaranazyrova/Deposits');

// Compute enrichment index E = (D_in/A_in) / (D_total/A_total)
var totalArea = sobelZones.geometry().area();
var totalDeposits = deposits.size();
var highGradient = sobelZones.gt(0.5);
var depositsIn = deposits.filterBounds(highGradient.geometry());
var enrichment = depositsIn.size().divide(totalDeposits)
  .divide(highGradient.geometry().area().divide(totalArea));

print('Enrichment index (E):', enrichment);

Map.centerObject(sobelZones);
Map.addLayer(highGradient, {palette:['red']}, 'Sobel High Gradient Zones');
Map.addLayer(deposits, {color:'yellow'}, 'Mineral Deposits');
