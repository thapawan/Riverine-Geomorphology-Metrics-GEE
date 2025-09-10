/**** GEE — Seasonal water masks for AL rivers (DW + AWEI cross-check; robust & error-free)
 * Refs:
 * - Dynamic World V1 (official ID): GOOGLE/DYNAMICWORLD/V1 (Brown et al., 2022)
 * https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_DYNAMICWORLD_V1
 * - JRC Global Surface Water (occurrence): validated ref for water (Pekel et al., 2016)
 * - Accuracy: confusion matrix & Kappa (Congalton, 1991); F1/IoU (ESSD practice)
 * - USGS WDFN for hydrograph-derived seasons: https://doi.org/10.5066/F7P55KJN
 ****/

// ---------------------- USER SETTINGS ----------------------
var EXPORT_FOLDER     = 'AL_River_WaterMasks'; // Google Drive folder
var BUFFER_METERS     = 500;                    // Buffer if river asset is line-like
var DW_PROB_THRESH    = 0.50;                   // DW 'water' probability threshold (tune if needed)
var JRC_OCC_THRESH    = 50;                     // JRC stable-water: occurrence > 50%
var VALID_N_PERCLASS  = 2000;                   // Validation samples per class (water/nonwater)
var YEARS             = [2016, 2020, 2024];
var SEASONS           = ['dry', 'wet'];

// AWEI (with-shadow) settings (Feyisa et al., 2014)
var USE_AWEI          = true;                   // set false if you want DW-only
var AWEI_VARIANT      = 'AWEI_sh';              // 'AWEI_sh' or 'AWEI_nsh'
var AWEI_THRESH       = 0.0;                    // AWEI > 0 => water (tune to ±0.05 if needed)

// List of all rivers to run
var RIVERS_TO_RUN = ['SipseyRiver', 'BlackWarriorRiver', 'CahabaRiver', 'CoosaRiver'];

// ---------------------- YOUR ASSETS ------------------------
var RiversFC = {
  'BlackWarriorRiver': ee.FeatureCollection('projects/skilful-boulder-440618-e7/assets/BlackWarriorRiver'),
  'CahabaRiver'      : ee.FeatureCollection('projects/skilful-boulder-440618-e7/assets/CahabaRiver'),
  'SipseyRiver'      : ee.FeatureCollection('projects/skilful-boulder-440618-e7/assets/SipseyRiver'),
  'CoosaRiver'       : ee.FeatureCollection('projects/skilful-boulder-440618-e7/assets/CoosaRiver')
};
var GAGES = ee.FeatureCollection('projects/skilful-boulder-440618-e7/assets/USGS_GaugeStations');

// ---------------------- DATASETS ---------------------------
var DW  = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1');
var JRC = ee.Image('JRC/GSW1_4/GlobalSurfaceWater');
var S2  = ee.ImageCollection('COPERNICUS/S2_SR');

// ---------------------- SEASON WINDOWS --------------------
var seasonWindows = {
  'SipseyRiver': {
    2016: {dry:['2016-08-01','2016-10-31'], wet:['2016-01-01','2016-04-30']},
    2020: {dry:['2020-08-01','2020-10-31'], wet:['2020-01-01','2020-04-30']},
    2024: {dry:['2024-08-01','2024-10-31'], wet:['2024-01-01','2024-04-30']}
  },
  'BlackWarriorRiver': {
    2016: {dry:['2016-08-01','2016-10-31'], wet:['2016-01-01','2016-04-30']},
    2020: {dry:['2020-08-01','2020-10-31'], wet:['2020-01-01','2020-04-30']},
    2024: {dry:['2024-08-01','2024-10-31'], wet:['2024-01-01','2024-04-30']}
  },
  'CahabaRiver': {
    2016: {dry:['2016-08-01','2016-10-31'], wet:['2016-01-01','2016-04-30']},
    2020: {dry:['2020-08-01','2020-10-31'], wet:['2020-01-01','2020-04-30']},
    2024: {dry:['2024-08-01','2024-10-31'], wet:['2024-01-01','2024-04-30']}
  },
  'CoosaRiver': {
    2016: {dry:['2016-08-01','2016-10-31'], wet:['2016-01-01','2016-04-30']},
    2020: {dry:['2020-08-01','2020-10-31'], wet:['2020-01-01','2020-04-30']},
    2024: {dry:['2024-08-01','2024-10-31'], wet:['2024-01-01','2024-04-30']}
  }
};

// ---------------------- HELPERS ---------------------------
function corridorGeometry(fc) {
  var geom = fc.geometry();
  var hasArea = ee.Number(geom.area(1)).gt(0);
  var corridor = ee.Algorithms.If(hasArea, geom, geom.buffer(BUFFER_METERS));
  return ee.Geometry(corridor).simplify(1);
}

function dwSeason(aoi, start, end) {
  return DW.filterBounds(aoi).filterDate(start, end);
}
function dwWaterComposite(dwCol) { return dwCol.select('water').median(); }
function dwBinary(probImg, thr)    { return probImg.gte(thr).rename('dw_water'); }

function jrcStable(aoi, occThr) {
  var occ = JRC.select('occurrence').clip(aoi);
  var water = occ.gte(occThr);
  var nonwater = occ.lte(1);
  return water.rename('jrc_water').updateMask(water.or(nonwater));
}

// ---- AWEI (Feyisa et al., 2014) ----
function maskS2SR(im){
  var scl = im.select('SCL');
  var mask = scl.neq(0).and(scl.neq(1)).and(scl.neq(3))
               .and(scl.neq(8)).and(scl.neq(9))
               .and(scl.neq(10)).and(scl.neq(11));
  return im.updateMask(mask);
}
function awei_sh(im){
  var G=im.select('B3').divide(10000), N=im.select('B8').divide(10000),
      SW1=im.select('B11').divide(10000), SW2=im.select('B12').divide(10000);
  return G.subtract(SW1).multiply(4).subtract(N.multiply(0.25).add(SW2.multiply(2.75))).rename('AWEI');
}
function awei_nsh(im){
  var B=im.select('B2').divide(10000), G=im.select('B3').divide(10000),
      N=im.select('B8').divide(10000), SW1=im.select('B11').divide(10000), SW2=im.select('B12').divide(10000);
  return B.add(G.multiply(2.5)).subtract(N.multiply(1.5).add(SW1).add(SW2.multiply(0.25))).rename('AWEI');
}
function buildAWEISeason(aoi, start, end, variant, thr){
  var s2 = S2.filterBounds(aoi).filterDate(start, end).map(maskS2SR);
  var n  = s2.size();
  // If empty (e.g., early 2016 L2A gaps), skip AWEI gracefully.
  var aweiMed = ee.Algorithms.If(n.gt(0),
    ee.ImageCollection(s2.map(function(im){ return (variant==='AWEI_sh'? awei_sh(im): awei_nsh(im))
                                           .copyProperties(im, im.propertyNames()); }))
      .median().clip(aoi),
    ee.Image(0).updateMask(ee.Image(0)).rename('AWEI') // empty image
  );
  var aweiMask = ee.Image(aweiMed).gte(thr).rename('awei_water');
  return {awei: ee.Image(aweiMed), mask: aweiMask, count: n};
}

function stratifiedValidate(aoi, pred01, ref01, nPerClass){
  var stack = ee.Image.cat(ref01.rename('ref'), pred01.rename('pred')).toInt();
  var samples = stack.stratifiedSample({
    numPoints: nPerClass,
    classBand: 'ref',
    region: aoi,
    scale: 10,
    geometries: true,
    seed: 42
  });
  var cm = samples.errorMatrix('ref','pred');
  return {samples:samples, cm:cm};
}
function exportMask(img, region, desc){
  Export.image.toDrive({image: img.toByte(), description: desc, folder: EXPORT_FOLDER, region: region, scale: 10, maxPixels: 1e13});
}
function exportCSV(fc, desc){
  Export.table.toDrive({collection: fc, description: desc, folder: EXPORT_FOLDER, fileFormat: 'CSV'});
}

// ---------------------- RUN ALL RIVERS ---------------------
// Use a function to process each river to keep the code clean and reusable.
var processRiver = function(RIVER_TO_RUN) {
  var fc = RiversFC[RIVER_TO_RUN];
  var aoi = corridorGeometry(fc);

  // Show & export gauges used for hydrographs (for the paper figure/provenance)
  var riverGages = GAGES.filterBounds(aoi);
  Map.addLayer(riverGages, {color:'red'}, RIVER_TO_RUN+' gauges');
  exportCSV(riverGages, RIVER_TO_RUN + '_USGS_Gauges');

  // Collect per-season metrics for summary CSV
  var metrics = [];

  YEARS.forEach(function(year){
    SEASONS.forEach(function(season){
      var win = seasonWindows[RIVER_TO_RUN][year][season];
      var start = ee.Date(win[0]), end = ee.Date(win[1]).advance(1,'day');

      // 1) Dynamic World seasonal
      var dwCol = dwSeason(aoi, start, end);
      var dwProb = dwWaterComposite(dwCol).clip(aoi);
      var dwMask = dwBinary(dwProb, DW_PROB_THRESH);

      // 2) Optional AWEI seasonal (skip silently if S2 L2A empty)
      var aweiOut = (USE_AWEI ? buildAWEISeason(aoi, start, end, AWEI_VARIANT, AWEI_THRESH) : null);
      var hasAWEI = (USE_AWEI ? ee.Number(aweiOut.count).gt(0) : ee.Number(0));

      // 3) JRC stable-water reference
      var jrcRef = jrcStable(aoi, JRC_OCC_THRESH);

      // 4) Validation (DW vs JRC) — per your note: “validation one is fine”
      var pred01 = dwMask.unmask(0).rename('pred');
      var ref01  = jrcRef.unmask(0).rename('ref');
      var val    = stratifiedValidate(aoi, pred01, ref01, VALID_N_PERCLASS);
      var cm     = val.cm;
      var arr    = cm.array();

      var oa     = ee.Number(cm.accuracy());
      var kappa  = ee.Number(cm.kappa()); // Congalton (1991)
      var recall1    = ee.Array(cm.producersAccuracy()).get([1]);
      var precision1 = ee.Array(cm.consumersAccuracy()).get([1]);
      var f1_1 = ee.Number(precision1).multiply(recall1).multiply(2).divide(ee.Number(precision1).add(recall1));
      var tp = ee.Number(arr.get([1,1])), fn = ee.Number(arr.get([1,0])), fp = ee.Number(arr.get([0,1]));
      var iou1 = tp.divide(tp.add(fn).add(fp)); // IoU (water)

      // 5) Exports
      var base = RIVER_TO_RUN + '_' + year + '_' + season;
      exportMask(dwMask, aoi, base + '_DW_WaterMask');
      exportCSV(val.samples, base + '_ValSamples_DW_vs_JRC');

      if (USE_AWEI){
        var aweiMask = ee.Image(aweiOut.mask);
        print('S2 L2A image count (for AWEI)', base, aweiOut.count);
        var maskToExport = ee.Image(aweiMask).updateMask(aweiMask);
        var disagree = dwMask.neq(aweiMask).selfMask().rename('DW_XOR_AWEI');
        exportMask(maskToExport, aoi, base + '_AWEI_WaterMask');
        exportMask(disagree,      aoi, base + '_DWxAWEI_Disagree');
      }

      var row = ee.Feature(null, {
        river: RIVER_TO_RUN, year: year, season: season,
        start: win[0], end: win[1],
        dw_prob_thresh: DW_PROB_THRESH, jrc_occ_thresh: JRC_OCC_THRESH,
        OA: oa, Kappa: kappa, Precision_water: precision1, Recall_water: recall1,
        F1_water: f1_1, IoU_water: iou1, nSamples: val.samples.size()
      });
      metrics.push(row);

      // Console for sanity check
      print('Confusion matrix (DW vs JRC)', base, cm);
      print('Metrics', base, {OA:oa, Kappa:kappa, Precision:precision1, Recall:recall1, F1:f1_1, IoU:iou1});
    });
  });

  // Export per-river summary metrics
  exportCSV(ee.FeatureCollection(metrics), RIVER_TO_RUN + '_DW_JRC_SummaryMetrics');

  // Map preview - only show the last river for preview
  Map.centerObject(aoi, 9);
  Map.addLayer(aoi, {color:'yellow'}, RIVER_TO_RUN + ' corridor');
  Map.addLayer(JRC.select('occurrence'), {min:0,max:100,palette:['000000','0000ff']}, 'JRC occurrence');
};

// Run the analysis for each river in the list
RIVERS_TO_RUN.forEach(processRiver);
