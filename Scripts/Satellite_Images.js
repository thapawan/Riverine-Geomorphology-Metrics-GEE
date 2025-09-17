/**** CONFIG ***************************************************************/
var BUFFER_M = 500;                     // Buffer (m) around line to define AOI
var DRIVE_FOLDER = 'GEE_River_Exports'; // Google Drive export folder
var VIS_TRUECOLOR = {min: 0.03, max: 0.35, bands: ['red','green','blue'], gamma: 1.1};
var FORCE_10M_VIS = false;              // Optional: force 10 m reproject for on-map display

/**** GEOMETRIES: DEFINE RIVER SEGMENTS FROM START/END COORDS **************/
// Black Warrior River (~100 km)
var BWR_START = ee.Geometry.Point([-87.8053833333, 32.8907361111]); // 32°53'26.65"N, 87°48'19.38"W
var BWR_END   = ee.Geometry.Point([-87.8525944444, 32.5323000000]); // 32°31'56.28"N, 87°51'9.34"W
var BWR_LINE  = ee.Geometry.LineString([BWR_START.coordinates(), BWR_END.coordinates()]);
var BWR_AOI   = BWR_LINE.buffer(BUFFER_M); // use the actual buffer polygon

// Cahaba River (~50 km)
var CHB_START = ee.Geometry.Point([-87.1979333333, 32.5333850000]); // 32°32'0.19"N, 87°11'52.56"W
var CHB_END   = ee.Geometry.Point([-87.0938777778, 32.3180250000]); // 32°19'4.89"N, 87° 5'37.96"W
var CHB_LINE  = ee.Geometry.LineString([CHB_START.coordinates(), CHB_END.coordinates()]);
var CHB_AOI   = CHB_LINE.buffer(BUFFER_M);

// Sanity check lengths (km)
print('Black Warrior length (km):', BWR_LINE.length().divide(1000));
print('Cahaba length (km):',       CHB_LINE.length().divide(1000));

/**** LANDSAT 2000: MERGE L5 TM + L7 ETM+ (C2 L2 SR), UPSAMPLE ON EXPORT ***/
function scaleMaskLandsatC2(img, sensor) {
  // Common SR scale for C2 L2: SR = DN * 0.0000275 - 0.2
  var srBands = (sensor === 'L5')
    ? ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7']
    : ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7']; // same mapping for L7

  var sr = img.select(srBands).multiply(0.0000275).add(-0.2)
              .rename(['blue','green','red','nir','swir1','swir2']);

  var qa = img.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 1).eq(0)   // dilated cloud
           .and(qa.bitwiseAnd(1 << 3).eq(0)) // cloud
           .and(qa.bitwiseAnd(1 << 4).eq(0)) // cloud shadow
           .and(qa.bitwiseAnd(1 << 5).eq(0)); // snow

  return sr.updateMask(mask).copyProperties(img, img.propertyNames());
}

function landsat2000Composite(aoi) {
  var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
              .filterBounds(aoi).filterDate('2000-01-01', '2000-12-31')
              .map(function(i){ return scaleMaskLandsatC2(i, 'L5'); });

  var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
              .filterBounds(aoi).filterDate('2000-01-01', '2000-12-31')
              .map(function(i){ return scaleMaskLandsatC2(i, 'L7'); });

  var merged = l5.merge(l7);
  print('Landsat count (2000):', merged.size());

  return merged.median().clip(aoi);  // upsample happens at export (10 m)
}

/**** SENTINEL-2: HYBRID SR+TOA TO FILL GAPS (2016 & 2024) *****************/
// SR path (uses SCL classes 4,5,6,7,11)
function maskS2_SCL(img) {
  var scl = img.select('SCL');
  var good = scl.eq(4).or(scl.eq(5)).or(scl.eq(6)).or(scl.eq(7)).or(scl.eq(11));
  var sr = img.select(['B2','B3','B4','B8']).divide(10000)
              .rename(['blue','green','red','nir']);
  return sr.updateMask(good).copyProperties(img, img.propertyNames());
}
// TOA fallback (QA60 bits 10/11)
function maskS2_QA60(img) {
  var qa = img.select('QA60');
  var cloud = qa.bitwiseAnd(1 << 10).neq(0).or(qa.bitwiseAnd(1 << 11).neq(0));
  var toa = img.select(['B2','B3','B4','B8']).divide(10000)
               .rename(['blue','green','red','nir']);
  return toa.updateMask(cloud.not()).copyProperties(img, img.propertyNames());
}

function s2HybridComposite(aoi, start, end) {
  var s2sr = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
      .filterBounds(aoi).filterDate(start, end)
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 80))
      .map(maskS2_SCL);

  var s2toa = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
      .filterBounds(aoi).filterDate(start, end)
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 80))
      .map(maskS2_QA60);

  // If SR exists, merge SR + TOA to fill any SR gaps; else use TOA only.
  var merged = ee.ImageCollection(ee.Algorithms.If(
      s2sr.size().gt(0),
      s2sr.merge(s2toa),
      s2toa
  ));
  print('Sentinel-2 images (SR):', s2sr.size(), ' (TOA):', s2toa.size(),
        ' (merged):', merged.size(), ' for ', start.slice(0,4));

  return merged.median().clip(aoi);
}

/**** BUILD COMPOSITES *****************************************************/
var bwr_l7_2000 = landsat2000Composite(BWR_AOI);
var bwr_s2_2016 = s2HybridComposite(BWR_AOI, '2016-01-01', '2016-12-31');
var bwr_s2_2024 = s2HybridComposite(BWR_AOI, '2024-01-01', '2024-12-31');

var chb_l7_2000 = landsat2000Composite(CHB_AOI);
var chb_s2_2016 = s2HybridComposite(CHB_AOI, '2016-01-01', '2016-12-31');
var chb_s2_2024 = s2HybridComposite(CHB_AOI, '2024-01-01', '2024-12-31');

/**** OPTIONAL: FORCE 10 m VISUALIZATION ***********************************/
function maybeReproject10m(img) {
  return FORCE_10M_VIS ? img.reproject({crs: 'EPSG:32616', scale: 10}) : img;
}

/**** MAP LAYERS ***********************************************************/
Map.centerObject(BWR_AOI, 10);

// AOI outlines
Map.addLayer(ee.Image().paint(BWR_AOI, 1, 2), {palette:['yellow']}, 'BWR AOI (buffer)', false);
Map.addLayer(ee.Image().paint(CHB_AOI, 1, 2), {palette:['yellow']}, 'Cahaba AOI (buffer)', false);

// Composites (with optional 10 m reproject for display)
Map.addLayer(maybeReproject10m(bwr_l7_2000), VIS_TRUECOLOR, 'BWR – Landsat 2000 (upsampled on export)', true, 1.0);
Map.addLayer(maybeReproject10m(bwr_s2_2016), VIS_TRUECOLOR, 'BWR – Sentinel‑2 2016 (hybrid)', true, 1.0);
Map.addLayer(maybeReproject10m(bwr_s2_2024), VIS_TRUECOLOR, 'BWR – Sentinel‑2 2024', true, 1.0);

Map.addLayer(maybeReproject10m(chb_l7_2000), VIS_TRUECOLOR, 'Cahaba – Landsat 2000 (upsampled on export)', true, 1.0);
Map.addLayer(maybeReproject10m(chb_s2_2016), VIS_TRUECOLOR, 'Cahaba – Sentinel‑2 2016 (hybrid)', true, 1.0);
Map.addLayer(maybeReproject10m(chb_s2_2024), VIS_TRUECOLOR, 'Cahaba – Sentinel‑2 2024', true, 1.0);

// Quick coverage masks (white = has data after masking)
function coverageMask(img) {
  return img.mask().reduce(ee.Reducer.min()).selfMask().visualize({palette:['ffffff'], opacity:0.6});
}
Map.addLayer(coverageMask(bwr_s2_2016), {}, 'BWR 2016 coverage mask', false);
Map.addLayer(coverageMask(chb_s2_2016), {}, 'Cahaba 2016 coverage mask', false);

/**** EXPORTS (to Google Drive) ********************************************/
function exportComposite(img, aoi, riverCode, sensor, year) {
  Export.image.toDrive({
    image: img.resample('bilinear'), // ensure bilinear before export
    description: riverCode + '_' + sensor + '_' + year + '_10m',
    fileNamePrefix: riverCode + '_' + sensor + '_' + year + '_10m',
    folder: DRIVE_FOLDER,
    region: aoi,
    crs: 'EPSG:32616', // UTM Zone 16N
    scale: 10,
    maxPixels: 1e13
  });
}

// Black Warrior exports
exportComposite(bwr_l7_2000, BWR_AOI, 'BWR', 'Landsat',   2000);
exportComposite(bwr_s2_2016, BWR_AOI, 'BWR', 'Sentinel2', 2016);
exportComposite(bwr_s2_2024, BWR_AOI, 'BWR', 'Sentinel2', 2024);

// Cahaba exports
exportComposite(chb_l7_2000, CHB_AOI, 'CHB', 'Landsat',   2000);
exportComposite(chb_s2_2016, CHB_AOI, 'CHB', 'Sentinel2', 2016);
exportComposite(chb_s2_2024, CHB_AOI, 'CHB', 'Sentinel2', 2024);
