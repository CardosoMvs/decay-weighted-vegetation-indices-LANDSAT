/**
 * Temporal Decay Weighted Vegetation Indices Processor
 * Applies exponential decay weights to historical vegetation indices (NDVI/EVI/SAVI)
 * Creates annual composites incorporating 7 years of temporal context
 * Outputs processed indices for full year, wet season, and dry season
 * 
 * Input: MapBiomas Landsat mosaics
 * Output: Decay-weighted vegetation indices (1985-2023)
 *
 * Contact: cardoso.mvs@gmail.com
 */

// 1. Define Area of Interest (Brazilian biomes)
var biomas = ee.FeatureCollection('projects/mapbiomas-workspace/AUXILIAR/biomas_IBGE_250mil');
var aoi = biomas;
var aoi_img = ee.Image().paint(aoi).eq(0);
var aoi_bounds = aoi.geometry().bounds();

// 2. Landsat satellite mapping by year
var year_sat = {
  1985: 'l5', 1986: 'l5', 1987: 'l5', 1988: 'l5', 1989: 'l5', 
  1990: 'l5', 1991: 'l5', 1992: 'l5', 1993: 'l5', 1994: 'l5',
  1995: 'l5', 1996: 'l5', 1997: 'l5', 1998: 'l5', 1999: 'l5',
  2000: 'l5', 2001: 'l7', 2002: 'l7', 2003: 'l5', 2004: 'l5',
  2005: 'l5', 2006: 'l5', 2007: 'l5', 2008: 'l5', 2009: 'l5',
  2010: 'l5', 2011: 'l7', 2012: 'l7', 2013: 'l8', 2014: 'l8',
  2015: 'l8', 2016: 'l8', 2017: 'l8', 2018: 'l8', 2019: 'l8',
  2020: 'l8', 2021: 'l8', 2022: 'l8', 2023: 'l8'
};

// 3. Processing years (1985-2023)
var years = [
  1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,
  1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,
  2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,
  2015,2016,2017,2018,2019,2020,2021,2022,2023
];

// 4. Load MapBiomas Landsat mosaics
var mb_col = ee.ImageCollection('projects/nexgenmap/MapBiomas2/LANDSAT/BRAZIL/mosaics-2');

// 5. Vegetation indices to process
var mb_indices_year = [
  'ndvi_median', 'ndvi_median_wet', 'ndvi_median_dry',
  'evi2_median', 'evi2_median_wet', 'evi2_median_dry',
  'savi_median', 'savi_median_wet', 'savi_median_dry'
];

// 6. Main processing function
var index_col = years.map(function(year) {
  // Initialize containers for each index
  var container_ndvi = ee.Image().select();
  var container_ndvi_wet = ee.Image().select();
  var container_ndvi_dry = ee.Image().select();
  var container_evi2 = ee.Image().select();
  var container_evi2_wet = ee.Image().select();
  var container_evi2_dry = ee.Image().select();
  var container_savi = ee.Image().select();
  var container_savi_wet = ee.Image().select();
  var container_savi_dry = ee.Image().select();

  // Exponential decay factors (alpha=0.7, sum ~= 1)
  var alfa_07 = [0.327, 0.229, 0.160, 0.112, 0.078, 0.055, 0.038];

  // Apply decay to previous 7 years
  alfa_07.forEach(function(factor) {
    var year_decay = year - alfa_07.indexOf(factor);
    
    // Handle edge cases for first/last years
    if (year_decay < years[0]) year_decay = years[0];
    
    // Get neighboring years for temporal blending
    var year_prev = year_decay - 1;
    var year_post = year_decay + 1;
    
    if (year_decay === years[0]) {
      year_prev = year_decay + 1;
      year_post = year_decay + 2;
    }
    if (year_decay === years.slice(-1)[0]) {
      year_prev = year_decay - 1;
      year_post = year_decay - 2;
    }

    // Process each vegetation index
    mb_indices_year.forEach(function(index) {
      // Load index data for target and adjacent years
      var index_band = mb_col.filter(ee.Filter.eq('year', year_decay))
                           .filter(ee.Filter.eq('satellite', year_sat[year_decay]))
                           .select([index])
                           .mosaic();
      
      var index_band_prev = mb_col.filter(ee.Filter.eq('year', year_prev))
                                .filter(ee.Filter.eq('satellite', year_sat[year_prev]))
                                .select([index])
                                .mosaic();
      
      var index_band_post = mb_col.filter(ee.Filter.eq('year', year_post))
                                .filter(ee.Filter.eq('satellite', year_sat[year_post]))
                                .select([index])
                                .mosaic();

      // Temporal blending of indices
      var index_blend = index_band_post
        .blend(index_band_prev)
        .blend(index_band);

      // Apply decay weight and add to container
      var index_mean = index_blend.multiply(factor);
      
      // Add to appropriate container based on index type
      switch (index) {
        case 'ndvi_median':
          container_ndvi = container_ndvi.addBands(index_mean);
          break;
        case 'ndvi_median_wet':
          container_ndvi_wet = container_ndvi_wet.addBands(index_mean);
          break;
        case 'ndvi_median_dry':
          container_ndvi_dry = container_ndvi_dry.addBands(index_mean);
          break;
        case 'evi2_median':
          container_evi2 = container_evi2.addBands(index_mean);
          break;
        case 'evi2_median_wet':
          container_evi2_wet = container_evi2_wet.addBands(index_mean);
          break;
        case 'evi2_median_dry':
          container_evi2_dry = container_evi2_dry.addBands(index_mean);
          break;
        case 'savi_median':
          container_savi = container_savi.addBands(index_mean);
          break;
        case 'savi_median_wet':
          container_savi_wet = container_savi_wet.addBands(index_mean);
          break;
        case 'savi_median_dry':
          container_savi_dry = container_savi_dry.addBands(index_mean);
          break;
      }
    });
  });

  // Sum weighted contributions and convert to byte (0-100 scale)
  container_ndvi = container_ndvi.reduce('sum').divide(100).byte().rename('mb_ndvi_median');
  container_ndvi_wet = container_ndvi_wet.reduce('sum').divide(100).byte().rename('mb_ndvi_median_wet');
  container_ndvi_dry = container_ndvi_dry.reduce('sum').divide(100).byte().rename('mb_ndvi_median_dry');
  container_evi2 = container_evi2.reduce('sum').divide(100).byte().rename('mb_evi2_median');
  container_evi2_wet = container_evi2_wet.reduce('sum').divide(100).byte().rename('mb_evi2_median_wet');
  container_evi2_dry = container_evi2_dry.reduce('sum').divide(100).byte().rename('mb_evi2_median_dry');
  container_savi = container_savi.reduce('sum').divide(100).byte().rename('mb_savi_median');
  container_savi_wet = container_savi_wet.reduce('sum').divide(100).byte().rename('mb_savi_median_wet');
  container_savi_dry = container_savi_dry.reduce('sum').divide(100).byte().rename('mb_savi_median_dry');

  // Combine all indices into single image
  var container = container_ndvi.addBands([
    container_ndvi_wet, container_ndvi_dry,
    container_evi2, container_evi2_wet, container_evi2_dry,
    container_savi, container_savi_wet, container_savi_dry
  ]);

  // Set image metadata
  var start = ee.Date('' + year + '-01-01').millis();
  var end = ee.Date('' + (year + 1) + '-01-01').millis();

  return container.set({
    year: year,
    'system:time_start': start,
    'system:time_end': end
  });
});

// 7. Create final ImageCollection
index_col = ee.ImageCollection(index_col);
print('Processed Vegetation Indices Collection', index_col);

// 8. Export results
years.forEach(function(year) {
  var image = index_col.filter(ee.Filter.eq('year', year)).first();
  var description = '' + year;
  var assetId = 'projects/mapbiomas-workspace/SOLOS/COVARIAVEIS/LANDSAT_MB_INDICES_DECAY/' + description;

  // Add to map for visualization
  Map.addLayer(image.select('mb_ndvi_median'), {min: 0, max: 100, palette: ['red', 'yellow', 'green']}, 'NDVI ' + description);

  // Export to Earth Engine Assets
  Export.image.toAsset({
    image: image,
    description: 'GT_MC_SOLO-LANDSAT_MB_INDICES_DECAY-' + description,
    assetId: assetId,
    pyramidingPolicy: 'median',
    region: aoi_bounds,
    scale: 30,
    maxPixels: 1e13
  });
});

// Add map controls
Map.centerObject(aoi, 4);
