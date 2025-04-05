# 🌱 Temporal Decay Weighted Vegetation Indices Processor

## 📝 Description  
This repository contains a Google Earth Engine (GEE) JavaScript script that processes vegetation indices (NDVI, EVI, and SAVI) from MapBiomas Landsat mosaics by applying temporal decay weights to create smoothed annual composites. The approach incorporates historical context by weighting previous years' data using an exponential decay function.

## 🔍 Key Features
- **Temporal Smoothing**: Applies exponential decay weights (alpha=0.7) to 7 years of historical data
- **Multi-Index Processing**: Handles NDVI, EVI, and SAVI indices simultaneously
- **Seasonal Analysis**: Generates separate composites for:
  - Full year
  - Wet season
  - Dry season
- **Sensor-Aware Processing**: Automatically handles Landsat 5, 7, and 8 sensor differences
- **Large-Scale Processing**: Covers all Brazilian biomes from 1985-2023
- **Optimized Output**: Converts results to efficient byte format (0-100 scale)

## 🛠️ Technical Details
- **Input Data**: MapBiomas Landsat mosaics (projects/nexgenmap/MapBiomas2/LANDSAT/BRAZIL/mosaics-2)
- **Output**: Annual composites stored as Earth Engine assets
- **Processing**: 
  - Temporal blending of adjacent years
  - Weighted summation of decayed values
  - Byte conversion for efficient storage

## 📊 Applications
The resulting decay-weighted indices are particularly useful for:
- Vegetation trend analysis
- Change detection algorithms
- Time-series modeling
- Ecological monitoring
- Land cover classification

## 🚀 Usage
The script is ready to run in the Google Earth Engine Code Editor. Simply copy and paste into a new script, and the results will be exported to your specified Earth Engine asset location.
