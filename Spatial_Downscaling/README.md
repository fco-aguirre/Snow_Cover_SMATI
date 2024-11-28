## This script is to apply the Spatial Downscaling into your files:

### Befor to run this script
1. Inside DATA/*<name_project>* folder you must replace the names of your spatial files indicated in the Parameters_Data, as is shown in the Brunswick example
  You need to prepare the next files:
  - DEM (GeoTiff) with the same areal extension that your MODIS download area. 
    - Put this file inside of DEM folder: ***Snow_cover_SMATI/DATA/name_project/DEM/***
  - The Water Mask (GeoTiff) with the same areal extension that your MODIS download area. 
    - Put this file inside of land_water_mask folder: ***Snow_cover_SMATI/name_project/DATA/land_water_mask/***
  - The watershed (shapefile) where the study will be performance (in UTM projection)
    - Put this file inside of land_water_mask folder: ***Snow_cover_SMATI/name_project/DATA/Watershed/***

### To apply this script:

**Run in the terminal:**

    Rscript Downscaling_parall_v6.r -s <Folder_name> -y <year>

1. -s DIR, source_file DIR
- Corresponds to the *<name_project>* of the folder where the downloaded files are ordered inside of the DATA directory
- To use the example dataset, this name is *'Brunswick'*

2. -y YEAR, year YEAR, the year that you need to process
- To use the example dataset, this year is *2015*