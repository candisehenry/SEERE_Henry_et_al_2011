# SEERE_Henry_et_al_2011

General SEERE set up:
1. Model should be run using run_SEERE_guatemala.py. It will call other scripts required.
2. Script d_setup_guatemala.py might need to be changed to call files from correct paths. This file is hardcoded to call specific file names related to the Guatemala analysis.
3. SEERE calls the other scripts in letter order (e.g., a_get_solarwinddata, then b_calc_solarpower, c_interp_solar, etc.). Scripts without lettering are supporting functions that get called by the letter scripts.



Setting up other data:

High-resolution population tiffs (e.g., grid-level population per km^2 data):
1. Open raster in QGIS.
2. Decrease grid scale (if high resolution) by:
    a) layer > add layer > raster
    b) raster > align rasters > add raster to align > resampling method: nearest neighbor > select "rescale values according to cell size" > change cell size
    c) right click raster layer > export > save as geotiff file
    d) right click raster layer > properties > information > write down coordinate extent
3. Open downsized tiff file in Python using b_calc_enedem.py > calc_energydemand2().
4. Update grid cell size (step_pop) in d_setup.py.
5. Update file names in d_setup.py and run_SEERE.py with correct folders and files.
6. Run run_SEERE.py.

Solar and wind data:
1. Open a_get_solarwinddata.py.
2. Change dataURL to correct bounding box coordinates.
3. Run import_nasa_api(dataURL).
4. If error, print(results) and click on variable to see error. (Most likely coordinate bounds too big for API, and must run in chunks of smaller spatial scales.)
5. If/once no error, copy output URL into internet browser. CSV should automatically download.
6. Update grid cell size (step) in d_setup.py.
7. Update file names in d_setup.py and run_SEERE.py with correct folders and files.
8. Run run_SEACART.py. (Use .npy outputs for data analysis in Python. Use .asc files for mapping in QGIS. [See next step for details.])
9. Mapping in QGIS:
    a) Load output .asc files as raster.
    b) Use country shapefile to clip raster extent.

Hydro flow, slope, catchment area data:
1. Open streamlines shapefile in QGIS.
2. Open flow/slope/catchment area csvs in QGIS.
3. Combine shapefile and csv by:
    a) right-click shapefile layer > properties > joins > join layer = csv > join field = COMID > target field = flow/slope/area > custom field name prefix = empty
    b) right-click shapefile layer > properties > fields > edit > manually clean up unused fields
    c) right-click shapefile layer > export > save vector layer as ESRI shapefile *** Make sure CRS is WGS 84!! ***
    d) load newly saved shapefile in QGIS *** Important because QGIS can't make raster from joined data!! ***
4. Convert to from shapefile/vector to raster by:
    a) raster > conversion > rasterize > input layer = newly loaded shapefile > burn-in value = flow/slope/area > output raster size units = georeferenced units > resolution = step > output extent = use layer extent 
    b) if error, go to view > panel > log messages to see error log.
    c) right-click raster layer > export > save as > format = geotiff > CRS = WGS 84 > extent = calculate from layer
    d) right click raster layer > properties > information > write down coordinate extent
5. Open tiff files in Python using a_get_hydrodata().
