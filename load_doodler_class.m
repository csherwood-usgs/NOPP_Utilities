% load_doodler_class - Reads the .nc file containing Doodler classes
%
% the coordinate system is in the attributes, but easy to read with ncdump
% 
% (CRS) C:\crs\src\NOPP_Utilities>ncdump naip3_crs_label.nc | more
% netcdf naip3_crs_label {
% dimensions:
%         band = 1 ;
%         x = 1024 ;
%         y = 1024 ;
% variables:
%         int band(band) ;
%         double x(x) ;
%                 x:_FillValue = NaN ;
%         double y(y) ;
%                 y:_FillValue = NaN ;
%         int spatial_ref ;
%                 spatial_ref:crs_wkt = "PROJCS[\"NAD83 / UTM zone 16N\",GEOGCS[\"NAD83\",DATUM[\"North_American_Datum_1983\",SPHEROID[\"GRS 1980\",6378137,298.257222101,AUTHORITY[\"EPSG\",\"7019\"]],AUTHORITY[\"EPSG\",\"6269\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4269\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",-87],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"26916\"]]" ;
%                 spatial_ref:semi_major_axis = 6378137. ;
%                 spatial_ref:semi_minor_axis = 6356752.31414036 ;
%                 spatial_ref:inverse_flattening = 298.257222101 ;
%                 spatial_ref:reference_ellipsoid_name = "GRS 1980" ;
%                 spatial_ref:longitude_of_prime_meridian = 0. ;
%                 spatial_ref:prime_meridian_name = "Greenwich" ;
%                 spatial_ref:geographic_crs_name = "NAD83" ;
%                 spatial_ref:horizontal_datum_name = "North American Datum 1983" ;
%                 spatial_ref:projected_crs_name = "NAD83 / UTM zone 16N" ;
%                 spatial_ref:grid_mapping_name = "transverse_mercator" ;
%                 spatial_ref:latitude_of_projection_origin = 0. ;
%                 spatial_ref:longitude_of_central_meridian = -87. ;
%                 spatial_ref:false_easting = 500000. ;
%                 spatial_ref:false_northing = 0. ;
%                 spatial_ref:scale_factor_at_central_meridian = 0.9996 ;
%                 spatial_ref:spatial_ref = "PROJCS[\"NAD83 / UTM zone 16N\",GEOGCS[\"NAD83\",DATUM[\"North_American_Datum_1983\",SPHEROID[\"GRS 1980\",6378137,298.257222101,AUTHORITY[\"EPSG\",\"7019\"]],AUTHORITY[\"EPSG\",\"6269\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4269\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",-87],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"26916\"]]" ;
%                 spatial_ref:GeoTransform = "653891.0 1.0 0.0 3295343.0 0.0 -1.0" ;
%         ubyte __xarray_dataarray_variable__(band, y, x) ;
%                 __xarray_dataarray_variable__:scale_factor = 1. ;
%                 __xarray_dataarray_variable__:add_offset = 0. ;
%                 __xarray_dataarray_variable__:grid_mapping = "spatial_ref" ;
% 
nc_name = 'naip3_crs_label.nc'
y = ncread('naip3_crs_label.nc','y');
x = ncread('naip3_crs_label.nc','x');
class = ncread('naip3_crs_label.nc','class');

% I need to get these into the nc file...these are in order
class_names = {'water','sand','brushy','marsh', 'peat','anthro'}
% for some reason, you have to transform the matrix
pcolor(x,y,class')
shading flat
shg
