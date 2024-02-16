# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 14:50:20 2022F

@author: pearsonra
"""

import numpy
import xarray
import rioxarray
import rioxarray.merge
import rasterio
import geopandas
import shapely
import pathlib
import shutil
import OSMPythonTools.overpass
import pandas

OSM_CRS = "EPSG:4326"
FABDEM_NODATA = 0

def load_dem(dem_path: pathlib):
    """ Load in a DEM - ensure values are float32. """
    dem = rioxarray.open_rasterio(dem_path, masked=True, chunks=True).squeeze("band", drop=True)
    dem = dem.rio.set_nodata(numpy.nan)
    dem = dem.rio.write_nodata(dem.rio.nodata)
    dem['x'] = dem.x.astype(numpy.float32)
    dem['y'] = dem.y.astype(numpy.float32)
    dem = dem.astype(numpy.float32)
    return dem

def resample_dem(dem_in: rioxarray, resolution: float, boundary: geopandas):
    """ A function to downscale a LiDAR DEM and combine with a FABDEM where there is no
    data to cover an island. """

    bounds = boundary.buffer(resolution).bounds
    if numpy.sign(float(dem_in.y[1] - dem_in.y[0])) > 0:
        new_y = numpy.arange(bounds.miny[0], bounds.maxy[0], resolution)
    else:
        new_y = numpy.arange(bounds.maxy[0], bounds.miny[0], -resolution)
    if numpy.sign(float(dem_in.x[1] - dem_in.x[0])) > 0:
        new_x = numpy.arange(bounds.minx[0], bounds.maxx[0], resolution)
    else:
        new_x = numpy.arange(bounds.maxx[0], bounds.minx[0], -resolution)
    new_y = new_y[new_y > min(float(dem_in.y[0]), float(dem_in.y[-1]))]
    new_y = new_y[new_y < max(float(dem_in.y[0]), float(dem_in.y[-1]))]
    new_x = new_x[new_x > min(float(dem_in.x[0]), float(dem_in.x[-1]))]
    new_x = new_x[new_x < max(float(dem_in.x[0]), float(dem_in.x[-1]))]

    dem_out = dem_in.interp(x=new_x, y=new_y, method="linear").rio.write_crs(
        dem_in.rio.crs
    )
    dem_out.rio.write_transform(dem_out.rio.transform(recalc=True), inplace=True)
    dem_out.rio.write_nodata(dem_in.rio.nodata, inplace=True)
    dem_out.rio.clip(boundary.geometry, all_touched=True, drop=False)

    return dem_out


def osm_query_islands(
    boundary_crs4326: geopandas.GeoDataFrame,
    element_type: str,
    selector: str,
    element_dict: {} = None,
):

    # Construct query for islands in bbox
    query = OSMPythonTools.overpass.overpassQueryBuilder(
        bbox=[
            boundary_crs4326.bounds.miny[0],
            boundary_crs4326.bounds.minx[0],
            boundary_crs4326.bounds.maxy[0],
            boundary_crs4326.bounds.maxx[0],
        ],
        elementType=element_type,
        selector=selector,
        out="body",
        includeGeometry=True,
    )

    # Perform query
    overpass = OSMPythonTools.overpass.Overpass()
    islands = overpass.query(query)

    # Create dictionary for results if not passed in
    if element_dict is None:
        element_dict = {
            "geometry": [],
            "OSM_id": [],
            "name": [],
            "type": [],
        }
    # Add results to the dictionary
    for element in islands.elements():
        tags = element.tags()
        element_dict["geometry"].append(element.geometry())
        element_dict["OSM_id"].append(element.id())
        element_dict["name"].append(
            tags["name"] if "name" in tags.keys() else "Unnamed"
        )
        element_dict["type"].append(
            tags["place"] if "place" in tags.keys() else "Not Specified"
        )
    return element_dict

def get_gadm_islands_in_boundary(island_dict: dict, all_land_path: pathlib.Path, output_path: pathlib.Path):
    """ All islands within the boundary - clip to boundary. """
    
    boundary_crs4326 = create_boundary_epsg4326(island_dict=island_dict, output_path=output_path)
    island_dict["boundary"] = boundary_crs4326
    local_crs = island_dict["crs"]
    
    land_boundary = geopandas.read_file(all_land_path).to_crs(local_crs)
    land_boundary.clip(boundary_crs4326.to_crs(local_crs), keep_geom_type=True)
    land_boundary = land_boundary.dissolve()
    label = "gadm"
    
    if len(land_boundary) == 0:
        land_boundary = get_osm_islands_in_boundary(
            boundary_crs4326=boundary_crs4326, local_crs=local_crs
        )
        land_boundary = land_boundary.dissolve()
        label = "osm"
    
    land_boundary.to_file(output_path / f"{label}_land_{island_dict['name']}.geojson")

    return land_boundary
    
    

def get_osm_islands_in_boundary(boundary_crs4326: geopandas, local_crs: str):
    """ All islands within the boundary - clip to boundary. Return islands defined by
    ways and relations. """
    # Query 'ways'
    island_dict = osm_query_islands(
        boundary_crs4326=boundary_crs4326,
        element_type="way",
        selector='"place"="island"',  #'"natural"="coastline"', <- using place instead of coastline as is always a polygon
    )
    island_dict = osm_query_islands(
        boundary_crs4326=boundary_crs4326,
        element_type="way",
        selector='"place"="islet"',
        element_dict=island_dict,
    )

    # Query "Relations"
    island_dict = osm_query_islands(
        boundary_crs4326=boundary_crs4326,
        element_type="relation",
        selector='"place"="island"',
        element_dict=island_dict,
    )

    island_dict = osm_query_islands(
        boundary_crs4326=boundary_crs4326,
        element_type="relation",
        selector='"place"="islet"',
        element_dict=island_dict,
    )

    # Query "Relations" reef
    reef_dict = osm_query_islands(
        boundary_crs4326=boundary_crs4326,
        element_type="relation",
        selector='"natural"="reef"',
    )
    # Query "way" reef
    reef_dict = osm_query_islands(
        boundary_crs4326=boundary_crs4326,
        element_type="way",
        selector='"natural"="reef"',
        element_dict=reef_dict,
    )
    # Create data frames
    islands = []
    if len(island_dict['geometry']) == 0 and len(reef_dict['geometry']) == 0:
        print("No island, islets or reef's found in the specified bounding box")
    if len(island_dict['geometry']) > 0:
        islands.append(geopandas.GeoDataFrame(island_dict, crs=OSM_CRS))
    # Add reefs if any exist - tidy up first
    if len(reef_dict['geometry']) > 0:
        reefs = geopandas.GeoDataFrame(reef_dict, crs=OSM_CRS)
        reefs = reefs[(reefs.geometry.type == "Polygon") | (reefs.geometry.type == "MultiPolygon")]
        # Tidy up reef polygons - split MultiPolygons
        reefs["geometry"] = reefs.apply(
            lambda row: shapely.geometry.MultiPolygon(
                shapely.geometry.Polygon(p.exterior) for p in row.geometry.geoms
            )
            if row.geometry.type == "MultiPolygon"
            else row.geometry,
            axis=1,
        )
        reefs["type"] = "reef"
        reefs = reefs.explode(index_parts=False)
        # Take the outside of Polygons
        reefs["geometry"] = reefs.apply(
            lambda row: shapely.geometry.Polygon(row.geometry.exterior), axis=1,
        )
        islands.append(reefs)
    # Combine
    islands = geopandas.GeoDataFrame(
        pandas.concat(islands, ignore_index=True), crs=OSM_CRS
    )

    # Combine any duplicate OSM ids
    islands = islands.dissolve(by="OSM_id")
    islands = islands.clip(boundary_crs4326).to_crs(local_crs)
    return islands


def construct_fab_folder_name(lat: int, lon: int, fab_version: str):
    "Constuct the FABDEM folder name from the lat and lon coordinates."

    if lon <= 0:
        east_west = [
            f"W{abs(int(numpy.floor((lon + 0.5) / 10.0)) * 10):02d}",
            f"W{abs(int(numpy.ceil((lon + 0.5) / 10.0)) * 10):02d}",
        ]
    elif lon >= 170:
        east_west = [
            f"E{int(numpy.floor((lon + 0.5) / 10.0)) * 10:02d}",
            f"W{int(numpy.ceil((lon + 0.5) / 10.0)) * 10:02d}",
        ]
    else:
        east_west = [
            f"E{int(numpy.floor((lon + 0.5) / 10.0)) * 10:02d}",
            f"E{int(numpy.ceil((lon + 0.5) / 10.0)) * 10:02d}",
        ]
    if lat >= 0:
        north_south = [
            f"N{int(numpy.floor((lat + 0.5) / 10.0)) * 10:02d}",
            f"N{int(numpy.ceil((lat + 0.5) / 10.0)) * 10:02d}",
        ]
    elif lat >= -10:
        north_south = [
            f"S{abs(int(numpy.floor((lat + 0.5) / 10.0)) * 10):02d}",
            f"N{abs(int(numpy.ceil((lat + 0.5) / 10.0)) * 10):02d}",
        ]
    else:
        north_south = [
            f"S{abs(int(numpy.floor((lat + 0.5) / 10.0)) * 10):02d}",
            f"S{abs(int(numpy.ceil((lat + 0.5) / 10.0)) * 10):02d}",
        ]
    folder_name = (
        f"{north_south[0]}{east_west[0]}-" f"{north_south[1]}{east_west[1]}_FABDEM_{fab_version}"
    )

    return folder_name


def construct_fab_name(lat: int, lon: int, fab_version: str):
    "Constuct the FABDEM folder name from the lat and lon coordinates."

    if lon <= 0:
        east_west = "W"
    else:
        east_west = "E"
    if lat >= 0:
        north_south = "N"
    else:
        north_south = "S"
    name = f"{north_south}{abs(lat):02d}" f"{east_west}{abs(lon):02d}_FABDEM_{fab_version}.tif"

    return name


def combine_fabs_in_boundary(
    island_dict: dict, fab_path: pathlib, fab_version: str = "V1-2"
):
    """ Combine all FAB lat and lon tiles into a single DEM. """

    bounds = island_dict["boundary"].bounds
    lats = numpy.arange(
        numpy.floor(bounds.miny[0]), numpy.floor(bounds.maxy[0]) + 1, 1, dtype=int
    )
    lons = numpy.arange(
        numpy.floor(bounds.minx[0]), numpy.floor(bounds.maxx[0]) + 1, 1, dtype=int
    )

    fabs = []
    for lat in lats:
        for lon in lons:
            fab_folder_name = construct_fab_folder_name(lat=lat, lon=lon, fab_version=fab_version)
            fab_name = construct_fab_name(lat=lat, lon=lon, fab_version=fab_version)
            if (fab_path / fab_version / fab_folder_name).exists():
                if (fab_path / fab_version / fab_folder_name / fab_name).exists():
                    fabs.append(
                        load_dem(fab_path / fab_version / fab_folder_name / fab_name)
                    )
                else:
                    print(f"/tWarning FABDEM tile {fab_name} doesn't exist in folder {fab_folder_name}")
            elif (fab_path / fab_version / f"{fab_folder_name}.zip").exists(): 
                # Not yet unzipped
                shutil.unpack_archive(str(fab_path / fab_version / f"{fab_folder_name}.zip"),
                                      str(fab_path / fab_version / fab_folder_name), 'zip')
                if (fab_path / fab_version / fab_folder_name / fab_name).exists():
                    fabs.append(
                        load_dem(fab_path / fab_version / fab_folder_name / fab_name)
                    )
                else:
                    print(f"/tWarning FABDEM tile {fab_name} doesn't exist in folder {fab_folder_name}")
            else:
                print(f"No {fab_folder_name}/{fab_name} FABDEM exists")
    fab = rioxarray.merge.merge_arrays(fabs)
    fab = (
        fab.rio.clip(island_dict["boundary"].geometry)
        .rio.reproject(island_dict["crs"])
        .rio.write_crs(island_dict["crs"])
    )
    fab.rio.write_nodata(fab.rio.nodata).to_netcdf(
        fab_path
        / f"{island_dict['name']}_FABDEM_crs{island_dict['crs']}.nc"
    )
    return fab


def create_boundary_epsg4326(island_dict: dict, output_path: pathlib.Path):
    """ Create a GeoDataFrame of the boundary of an island group from lats and lons. """

    boundary_crs4326 = geopandas.GeoDataFrame(
        geometry=[
            shapely.geometry.Polygon(zip(island_dict["lons"], island_dict["lats"]))
        ],
        crs=OSM_CRS,
    )
    boundary_crs4326.to_file(
        output_path / f"boundary_{island_dict['name']}.geojson"
    )
    return boundary_crs4326


def resample_lidar(lidar_dem: xarray.DataArray, resolution: float, clipping_boundary: geopandas.GeoDataFrame) -> xarray.DataArray:
    if clipping_boundary is not None:
        lidar_dem = lidar_dem.rio.clip(
            clipping_boundary.buffer(resolution).geometry,
            drop=True, all_touched=False
        )
    else:
        clipping_boundary = get_bbox_from_raster(lidar_dem)
    resampled_lidar = resample_dem(
        dem_in=lidar_dem,
        resolution=resolution,
        boundary=clipping_boundary.buffer(resolution),
    )
    return resampled_lidar

def dem_data_extents(dem: xarray.DataArray):
    resolution = max(dem.rio.resolution())
    
    extents = [shapely.geometry.shape(polygon[0]) for polygon in 
               rasterio.features.shapes(numpy.uint8(dem.notnull().values)) if polygon[1] == 1.0 ]
    extents = shapely.ops.unary_union(extents)

    # Remove internal holes for select types as these may cause self-intersections
    if type(extents) is shapely.geometry.Polygon:
        extents = shapely.geometry.Polygon(extents.exterior)
    elif type(extents) is shapely.geometry.MultiPolygon:
        extents = shapely.geometry.MultiPolygon( [shapely.geometry.Polygon(polygon.exterior) for polygon in extents.geoms ] )
    # Convert into a Geopandas dataframe
    extents = geopandas.GeoDataFrame( {"geometry": [extents]}, crs=dem.rio.crs, )

    # Move from image to the dem space & buffer(0) to reduce self-intersections
    transform = dem.rio.transform()
    extents = extents.affine_transform(
        [ transform.a, transform.b, transform.d, transform.e, transform.xoff, transform.yoff, ]
    ).buffer(0)

    # And make our GeoSeries into a GeoDataFrame
    extents = geopandas.GeoDataFrame(geometry=extents)
    
    # Buffer by 2x resolution and dissolve
    extents = geopandas.GeoDataFrame(geometry=extents.buffer(resolution)).dissolve().buffer(-resolution)
    return extents

def create_datasource_layer(fab_dem: xarray.DataArray, lidar_dem: xarray.DataArray,
                            clipping_boundary: geopandas.GeoDataFrame = None):
    # Create, combine and save data layer
    dems = []
    if lidar_dem is not None:
        data_source_lidar = xarray.zeros_like(lidar_dem, dtype=numpy.int32)
        data_source_lidar.rio.set_nodata(0, inplace=True)
        data_source_lidar.rio.write_nodata(data_source_lidar.rio.nodata, inplace=True)
        data_source_lidar = data_source_lidar.where(lidar_dem.isnull().data, 1)
        dems.append(data_source_lidar)
    if fab_dem is not None:
        data_source_fab = xarray.zeros_like(fab_dem, dtype=numpy.int32)
        data_source_fab.rio.set_nodata(0, inplace=True)
        data_source_fab.rio.write_nodata(data_source_fab.rio.nodata, inplace=True)
        data_source_fab = data_source_fab.where(fab_dem.isnull(), 2)
        dems.append(data_source_fab)
    if len(dems) > 1:
        data_source = rioxarray.merge.merge_arrays(dems, method="first",)
    else:
        data_source = dems[0]
    if clipping_boundary is not None:
        data_source = data_source.rio.clip(clipping_boundary.geometry, drop=True, all_touched=True)
    data_source.rio.write_crs(data_source.rio.crs, inplace=True)
    data_source = data_source.assign_attrs(title="The data source for each pixel", description="0 if no data, 1 if LiDAR, 2 if FABDEM")
    for attribute in ["AREA_OR_POINT", "CHANGELOG", "INSTITUTE"]:
        if attribute in data_source.attrs:
            data_source.attrs.pop(attribute)
    return data_source

def get_bbox_from_geometry(extents: geopandas.GeoDataFrame):
    bounds = extents.bounds
    bbox = geopandas.GeoDataFrame(geometry=[shapely.geometry.Polygon([[bounds.minx, bounds.miny], [bounds.minx, bounds.maxy],
                                                           [bounds.maxx, bounds.maxy], [bounds.maxx, bounds.miny]])],
                                 crs=extents.crs)
    return bbox

def get_bbox_from_raster(raster: xarray.DataArray):
    bounds = raster.rio.bounds()
    bbox = geopandas.GeoDataFrame(
        geometry=[shapely.geometry.Polygon([[bounds[0], bounds[1]], [bounds[2], bounds[1]],
                                            [bounds[2], bounds[3]], [bounds[0], bounds[3]]])],
        crs=raster.rio.crs,
    )

    return bbox

def creating_dems_all_islands(
    island_groups: dict, resolution: float, output_path: pathlib.Path
):
    """ Loop through each island group and save out the LiDAR DEM resampled to a coarser resolution.
    No triming is applied to the LiDAR. """
    
    for island_name, island_values in island_groups.items():
        print(f"Producing DEM(s) for {island_name} at {resolution}m.")
        label = ""
        land = island_values["land"]
        
        if "lidar" in island_values:
            # Combine all LiDAR DEMs
            if type(island_values["lidar"]) is not list:
                island_values["lidar"] = [island_values["lidar"]]
                
            lidar_list_5m = []
            for index, lidar_dem in enumerate(island_values["lidar"]):    
                print(f"\tReample LiDAR DEM {index + 1} of {len(island_values['lidar'])}")
                lidar_5m = resample_lidar(lidar_dem=lidar_dem,
                                          clipping_boundary=None,
                                          resolution=resolution)
                lidar_list_5m.append(lidar_5m)
            if len(lidar_list_5m) == 1:
                lidar_5m = lidar_list_5m[0]
            else:
                print(f"\tCombining {len(lidar_list_5m)} DEMs.")
                lidar_5m = rioxarray.merge.merge_arrays(lidar_list_5m, method="first")

                print("\tDefine extents of LiDAR DEMs.")
                lidar_extents = dem_data_extents(lidar_5m)
                lidar_extents.to_file(output_path / f"{resolution}m_dem_{island_name}_lidar_only.geojson")

                print(f"\tInterpolate small gaps in the LiDAR DEMs.")
                lidar_5m = lidar_5m.rio.interpolate_na(method="nearest")
                lidar_5m = lidar_5m.rio.clip(lidar_extents.geometry, drop=True, all_touched=False)
        else:
            lidar_5m = None
            
        # If not only LiDAR bring in FABDEM
        if "lidar_only" in island_values:
            fab_5m = None
            label += "_lidar_only"
        else:
            # Trim FAB then resample
            trimmed_fab = island_values["fab"].rio.clip(
                land.buffer(max(island_values["fab"].rio.resolution())).geometry,
                all_touched=True)
            fab_5m = resample_dem(
                dem_in=trimmed_fab,
                resolution=resolution,
                boundary=land,
            )
        if lidar_5m is not None and fab_5m is not None:
            print("\tMergining FAB and LiDAR DEMs.")
            label = "_including_lidar"
            merged_dem_5m = rioxarray.merge.merge_arrays([lidar_5m, fab_5m], method="first",)
        elif lidar_5m is not None:
            merged_dem_5m = lidar_5m
        elif fab_5m is not None:
            merged_dem_5m = fab_5m
        else:
            print("\tNo DEMs. Exiting.")
            return
        
        print("\tSaving DEMs.")
        merged_dem_5m.to_dataset(name="dem").to_netcdf(output_path / f"{resolution}m_dem_{island_name}{label}.nc",
                                                       encoding={"dem":  {'zlib': True, 'complevel': 1, }})
        merged_dem_5m.rio.to_raster(output_path / f"{resolution}m_dem_{island_name}{label}.tif", compress='deflate') #'zlib', 'deflate', "lzw"
        
        # Create, combine and save data layer
        data_source = create_datasource_layer(fab_dem=fab_5m, lidar_dem=lidar_5m, clipping_boundary=None)
        data_source.rio.to_raster(output_path / f"{resolution}m_dem_data_source_{island_name}{label}.tif", compress='deflate')
        
        if lidar_5m is not None:
            print("\tClip then save DEMs to land.")
            merged_dem_5m = merged_dem_5m.rio.clip(land.geometry, all_touched=False)
            merged_dem_5m.to_dataset(name="dem").to_netcdf(output_path / f"{resolution}m_dem_{island_name}{label}_land.nc",
                                                           encoding={"dem":  {'zlib': True, 'complevel': 1, }})
            merged_dem_5m.rio.to_raster(output_path / f"{resolution}m_dem_{island_name}{label}_land.tif", compress='deflate') #'zlib', 'deflate', "lzw"
            data_source = create_datasource_layer(fab_dem=fab_5m, lidar_dem=lidar_5m, clipping_boundary=land)
            data_source.rio.to_raster(output_path / f"{resolution}m_dem_data_source_{island_name}{label}_land.tif", compress='deflate')