# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 14:50:20 2022F

@author: pearsonra
"""

import numpy
import xarray
import rioxarray
import rioxarray.merge
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
            if (fab_path / fab_version / fab_folder_name / fab_name).exists():
                fabs.append(
                    load_dem(fab_path / fab_version / fab_folder_name / fab_name)
                )
            elif (fab_path / fab_version / f"{fab_folder_name}.zip").exists(): 
                # Not yet unzipped
                shutil.unpack_archive(str(fab_path / fab_version / f"{fab_folder_name}.zip"),
                                      str(fab_path / fab_version / fab_folder_name), 'zip')
                fabs.append(
                    load_dem(fab_path / fab_version / fab_folder_name / fab_name)
                )
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


def loop_through_islands_creating_dems(
    island_groups: dict, resolution: float, output_path: pathlib.Path, buffer: bool = False
):
    """ Loop through each island group and merge then save each out. If buffer is set the DEMs
    are created around the islands at a distance of 15 x resolution. """
    for island_name, island_values in island_groups.items():
        print(f"Combining DEMs for {island_name}")
        islands = get_osm_islands_in_boundary(
            boundary_crs4326=island_values["boundary"], local_crs=island_values["crs"]
        )
        islands.to_file(output_path / f"osm_islands_{island_name}.geojson")
        islands = islands.dissolve()
        # Add a buffer (i.e. 75m say) around each island before clipping
        clipping_boundary = islands.buffer(resolution * 15) if buffer else islands

        # Trim FAB then resample
        trimmed_fab = island_values["fab"].rio.clip(
            clipping_boundary.buffer(max(island_values["fab"].rio.resolution())).geometry,
            all_touched=True)
        fab_5m = resample_dem(
            dem_in=trimmed_fab,
            resolution=resolution,
            boundary=clipping_boundary,
        )

        # Trim LiDAR DEM - if it exists
        if "lidar" in island_values.keys():
            trimmed_lidar = island_values["lidar"].rio.clip(
                clipping_boundary.buffer(resolution).geometry,
                drop=True, all_touched=True
            )
            lidar_5m = resample_dem(
                dem_in=trimmed_lidar,
                resolution=resolution,
                boundary=clipping_boundary.buffer(resolution),
            )

            # Combine
            merged_dem_5m = rioxarray.merge.merge_arrays(
                [lidar_5m, fab_5m], method="first",
            )
            
        else:
            merged_dem_5m = fab_5m

        # Save files
        merged_dem_5m = merged_dem_5m.rio.clip(clipping_boundary.geometry,
                                               drop=True, all_touched=True
                                              )
        merged_dem_5m.to_dataset(name="dem").to_netcdf(output_path / f"{resolution}m_dem_{island_name}.nc",
                                                       encoding={"dem":  {'zlib': True, 'complevel': 1, }})

        merged_dem_5m.rio.to_raster(output_path / f"{resolution}m_dem_{island_name}.tif", compress='deflate')

        
        # Create and save data layer
        data_source = xarray.zeros_like(merged_dem_5m, dtype=numpy.int32)
        data_source.rio.set_nodata(-1, inplace=True)
        data_source.rio.write_nodata(data_source.rio.nodata, inplace=True)
        fab_5m = fab_5m.rio.clip(clipping_boundary.geometry, drop=True, all_touched=True)
        data_source = data_source.where(fab_5m.isnull(), 2)
        if "lidar" in island_values.keys():
            lidar_5m = lidar_5m.rio.clip(clipping_boundary.geometry, drop=True, all_touched=True)
            data_source = data_source.where(lidar_5m.isnull().data, 1)
        data_source.rio.write_crs(data_source.rio.crs, inplace=True)
        data_source = data_source.assign_attrs(title="The data source for each pixel", description="0 if no data, 1 if LiDAR, 2 if FABDEM")
        for attribute in ["AREA_OR_POINT", "CHANGELOG", "INSTITUTE"]:
            if attribute in data_source.attrs:
                merged_dem_5m.attrs.pop(attribute)
        data_source.rio.to_raster(output_path / f"{resolution}m_dem_data_source_{island_name}.tif", compress='deflate')
        
def loop_through_islands_creating_dems_from_lidar_only(
    island_groups: dict, resolution: float, output_path: pathlib.Path
):
    """ Loop through each island group and save out the LiDAR DEM resampled to a coarser resolution.
    No triming is applied to the LiDAR. """
    
    for island_name, island in island_groups.items():
        print(f'Processing Island {island_name}')

        bounds = island.rio.bounds()
        boundary = geopandas.GeoDataFrame(
            geometry=[
                shapely.geometry.Polygon([[bounds[0], bounds[1]], [bounds[2], bounds[1]],
                                          [bounds[2], bounds[3]], [bounds[0], bounds[3]]])
            ],
            crs=island.rio.crs,
        )
        dem_5m = resample_dem(
                    dem_in=island,
                    resolution=resolution,
                    boundary=boundary,
                )
        dem_5m.to_dataset(name="dem").to_netcdf(output_path / f"{resolution}m_dem_{island_name}.nc",
                                                encoding={"dem":  {'zlib': True, 'complevel': 1, }})
        dem_5m.rio.to_raster(output_path / f"{resolution}m_dem_{island_name}.tif", compress='deflate') #'zlib', 'deflate', "lzw"
