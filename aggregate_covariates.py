"""Process global datasets related to Nitrate.

The purpose here is to aggregate spatial data from various global datasets into
scalar predictor variables to be included in an analysis of the relationship
between catchment characteristics like landcover, topography, and
precipitation, and nitrite/nitrate loadings in surface water.

Major tasks:
 - summarize covariate datasets within catchment basins

"""
import os

from osgeo import gdal
from osgeo import osr
from osgeo import ogr
import pygeoprocessing
import pygeoprocessing.routing
# import taskgraph

# shapefile containing locations of stations with NOxN observations
_STATION_SHP_PATH = "C:/Users/ginge/Documents/NatCap/GIS_local/NCI_NDR/WB_surface_stations_noxn_obs_sa_bas.shp"

# shapefile containing watershed area corresponding to WB stations
_WB_BASIN_SHP_PATH = "F:/NCI_NDR/Watersheds_DRT/Rafa_watersheds_v2/WB_surface_stations_noxn_for_snapping_ContribArea.shp"


def map_FID_to_field(shp_path, field):
    """Map FID of each feature, according to GetFID(), to the given field.

    This allows for mapping of a dictionary of zonal statistics, where keys
    correspond to FID according to GetFID(), to another field that is preferred
    to identify features.

    Parameters:
        shp_path (string): path to shapefile
        field (string): the field to map to FID

    Returns:
        dictionary indexed by the FID of each feature retrieved with GetFID(),
            and values are the value of `field` for the feature

    """
    vector = gdal.OpenEx(shp_path, gdal.OF_VECTOR)
    layer = vector.GetLayer()
    FID_to_field = {
        feature.GetFID(): feature.GetField(field) for feature in layer}

    # clean up
    vector = None
    layer = None
    return FID_to_field


def main():
    """Program entry point."""



if __name__ == '__main__':
    main()
