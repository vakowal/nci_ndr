"""Process global datasets related to Nitrate.

The purpose here is to aggregate spatial data from various global datasets into
scalar predictor variables to be included in an analysis of the relationship
between catchment characteristics like landcover, topography, and
precipitation, and nitrite/nitrate loadings in surface water.

Major tasks:
 - delineate catchment basins for each surface water observation point
 - summarize covariate datasets within catchment basins

"""
from osgeo import gdal
from osgeo import ogr
import pygeoprocessing
import pygeoprocessing.routing
import taskgraph

# shapefile containing locations of stations with NOxN observations
_STATION_SHP_PATH = "F:/NCI_NDR/Data worldbank/station_data/station_df_noxn_surface_gte_1990.shp"
_FLOW_DIR_PATH = "F:/NCI_NDR/Data flow direction DRT/globe_fdr_shedsandh1k.asc"


def delineate_basins(outlet_shapefile_path):
    """Delineate watersheds from points representing watershed outlets.

    Steal code from delineateit.py to hack a custom watershed delineation
    routine.
    The goal is, for each point in the outlet shapefile, to identify pixels
    that contribute topographically to the pixel containing that point.
    Use the flow direction raster used by Barbarossa et al (CITE) so that we
    can as much as possible match the basins contributing to streamflow
    estimates in the FLO1K dataset.

    """


def main():
    """Program entry point."""
    delineate_basins(_STATION_SHP)


if __name__ == '__main__':
    main()
