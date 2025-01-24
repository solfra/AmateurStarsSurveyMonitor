################################################################################
#
#           Autor: Frank Soldano
# 
#  This is a librairy for works with Gaia catalogs
# 
################################################################################

from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astropy.wcs import WCS

def gaia_catalog_names_coord(coords,serach_radius):
    """
    Returne the name of the catalogue based on the coordinate uery and the cone radius.
    Must be the catalogue name for Gaia.cone_search_async

    Input:
        - coord: A SkyCoord object
        - radius: the radius used for cone search, in Quantity units
    
    Return:
        - cat_name: Catalogue name
    """
    ra  = coords.ra.to_string(alwayssign=True)
    dec = coords.dec.to_string(alwayssign=True)
    
    radius = serach_radius.to_string().split(' ')
    radius_value = radius[0]
    radius_unit  = radius[1]

    cat_name = f'GaiaCat_ra_{ra}_dec_{dec}_in_{radius_value}_{radius_unit}.fits.gz'

    return cat_name

