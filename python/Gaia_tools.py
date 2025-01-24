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
    
    Output:
        - cat_name: Catalogue name
    """
    ra  = coords.ra.to_string(alwayssign=True)
    dec = coords.dec.to_string(alwayssign=True)
    
    radius = serach_radius.to_string().split(' ')
    radius_value = radius[0]
    radius_unit  = radius[1]

    cat_name = f'GaiaCat_ra_{ra}_dec_{dec}_in_{radius_value}_{radius_unit}.fits.gz'

    return cat_name

def coord_from_cat_name(cat_name):
    """
    Extract coordinate from the name of the Gaia catalogue

    Input: 
        - cat_name: the catalogue name, in the format of the results of gaia_catalog_names_coord()

    Output:
        - coord: A skyCoord object 
        - search_radius: the radius used for cone search, in Quantity units
    """

    cat_name_list = cat_name.split('_')
    ra  = cat_name_list[cat_name_list.index('ra') +1]
    dec = cat_name_list[cat_name_list.index('dec')+1]

    radius_value = cat_name_list[cat_name_list.index('in') +1]
    radius_unit  = cat_name_list[cat_name_list.index('in') +2]

    coord = SkyCoord(ra,dec)
    search_radius = u.Quantity(radius_value, u.Unit(radius_unit))
    
    return coord, search_radius