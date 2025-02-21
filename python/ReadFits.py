################################################################################
#
# Author: Frank Soldano
#
# Read fits file and print some information 
#
#
################################################################################

from astropy.io import fits
from astropy.table import Table
import argparse
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import ZScaleInterval, ImageNormalize
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astropy.wcs import WCS
from tqdm import trange
from Gaia_tools import gaia_catalog_names_coord
from General_tools import init_dict_stat

parser = argparse.ArgumentParser(description='Read a Fits file')
parser.add_argument('file',type=str, help='File to read')
parser.add_argument('--gaiaCat',type=str,help='Name of the gaia catalog to read')
parser.add_argument('--show',help='show plot',action='store_true')
args = parser.parse_args()

file = args.file

fitsFile = fits.open(file)
img_data = fitsFile[0].data
img_header = fitsFile[0].header
wcs = WCS(img_header) 

coord = SkyCoord(img_header['CRVAL1']*u.deg,img_header['CRVAL2']*u.deg)

coneGaia = u.Quantity(3, u.deg)
if not args.gaiaCat:
    table_name = gaia_catalog_names_coord(coord,coneGaia)
    print('Query gaia stars...')
    Gaia.ROW_LIMIT = -1
    Gaia.cone_search_async(coordinate=coord, radius=coneGaia, output_file=table_name,output_format='fits',dump_to_file=True)

else:
    table_name = args.gaiaCat

gaiaTable = fits.open(table_name)[1].data
print(f'Table {table_name} read!')

gaia_stars = SkyCoord(gaiaTable['ra']*u.deg,gaiaTable['dec']*u.deg)
gaia_stars_pixel = wcs.world_to_pixel(gaia_stars)

starsImage = np.where( (0<gaia_stars_pixel[1]) & (gaia_stars_pixel[1]<img_data.shape[0]) & 
                       (0<gaia_stars_pixel[0]) & (gaia_stars_pixel[0]<img_data.shape[1]) & 
                       (gaiaTable['phot_g_mean_mag'] < 16))

x_gaia = gaia_stars_pixel[0][starsImage]
y_gaia = gaia_stars_pixel[1][starsImage]

print('Number of stars found in img:',len(x_gaia))

norm = ImageNormalize(img_data, interval=ZScaleInterval())
plt.figure(figsize=(15,15))
plt.imshow(img_data,norm=norm,origin='lower',cmap='gray')
plt.plot(x_gaia,y_gaia,'r.')
plt.savefig(file+'.png')
if args.show:
    plt.show()
plt.close()

# Extract catalogue with only stars of the image
# Column from the Gaia cat
column_name_Gaia = ['designation','source_id',
                    'ra','dec','parallax','pmra','pmdec',
                    'phot_bp_mean_mag','phot_rp_mean_mag',
                    'phot_g_mean_mag','bp_rp','bp_g','g_rp',
                    'phot_variable_flag']

gaia_export_cat = init_dict_stat(column_name_Gaia)

for column in gaia_export_cat :
    gaia_export_cat[column] = gaiaTable[column][starsImage]

# Add column from this code
column_name_img =['X_IMAGE','Y_IMAGE']

img_cat = init_dict_stat(column_name_img)
img_cat['X_IMAGE'] = x_gaia
img_cat['Y_IMAGE'] = y_gaia

final_table = Table(gaia_export_cat | img_cat)
final_table.write(file+'_cat.fits')