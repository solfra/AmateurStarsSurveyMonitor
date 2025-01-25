################################################################################
#
# Author: Frank Soldano
#
# Read fits file and print some information 
#
#
################################################################################

from astropy.io import fits
import argparse
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval, ImageNormalize
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astropy.wcs import WCS
from tqdm import trange
from Gaia_tools import gaia_catalog_names_coord

parser = argparse.ArgumentParser(description='Read a Fits file')
parser.add_argument('file',type=str, help='File to read')
args = parser.parse_args()

file = args.file

fitsFile = fits.open(file)
img_data = fitsFile[0].data
img_header = fitsFile[0].header

coord = SkyCoord(img_header['CRVAL1']*u.deg,img_header['CRVAL2']*u.deg)

coneGaia = u.Quantity(3, u.deg)
table_name = gaia_catalog_names_coord(coord,coneGaia)

print('Query gaia stars...')
Gaia.ROW_LIMIT = -1
Gaia.cone_search_async(coordinate=coord, radius=coneGaia, output_file=table_name,output_format='fits',dump_to_file=True)

wcs = WCS(img_header) 

gaiaTable = fits.open(table_name)[1].data

x_gaia, y_gaia = [], []

for i in trange(len(gaiaTable['source_id'])): 
    stars = wcs.world_to_pixel(SkyCoord(gaiaTable['ra'][i]*u.deg,gaiaTable['dec'][i]*u.deg))

    if 0 < stars[1] < img_data.shape[0] and 0 < stars[0] < img_data.shape[1] and gaiaTable['phot_g_mean_mag'][i] < 16:
        x_gaia.append(stars[0])
        y_gaia.append(stars[1])
print('Number of stars found in img:',len(x_gaia))

norm = ImageNormalize(img_data, interval=ZScaleInterval())
plt.figure(figsize=(15,15))
plt.imshow(img_data,norm=norm,origin='lower',cmap='gray')
plt.plot(x_gaia,y_gaia,'r.')
plt.savefig(file+'.png')