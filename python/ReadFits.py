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

parser = argparse.ArgumentParser(description='Read a Fits file')
parser.add_argument('file',type=str, help='File to read')
args = parser.parse_args()

file = args.file

fitsFile = fits.open(file)
img_data = fitsFile[0].data
img_header = fitsFile[0].header

for key in img_header:
    print(f' {key}: {img_header[key]}')

norm = ImageNormalize(img_data, interval=ZScaleInterval())
plt.figure(figsize=(15,15))
plt.imshow(img_data,norm=norm,origin='lower',cmap='gray')
plt.show()

coord = SkyCoord(img_header['CRVAL1']*u.deg,img_header['CRVAL2']*u.deg)

coneGaia = u.Quantity(5, u.deg)

Gaia.cone_search_async(coordinate=coord, radius=coneGaia, output_file=f'{file}_gaia_cone.fits',output_format='fits',dump_to_file=True)