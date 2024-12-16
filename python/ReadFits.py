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

parser = argparse.ArgumentParser(description='Read a Fits file')
parser.add_argument('file',type=str, help='File to read')
args = parser.parse_args()

file = args.file

fitsFile = fits.open(file)
for key in fitsFile[0].header:
    print(f' {key}: {fitsFile[0].header[key]}')

img_data = fitsFile[0].data
norm = ImageNormalize(img_data, interval=ZScaleInterval())
plt.figure(figsize=(15,15))
plt.imshow(img_data,norm=norm,origin='lower',cmap='gray')
plt.show()