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

parser = argparse.ArgumentParser(description='Read a Fits file')
parser.add_argument('file',type=str, help='File to read')
args = parser.parse_args()

file = args.file

fitsFile = fits.open(file)
for key in fitsFile[0].header:
    print(f' {key}: {fitsFile[0].header[key]}')
