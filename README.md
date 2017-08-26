#Void - tSZ cross correlation

This repository contains the full analysis pipeline used in INSERT_ARXIV_ID_HERE.

To reproduce the results in the paper:
1. Compile the stacking code by typing 'make' (this may require modifying Makefile to include the correct paths to the FITS and HEALPix libraries).
2. Run ./pipeline.sh. This script will run the full pipeline, including downloading all relevant datasets, computing signal and covariance matrix, estimating the theory prediction, fitting it to the data and generating all the plots in the paper. Some of the commands (e.g. to launch parallel jobs) in this script are platform specific.


