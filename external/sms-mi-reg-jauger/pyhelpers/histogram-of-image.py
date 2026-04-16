#!/usr/bin/env python3

#
# pip3 install SimpleITK
# pip3 install NumPy
# pip3 install matplotlib
# sudo pip3 install scipy
#

import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import sys
import argparse

parser = argparse.ArgumentParser(description="Display the histogram of an image.")
parser.add_argument("--input", required=True)
parser.add_argument("--outputHistogram", required=True)
args = parser.parse_args()

imageFileName = args.input
outputPlotFileName = args.outputHistogram

reader = sitk.ImageFileReader()
reader.SetFileName( imageFileName )
image = reader.Execute();
 
nda = sitk.GetArrayFromImage( image )

x = nda.flatten()
max_ele = np.amax( x )
min_ele = np.amin( x )
print('Image array min is ' + str(min_ele) + ', max is : ' + str(max_ele))

nphistbins, bin_edges = np.histogram( x, bins = 64 )

print('Histogram is : ')
print(nphistbins)

print('Size of histogram is : ' + str( len( nphistbins )) )
print('Shape of histogram is : ')
print( nphistbins.shape )

trimmed_histogram = nphistbins[1:-1]
print(trimmed_histogram)
xaxis = np.linspace(bin_edges[0], bin_edges[-1], len(trimmed_histogram))
plt.plot(xaxis, trimmed_histogram)
plt.savefig( outputPlotFileName )
plt.show()

exit(1)

# density = False plots counts
# density = True plots a density.
plt.hist(x, density=False, bins=64, label="Data") 
mn, mx = plt.xlim()
print('Limits of density are : ' + str(mn) + ', ' + str(mx) )
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = st.gaussian_kde(x)
plt.plot(kde_xs, kde.pdf(kde_xs), label="Probability Density Function" )
plt.legend(loc="upper left")
plt.ylabel('Probability')
plt.xlabel('Data')
plt.title('Histogram')
plt.savefig( outputPlotFileName )
plt.show()

quit(0)

