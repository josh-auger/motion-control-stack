#!/usr/bin/env python3

# pip3 install SimpleITK
# pip3 install NumPy
# pip3 install matplotlib
# sudo pip3 install scipy

import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import sys

if len(sys.argv) < 4:
  print("Usage: " + sys.argv[0] + " image1FileName image2FileName outPlot")
  quit(1)

image1FileName = sys.argv[1]
image2FileName = sys.argv[2]
outputPlotFileName = sys.argv[3]

print(f"Arguments count: {len(sys.argv)}")
for i, arg in enumerate(sys.argv):
    print(f"Argument {i:>6}: {arg}")

reader1 = sitk.ImageFileReader()
reader1.SetFileName( image1FileName )
image1 = reader1.Execute();
nda1 = sitk.GetArrayFromImage( image1 )
print('nda1 min, max: ' + str(np.min(nda1)) + ', ' + str(np.max(nda1)) )

reader2 = sitk.ImageFileReader()
reader2.SetFileName( image2FileName )
image2 = reader2.Execute();
nda2 = sitk.GetArrayFromImage( image2 )
print('nda2 min, max: ' + str(np.min(nda2)) + ', ' + str(np.max(nda2)) )

x1 = nda1.flatten()
x2 = nda2.flatten()
# density = False would plot counts
plt.scatter(x1, x2) 
#plt.legend(loc="upper left", title="Scatterplot")
plt.ylabel('x2')
plt.xlabel('x1')
plt.savefig( outputPlotFileName )
#plt.show()
plt.clf()

# density = False would plot counts
plt.hist(x1, density=True, bins=64, label="Data") 
mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = st.gaussian_kde(x1)
plt.plot(kde_xs, kde.pdf(kde_xs), label="x1 Probability Density Function" )
plt.legend(loc="upper left")
plt.ylabel('Probability')
plt.xlabel('Data')
plt.title('Histogram')
plt.savefig( 'hist1' + outputPlotFileName )
plt.clf()

# density = False would plot counts
plt.hist(x2, density=True, bins=64, label="Data") 
mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = st.gaussian_kde(x2)
plt.plot(kde_xs, kde.pdf(kde_xs), label="x2 Probability Density Function" )
plt.legend(loc="upper left")
plt.ylabel('Probability')
plt.xlabel('Data')
plt.title('Histogram')
plt.savefig( 'hist2' + outputPlotFileName )
plt.clf()
plt.close('all')

quit(0)

