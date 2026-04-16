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

reader = sitk.ImageFileReader()
#reader.SetImageIO("NRRDImageIO")
reader.SetFileName( "initialHistogram.nrrd" )
image = reader.Execute();
 
ssfilter = sitk.ShiftScaleImageFilter()
ssfilter.SetShift(1.0)
ssfilter.SetScale(1.0)
image = ssfilter.Execute( image )

logfilter = sitk.LogImageFilter()
image = logfilter.Execute(image) 

writer = sitk.ImageFileWriter()
writer.SetFileName( "scaled-initialHistogram.nrrd" )
writer.Execute(image)

nda = sitk.GetArrayFromImage( image )
#plt.imshow( nda , cmap = 'gray', interpolation='bilinear' )

#
# FIX: Now read this in from the marginal histograms
# and from the input data.
#

#np.random.seed( 42 )
#x = np.random.normal( size = 1000 )
x = nda.flatten()
# density = False would plot counts
plt.hist(x, density=True, bins=64, label="Data") 
mn, mx = plt.xlim()
plt.xlim(mn, mx)
kde_xs = np.linspace(mn, mx, 300)
kde = st.gaussian_kde(x)
plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")
plt.legend(loc="upper left")
plt.ylabel('Probability')
plt.xlabel('Data')
plt.title('Histogram')
plt.show()

