#!/usr/bin/env python3

import SimpleITK as sitk
import numpy as np
import argparse

#
# This program reads in an affine transform representing a slice displacement
# and computes the slice displacement.

# We believe that angles are always in radians in SimpleITK.

# Parse the arguments.
parser = argparse.ArgumentParser(description='compute slice displacement.')
parser.add_argument("--transformfile", required=True)
parser.add_argument("--outputfile", required=True)
parser.add_argument("--radius", required=False, type=float, default=50.0)
args = parser.parse_args()

radius = args.radius

read_transform = sitk.ReadTransform( args.transformfile )

# First three parameters are rotation angles in radians.
# Second three parameters are translations.

euler3d = sitk.Euler3DTransform()
euler3d.SetCenter(read_transform.GetCenter())
euler3d.SetTranslation(read_transform.GetTranslation())
euler3d.SetMatrix(read_transform.GetMatrix())

# Compute the displacement:
parms = np.asarray( euler3d.GetParameters() )
displacement = abs(parms[0]*radius) + abs(parms[1]*radius) + \
    abs(parms[2]*radius) + abs(parms[3]) + abs(parms[4]) + abs(parms[5])

with open( args.outputfile, 'a') as f:
    f.write( str(displacement) + '\n' )

exit(1)

