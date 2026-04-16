# Experiment with manipulating the image metadata to position it in space.
## Align an image to a rotated version of itself.
- Select an input image.
  - Display the image geometry information
    - crl-image-info.py
- Modify the input image to have axes in the LPS+ frame.
  - This may modify the origin and the order of voxels in the file.
  - The geometric location of each voxel will not change.
  - crl-orient-image.py
- Modify the input image to set the origin to be (0.0, 0.0, 0.0)
  - crl-translate-origin-to-zero.py
- Make an affine transform from the space directions
  - crl-space-directions-to-transform.py 

- Modify the image to have identity space directions
  - crl-identity-space-directions.py

crl-reorient-image.py --input t1w.nhdr --output t1w-reorient.nhdr
crl-translate-origin-to-zero.py --input t1w-reorient.nhdr \
  --output t1w-reorient-ozero.nhdr
crl-identity-space-directions.py --input t1w-reorient-ozero.nhdr --output t1w-reorient-ozero-directions-identity.nhdr

crl-space-directions-to-transform.py --input t1w-reorient-ozero.nhdr --transformfile directions.tfm

docker run --volume `pwd`:/data --rm -it --init crl/sms-mi-reg sms-mi-reg  t1w-reorient-ozero.nhdr  identity.tfm  align  t1w-reorient-ozero-directions-identity.nhdr


