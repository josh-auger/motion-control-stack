
# REQUIREMENT: shared-base-build image is a pre-requisite for the motion-control-stack container build
# docker build --no-cache -t shared-base-build:latest -f ./docker/Dockerfile.base .

# docker build --no-cache -t jauger/motion-control-stack:dev -f ./docker/Dockerfile.unified .

# Temporarily stage SLIMM-v3 for Docker build
rm -rf external/SLIMM-v3
cp -r ../SLIMM-v3 external/SLIMM-v3

# Build the motion-control-stack container with CUDA support for SLIMM-v3
docker build -t jauger/motion-control-stack:cuda -f ./docker/Dockerfile.unified .

# Remove staged copy
rm -rf external/SLIMM-v3
