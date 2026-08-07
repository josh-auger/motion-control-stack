
# REQUIREMENT: shared-base-build image is a pre-requisite for the motion-control-stack container build
# docker build --no-cache -t shared-base-build:latest -f ./docker/Dockerfile.base .

# docker build --no-cache -t jauger/motion-control-stack:dev -f ./docker/Dockerfile.unified .

# Build the motion-control-stack container with CUDA support for SLIMM-v3
docker build \
    --build-context slimm=../SLIMM-v3 \
    -t jauger/motion-control-stack:cuda \
    -f ./docker/Dockerfile.unified \
    .
