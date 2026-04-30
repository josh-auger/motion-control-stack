
# REQUIREMENT: shared-base-build image is a pre-requisite for the motion-control-stack container build
#docker build --no-cache -t shared-base-build:latest -f ./docker/Dockerfile.base .

#docker build --no-cache -t jauger/motion-control-stack:dev -f ./docker/Dockerfile.unified .
docker build --rm -t jauger/motion-control-stack:dev -f ./docker/Dockerfile.unified .