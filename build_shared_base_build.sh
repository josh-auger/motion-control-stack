# docker build --no-cache -t shared-base-build:latest -f ./docker/Dockerfile.base .

docker build --no-cache --target build-base   -t shared-base-build:devel   -f ./docker/Dockerfile.base .
docker build --no-cache --target runtime-base -t shared-base-build:runtime -f ./docker/Dockerfile.base .
