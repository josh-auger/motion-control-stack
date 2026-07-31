# docker build --no-cache -t shared-base-build:latest -f ./docker/Dockerfile.base .

docker build --target build-base   -t shared-base-build:devel   -f dockerfile.base .
docker build --target runtime-base -t shared-base-build:runtime -f dockerfile.base .
