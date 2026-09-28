#!/usr/bin/env bash

# REQUIREMENT: shared-base-build image is a pre-requisite for the motion-control-stack container build
# docker build --no-cache -t shared-base-build:latest -f ./docker/Dockerfile.base .

# docker build --no-cache -t jauger/motion-control-stack:dev -f ./docker/Dockerfile.unified .

# Build the motion-control-stack container with CUDA support for SLIMM-v3
IMAGE_TAG=jauger/motion-control-stack:cuda

if [ "$#" -gt 1 ]; then
    echo "Usage: $0 [--mars]" >&2
    exit 2
fi

case "${1:-}" in
    --mars)
        IMAGE_TAG=jauger/motion-control-stack:mars-cuda-11.6-sm86
        set -- --target runtime-mars
        ;;
    "")
        set --
        ;;
    *)
        echo "Usage: $0 [--mars]" >&2
        exit 2
        ;;
esac

docker build \
    "$@" \
    --build-context slimm=../SLIMM-v3 \
    --build-context sms_mi_reg=../sms-mi-reg \
    -t "$IMAGE_TAG" \
    -f ./docker/Dockerfile.unified \
    .
