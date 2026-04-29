#!/bin/bash

# DOCKER-COMPOSE RUN
# ================================
# Bash script to start stack of motion monitoring and control containers ("motion_control_stack"):
#   fire-server (detached) = python fire server for image data send/receive with scanner
#   queue-processor (detached) = execute image registration of compiled image data
#   motion-monitor (with foreground terminal log) = characterize and live stream motion results
#
# This requires first making the bash script executable:
#   chmod +x start_compose_motion_control_stack.sh
# Then run with:
#   ./start_compose_motion_control_stack.sh
#
# To gracefully stop the composed services:
#   ctrl+c (possibly twice) in the same shell instance

set -e
export HOST_UID=$(id -u)
export HOST_GID=$(id -g)
export DATA_DIR="${DATA_DIR:-$(pwd)/data}"

if [ ! -d "$DATA_DIR" ]; then
  echo "ERROR: DATA_DIR does not exist: $DATA_DIR"
  exit 1
fi

echo "Starting motion control stack..."
echo "  DATA_DIR: $DATA_DIR"
echo

cleanup() {
  echo
  echo "Stopping motion control stack..."
  docker compose down
}

#trap cleanup EXIT
trap cleanup SIGINT SIGTERM
docker compose up -d
#docker compose ps
docker compose logs -f motion-monitor