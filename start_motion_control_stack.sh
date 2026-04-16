
#!/bin/bash
set -e

DATA_DIR="$(pwd)/data"

MOCO_FLAG="${MOCO_FLAG:-off}"
REG_TYPE="${REG_TYPE:-smsgroup}"
HEAD_RADIUS="${HEAD_RADIUS:-50}"
MOTION_THRESH="${MOTION_THRESH:-0.3}"

echo "Starting motion-control-stack cluster..."

# Fire Server
echo "Launching fire-server..."
docker run -d --rm \
  -u $(id -u):$(id -g) \
  -p 9002:9002 \
  -v "$DATA_DIR:/data" \
  -e MOCO_FLAG="$MOCO_FLAG" \
  -e REG_TYPE="$REG_TYPE" \
  --name fire-server \
  jauger/motion-control-stack:dev fire-server
echo "Fire-server exit code: $?"

# Queue Processor
echo "Launching queue-processor..."
docker run -d --rm \
  -u $(id -u):$(id -g) \
  -v "$DATA_DIR:/data" \
  -e MOCO_FLAG="$MOCO_FLAG" \
  --name queue-processor \
  jauger/motion-control-stack:dev queue-processor
echo "Queue-processor exit code: $?"

# Motion Monitor
echo "Launching motion-monitor..."
docker run -d --rm \
  -u $(id -u):$(id -g) \
  -p 8080:8080 \
  -v "$DATA_DIR:/data" \
  -e HEAD_RADIUS="$HEAD_RADIUS" \
  -e MOTION_THRESH="$MOTION_THRESH" \
  --name motion-monitor \
  jauger/motion-control-stack:dev motion-monitor
echo "Motion-monitor exit code: $?"

echo "All services started."
echo "Use: docker logs -f fire-server"