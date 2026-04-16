#!/bin/bash
set -e

MODE="$1"
shift || true

echo "[motion-control-stack] mode = $MODE"

# ============================================================
# FIRE SERVER
# ============================================================
if [ "$MODE" = "fire-server" ]; then

  MOCO_FLAG="${MOCO_FLAG:-off}"
  REG_TYPE="${REG_TYPE:-smsgroup}"

  echo "Starting fire-server"
  echo "MOCO_FLAG=$MOCO_FLAG"
  echo "REG_TYPE=$REG_TYPE"

  exec python3 /opt/apps/fire_server/main.py \
    -v \
    -H=0.0.0.0 \
    -p=9002 \
    -S /saved_data \
    --moco="$MOCO_FLAG" \
    --regtype="$REG_TYPE"

# ============================================================
# QUEUE PROCESSOR
# ============================================================
elif [ "$MODE" = "queue-processor" ]; then

  MOCO_FLAG="${MOCO_FLAG:-off}"

  echo "Starting queue-processor"
  echo "MOCO_FLAG=$MOCO_FLAG"

  exec python3 /opt/apps/queue_processor/process_queue_directory.py \
    /data \
    --moco "$MOCO_FLAG"

# ============================================================
# MOTION MONITOR
# ============================================================
elif [ "$MODE" = "motion-monitor" ]; then

  HEAD_RADIUS="${HEAD_RADIUS:-50}"
  MOTION_THRESH="${MOTION_THRESH:-0.3}"

  echo "Starting motion-monitor"
  echo "HEAD_RADIUS=$HEAD_RADIUS"
  echo "MOTION_THRESH=$MOTION_THRESH"

  exec python3 /opt/apps/motion_monitor/monitor_directory.py \
    -p=8080 \
    /working \
    --head_radius "$HEAD_RADIUS" \
    --motion_threshold "$MOTION_THRESH"

# ============================================================
# HELP
# ============================================================
else
  echo ""
  echo "Usage:"
  echo "  fire-server"
  echo "  motion-monitor"
  echo "  queue-processor"
  echo ""
  exit 1
fi