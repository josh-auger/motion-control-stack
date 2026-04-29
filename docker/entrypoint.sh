#!/bin/bash
set -e

MODE="$1"
shift || true

echo "[motion-control-stack] mode = $MODE"

# ============================================================
# SIMULTANEOUS MULTI-SERVICE MODE
# ============================================================
if [ "$MODE" = "all" ]; then
  MOCO_FLAG="${MOCO_FLAG:-off}"
  REG_TYPE="${REG_TYPE:-smsgroup}"
  FIFO_FLAG="${FIFO_FLAG:-on}"
  HEAD_RADIUS="${HEAD_RADIUS:-50}"
  MOTION_THRESH="${MOTION_THRESH:-0.3}"
  echo "Starting ALL services..."

  # ---- Start services in background ----
  echo "MOCO_FLAG=$MOCO_FLAG"
  echo "REG_TYPE=$REG_TYPE"
  python3 /opt/apps/fire_server/main.py \
    -v -H=0.0.0.0 -p=9002 -S /data \
    --moco="$MOCO_FLAG" \
    --regtype="$REG_TYPE" \
    > /dev/null 2>&1 &    # suppress logs in terminal
  FIRE_PID=$!

  echo "FIFO_FLAG=$FIFO_FLAG"
  python3 /opt/apps/queue_processor/process_queue_directory.py \
    /data \
    --fifo "$FIFO_FLAG" \
    > /dev/null 2>&1 &    # suppress logs in terminal
  QUEUE_PID=$!

  echo "HEAD_RADIUS=$HEAD_RADIUS"
  echo "MOTION_THRESH=$MOTION_THRESH"
  python3 /opt/apps/motion_monitor/monitor_directory.py \
    -p=8080 \
    /data \
    --head_radius "$HEAD_RADIUS" \
    --motion_threshold "$MOTION_THRESH" &
  MONITOR_PID=$!

  echo "PIDs:"
  echo "  fire-server: $FIRE_PID"
  echo "  queue-processor: $QUEUE_PID"
  echo "  motion-monitor: $MONITOR_PID"

  # ---- Signal handling ----
  shutdown() {
    echo "Shutting down services..."
    kill -TERM $FIRE_PID $QUEUE_PID $MONITOR_PID 2>/dev/null || true
    wait
    exit 0
  }

  trap shutdown SIGINT SIGTERM

  # ---- Wait for any process to exit ----
  wait -n
  echo "One service exited. Shutting down others..."
  shutdown


# ============================================================
# SINGLE SERVICE MODES
# ============================================================
# FIRE SERVER
elif [ "$MODE" = "fire-server" ]; then
  MOCO_FLAG="${MOCO_FLAG:-off}"
  REG_TYPE="${REG_TYPE:-smsgroup}"
  echo "Starting fire-server"
  echo "MOCO_FLAG=$MOCO_FLAG"
  echo "REG_TYPE=$REG_TYPE"

  exec python3 /opt/apps/fire_server/main.py \
    -v \
    -H=0.0.0.0 \
    -p=9002 \
    -S /data \
    --moco="$MOCO_FLAG" \
    --regtype="$REG_TYPE"


# QUEUE PROCESSOR
elif [ "$MODE" = "queue-processor" ]; then
  FIFO_FLAG="${FIFO_FLAG:-on}"
  echo "Starting queue-processor"
  echo "FIFO_FLAG=$FIFO_FLAG"

  exec python3 /opt/apps/queue_processor/process_queue_directory.py \
    /data \
    --fifo "$FIFO_FLAG"


# MOTION MONITOR
elif [ "$MODE" = "motion-monitor" ]; then
  HEAD_RADIUS="${HEAD_RADIUS:-50}"
  MOTION_THRESH="${MOTION_THRESH:-0.3}"
  echo "Starting motion-monitor"
  echo "HEAD_RADIUS=$HEAD_RADIUS"
  echo "MOTION_THRESH=$MOTION_THRESH"

  exec python3 /opt/apps/motion_monitor/monitor_directory.py \
    -p=8080 \
    /data \
    --head_radius "$HEAD_RADIUS" \
    --motion_threshold "$MOTION_THRESH"


# HELP
else
  echo ""
  echo "Usage:"
  echo "  fire-server"
  echo "  motion-monitor"
  echo "  queue-processor"
  echo ""
  exit 1
fi