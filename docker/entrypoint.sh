#!/bin/bash
set -e

MODE="$1"
shift || true

# Set working directory for all generated outputs
# Docker launcher should override default to use WORKDIR=/data
# MARS chroot launch will default to use WORKDIR=/tmp/share
WORKDIR="${WORKDIR:-/tmp/share}"

# Verify working directory before attempting to launch services
if [ ! -d "$WORKDIR" ]; then
  echo "ERROR: Working directory does not exist: $WORKDIR"
  exit 1
fi

# Set up entrypoint logging
LOGDIR="$WORKDIR/logs_motioncontrolstack"
mkdir -p "$LOGDIR"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOGFILE="$LOGDIR/log_entrypoint_${TIMESTAMP}.log"
log() {
    echo "$@"
    echo "$@" >> "$LOGFILE"
}

log "========================================"
log "Motion Control Stack"
log "Started: $(date)"
log "MODE=$MODE"
log "WORKDIR=$WORKDIR"
log "========================================"


# ============================================================
# SIMULTANEOUS MULTI-SERVICE MODE
# ============================================================
if [ "$MODE" = "all" ]; then
  MOCO_FLAG="${MOCO_FLAG:-off}"
  REG_TYPE="${REG_TYPE:-smsgroup}"
  FIFO_FLAG="${FIFO_FLAG:-on}"
  HEAD_RADIUS="${HEAD_RADIUS:-50}"
  MOTION_THRESH="${MOTION_THRESH:-0.3}"
  STREAM_FLAG="${STREAM_FLAG:-off}"
  
  log "Starting ALL services..."
  log "  WORKDIR=$WORKDIR"
  log "  MOCO_FLAG=$MOCO_FLAG"
  log "  REG_TYPE=$REG_TYPE"
  log "  FIFO_FLAG=$FIFO_FLAG"
  log "  HEAD_RADIUS=$HEAD_RADIUS"
  log "  MOTION_THRESH=$MOTION_THRESH"
  log "  STREAM_FLAG=$STREAM_FLAG"

  # ---- Start services in background ----
  log "Launching fire-server..."
  python3 /opt/apps/fire_server/main.py \
    -v \
    -H=0.0.0.0 \
    -p=9002 \
    -S "$WORKDIR" \
    --moco="$MOCO_FLAG" \
    --regtype="$REG_TYPE" \
    > "$LOGDIR/STDOUT_fire_server_${TIMESTAMP}.log" 2>&1 &   # redirect stdout/stderr to log file
  FIRE_PID=$!

  sleep 0.5  # give fire-server a moment to start, then check if it exited immediately
  if kill -0 "$FIRE_PID" 2>/dev/null; then
      log "  fire-server started successfully (PID=$FIRE_PID)"
  else
      log "  ERROR: fire-server exited immediately"
  fi


  log "Launching queue-processor..."
  python3 /opt/apps/queue_processor/process_queue_directory.py \
    "$WORKDIR" \
    --fifo "$FIFO_FLAG" \
    > "$LOGDIR/STDOUT_queue_processor_${TIMESTAMP}.log" 2>&1 &   # redirect stdout/stderr to log file
  QUEUE_PID=$!

  sleep 0.5  # give queue-processor a moment to start, then check if it exited immediately
  if kill -0 "$QUEUE_PID" 2>/dev/null; then
      log "  queue-processor started successfully (PID=$QUEUE_PID)"
  else
      log "  ERROR: queue-processor exited immediately"
  fi


  log "Launching motion-monitor..."
  python3 /opt/apps/motion_monitor/monitor_directory.py \
    -p=8080 \
    "$WORKDIR" \
    --webstream "$STREAM_FLAG" \
    --head_radius "$HEAD_RADIUS" \
    --motion_threshold "$MOTION_THRESH" \
    > "$LOGDIR/STDOUT_motion_monitor_${TIMESTAMP}.log" 2>&1 &    # redirect stdout/stderr to log file
  MONITOR_PID=$!

  sleep 0.5  # give motion-monitor a moment to start, then check if it exited immediately
  if kill -0 "$MONITOR_PID" 2>/dev/null; then
      log "  motion-monitor started successfully (PID=$MONITOR_PID)"
  else
      log "  ERROR: motion-monitor exited immediately"
  fi


  log "PIDs:"
  log "  fire-server: $FIRE_PID"
  log "  queue-processor: $QUEUE_PID"
  log "  motion-monitor: $MONITOR_PID"

  # ---- Signal handling ----
  shutdown() {
    STATUS=$?
    log ""
    log "========================================"
    log "Shutdown initiated"
    log "Exit status: $STATUS"
    log "========================================"
    kill -TERM $FIRE_PID $QUEUE_PID $MONITOR_PID 2>/dev/null || true
    wait
    exit 0
  }

  trap shutdown SIGINT SIGTERM

  # ---- Wait for any process to exit ----
  wait -n
  log "One service exited. Shutting down others..."
  shutdown



# ============================================================
# SINGLE SERVICE MODES
# ============================================================
# FIRE SERVER
elif [ "$MODE" = "fire-server" ]; then
  MOCO_FLAG="${MOCO_FLAG:-off}"
  REG_TYPE="${REG_TYPE:-smsgroup}"
  log "Starting fire-server"
  log "  WORKDIR=$WORKDIR"
  log "  MOCO_FLAG=$MOCO_FLAG"
  log "  REG_TYPE=$REG_TYPE"

  exec python3 /opt/apps/fire_server/main.py \
    -v \
    -H=0.0.0.0 \
    -p=9002 \
    -S "$WORKDIR" \
    --moco="$MOCO_FLAG" \
    --regtype="$REG_TYPE"


# QUEUE PROCESSOR
elif [ "$MODE" = "queue-processor" ]; then
  FIFO_FLAG="${FIFO_FLAG:-on}"
  log "Starting queue-processor"
  log "  WORKDIR=$WORKDIR"
  log "  FIFO_FLAG=$FIFO_FLAG"

  exec python3 /opt/apps/queue_processor/process_queue_directory.py \
    "$WORKDIR" \
    --fifo "$FIFO_FLAG"


# MOTION MONITOR
elif [ "$MODE" = "motion-monitor" ]; then
  HEAD_RADIUS="${HEAD_RADIUS:-50}"
  MOTION_THRESH="${MOTION_THRESH:-0.3}"
  STREAM_FLAG="${STREAM_FLAG:-off}"
  log "Starting motion-monitor"
  log "  WORKDIR=$WORKDIR"
  log "  HEAD_RADIUS=$HEAD_RADIUS"
  log "  MOTION_THRESH=$MOTION_THRESH"
  log "  STREAM_FLAG=$STREAM_FLAG"

  exec python3 /opt/apps/motion_monitor/monitor_directory.py \
    -p=8080 \
    "$WORKDIR" \
    --webstream "$STREAM_FLAG" \
    --head_radius "$HEAD_RADIUS" \
    --motion_threshold "$MOTION_THRESH"


# HELP
else
  log "Usage:"
  log "  fire-server"
  log "  motion-monitor"
  log "  queue-processor"
  log ""
  exit 1
fi