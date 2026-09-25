
# User settings
# ================================
# Head radius assumption (mm) for displacement calculations
HEAD_RADIUS=50
# Threshold for framewise displacement (mm) to flag motion
MOTION_THRESH=0.3
# Periodic motion-dashboard interval in monotonic elapsed seconds
MOTION_DASHBOARD_INTERVAL_SEC=5.0

# Local host directory for output files
DATA_DIR=$(pwd)"/data"


# Docker run command
# ================================
docker run --rm -it \
  -u $(id -u):$(id -g) \
  -p 8080:8080 \
  -v $DATA_DIR:/data \
  -e HEAD_RADIUS=$HEAD_RADIUS \
  -e MOTION_THRESH=$MOTION_THRESH \
  -e MOTION_DASHBOARD_INTERVAL_SEC=$MOTION_DASHBOARD_INTERVAL_SEC \
  jauger/motion-control-stack:dev motion-monitor
