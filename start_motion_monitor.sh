
# User settings
# ================================
# Head radius assumption (mm) for displacement calculations
HEAD_RADIUS=50
# Threshold for framewise displacement (mm) to flag motion
MOTION_THRESH=0.3

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
  jauger/motion-control-stack:dev motion-monitor