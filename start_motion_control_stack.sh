# SIMULTANEOUS MULTI-SERVICE RUN
# ================================
# Bash script to start stack of motion monitoring and control services ("motion_control_stack"):
#   fire-server = python fire server for image data send/receive with scanner
#   queue-processor = execute image registration of compiled image data
#   motion-monitor = characterize and live stream motion results
#
# Example run command:
#   sh start_motion_control_stack.sh


# User settings
# ================================
# Local host directory for output files
DATA_DIR=$(pwd)"/data"
# Moco flag ("on", "off") toggles sending moco feedback to the scanner for prospective motion correction
MOCO_FLAG="off"
# Registration type (by "slice", "smsgroup", or "volume") determines the grouping of image data
REG_TYPE="smsgroup"
# Toggles processing data sequentially ("first-in-first-out" = "on") or only the most recent data (FIFO = "off")
FIFO_FLAG="on"
# Head radius assumption (mm) for displacement calculations
HEAD_RADIUS=50
# Threshold for framewise displacement (mm) to flag motion
MOTION_THRESH=0.3

docker run --rm -it \
  -u $(id -u):$(id -g) \
  -p 9002:9002 \
  -p 8080:8080 \
  -v $(pwd)/data:/data \
  -e MOCO_FLAG=off \
  -e FIFO_FLAG=on \
  -e REG_TYPE=smsgroup \
  -e HEAD_RADIUS=50 \
  -e MOTION_THRESH=0.3 \
  jauger/motion-control-stack:dev all