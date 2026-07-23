# SIMULTANEOUS MULTI-SERVICE RUN
# =====================================
# Bash script to start stack of motion monitoring and control services ("motion_control_stack"):
#   fire-server = python fire server for image data send/receive with scanner
#   queue-processor = execute image registration of compiled image data
#   motion-monitor = characterize and live stream motion results
#
# Example run command:
#   sh start_motion_control_stack.sh


# User-defined Configuration Parameters
# =====================================
# Local host directory for output files
HOST_DATA_DIR=$(pwd)"/data"

# Working directory inside container for output files (Docker uses /data, MARS chroot uses /tmp/share
WORKDIR="/data"

# Toggle sending moco feedback to the scanner for prospective motion correction ("on", "off")
MOCO_FLAG="off"

# Registration type determines the grouping of image data (by "slice", "smsgroup", or "volume") 
REG_TYPE="smsgroup"

# Toggle processing data sequentially ("first-in-first-out" = "on") or only the most recent data (FIFO = "off" = "last-in-first-out")
FIFO_FLAG="on"

# Head radius assumption (mm) for displacement calculations
HEAD_RADIUS=50

# Threshold for framewise displacement (mm) to flag motion
MOTION_THRESH=0.3

# Toggle web streaming in motion-monitor ("on", "off")
STREAM_FLAG="on"


# Docker Run Command
# =====================================
docker run --rm -it \
  -u $(id -u):$(id -g) \
  -p 9002:9002 \
  -p 8080:8080 \
  -v $HOST_DATA_DIR:"$WORKDIR" \
  -e WORKDIR="$WORKDIR" \
  -e MOCO_FLAG="$MOCO_FLAG" \
  -e FIFO_FLAG="$FIFO_FLAG" \
  -e REG_TYPE="$REG_TYPE" \
  -e HEAD_RADIUS="$HEAD_RADIUS" \
  -e MOTION_THRESH="$MOTION_THRESH" \
  -e STREAM_FLAG="$STREAM_FLAG" \
  jauger/motion-control-stack:dev all