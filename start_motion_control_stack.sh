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
HOST_DATA_DIR="$(pwd)/data"

# Working directory inside container for output files (Docker uses /data, MARS chroot uses /tmp/share
WORKDIR="/data"

# Toggle sending moco feedback to the scanner for prospective motion correction ("on", "off")
MOCO_FLAG="on"

# Registration engine determines the backend used for image registration ("sms-mi-reg", "cuda")
REG_ENGINE="cuda"

# CUDA process lifetime ("standalone", "persistent")
CUDA_EXECUTION_MODE="persistent"

# Registration type determines the grouping of image data (by "slice", "smsgroup", "volume", "by2slices") 
REG_TYPE="by2slices"

# Toggle processing data sequentially ("first-in-first-out" = "on") or only the most recent data (FIFO = "off" = "last-in-first-out")
FIFO_FLAG="off"

# Queue-processor lifecycle profiling ("on", "off"); leave off for normal runs
QUEUE_PROFILE_FLAG="on"

# Head radius assumption (mm) for displacement calculations
HEAD_RADIUS=50

# Threshold for framewise displacement (mm) to flag motion
MOTION_THRESH=0.3

# Toggle web streaming in motion-monitor ("on", "off")
STREAM_FLAG="off"

# Periodic motion-dashboard interval in monotonic elapsed seconds
MOTION_DASHBOARD_INTERVAL_SEC=5.0

# Motion-dashboard generation profiling ("on", "off"); leave off for normal runs
MOTION_DASHBOARD_PROFILE_FLAG="on"

# tSNR slice assembly ("physical", "reference_grid"). The conservative default
# preserves saved physical geometry; reference_grid measures scanner-visible pixels.
TSNR_ASSEMBLY_MODE="reference_grid"

# Toggle sending motion report back to scanner in-line display ("on", "off")
SEND_DASHBOARD_FLAG="off"


# Docker Run Command
# =====================================
# set -x
docker run --rm -it \
  -u $(id -u):$(id -g) \
  -p 9002:9002 \
  -p 8080:8080 \
  -v "$HOST_DATA_DIR:$WORKDIR" \
  -e WORKDIR="$WORKDIR" \
  -e MOCO_FLAG="$MOCO_FLAG" \
  -e FIFO_FLAG="$FIFO_FLAG" \
  -e QUEUE_PROFILE_FLAG="$QUEUE_PROFILE_FLAG" \
  -e REG_ENGINE="$REG_ENGINE" \
  -e CUDA_EXECUTION_MODE="$CUDA_EXECUTION_MODE" \
  -e REG_TYPE="$REG_TYPE" \
  -e HEAD_RADIUS="$HEAD_RADIUS" \
  -e MOTION_THRESH="$MOTION_THRESH" \
  -e STREAM_FLAG="$STREAM_FLAG" \
  -e MOTION_DASHBOARD_INTERVAL_SEC="$MOTION_DASHBOARD_INTERVAL_SEC" \
  -e MOTION_DASHBOARD_PROFILE_FLAG="$MOTION_DASHBOARD_PROFILE_FLAG" \
  -e TSNR_ASSEMBLY_MODE="$TSNR_ASSEMBLY_MODE" \
  -e SEND_DASHBOARD_FLAG="$SEND_DASHBOARD_FLAG" \
  --gpus all \
  jauger/motion-control-stack:cuda all
