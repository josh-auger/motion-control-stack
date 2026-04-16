
# User settings
# ================================
# Moco flag ("on", "off") toggles the sending of moco feedback to the scanner for prospective motion correction
MOCO_FLAG="off"

# Registration type (by "slice", "smsgroup", or "volume") determines the grouping of image data
REG_TYPE="smsgroup"

# Local host directory for output files
DATA_DIR=$(pwd)"/data"


# Docker run command
# ================================
docker run --rm -it \
  -u $(id -u):$(id -g) \
  -p 9002:9002 \
  -v $DATA_DIR:/saved_data \
  -e MOCO_FLAG="$MOCO_FLAG" \
  -e REG_TYPE="$REG_TYPE" \
  jauger/motion-control-stack:dev fire-server