
# User settings
# ================================
# Moco flag ("on", "off") toggles processing every file sequentially (moco "off") or only the most recent file (moco "on")
MOCO_FLAG="on"

# Local host directory for output files
DATA_DIR=$(pwd)"/data"


# Docker run command
# ================================
docker run --rm -it \
  -u $(id -u):$(id -g) \
  -p 9002:9002 \
  -v $DATA_DIR:/data \
  -e MOCO_FLAG="$MOCO_FLAG" \
  jauger/motion-control-stack:dev queue-processor