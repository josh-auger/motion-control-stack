
# User settings
# ================================
# Toggles processing every file sequentially ("first-in-first-out" = "on") or only the most recent file (FIFO = "off")
FIFO_FLAG="on"

# Local host directory for output files
DATA_DIR=$(pwd)"/data"


# Docker run command
# ================================
docker run --rm -it \
  -u $(id -u):$(id -g) \
  -v $DATA_DIR:/data \
  -e FIFO_FLAG="$FIFO_FLAG" \
  jauger/motion-control-stack:dev queue-processor