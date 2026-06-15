# Motion-Control-Stack

The **motion-control-stack** implements a real-time, multi-process pipeline for image-based motion estimation and prospective motion correction (MOCO).  

The system operates as a **closed-loop control pipeline**, continuously ingesting image data, estimating motion, and optionally feeding corrections back to the scanner.

---
## 🧱 Core Components

### 🔥 Fire-Server
- Listens for incoming MRD image data stream (port `9002`)
- Compiles streamed image data into:
  - Series metadata (`.json`)
  - NRRD format detached image header (`.nhdr`)
  - Raw byte string of image data (`.raw`)
  - Pointer file denoting acquisition groups (`.txt`)

#### When motion correction (MoCo) is enabled:
- Monitors data directory for new transform files (`.tfm`)
- Converts alignment transforms into scanner's global coordinate frame
- Sends MoCo feedback message to scanner to update imaging field-of-view

---
### 🔁 Queue-Processor
- Monitors data directory for new pointer files (`.txt`) from fire-server
- Loads image data (`.nhdr` + `.raw`) referenced in the pointer file
- Executes image registration using `sms-mi-reg`
- Outputs alignment transforms (`.tfm`) that map the reference volume to the new target slice(s)

#### When first-in-first-out (FIFO) is enabled:
- First-in-first-out (FIFO)
- Processes all incoming data sequentially as acquired

#### When FIFO is disabled:
- Last-in-first-out
- Drops stale pointer files to process only the most recent image data received
- Minimizes latency for real-time correction to the latest known position

---
### 📈 Motion-Monitor
- Monitors data directory for new alignment transform files (`.tfm`) from queue-processor registration
- Maintains a running ledger of alignment transforms
- Computes **framewise displacement** between sequential transforms
- Classifies motion events based on set motion threshold
- Generates a motion summary figure and livestreams to a web interface (`http://<host IP address>:8080/stream.mjpg`)

---
## 📂 Data Flow
All components communicate through a shared filesystem (typically mounted as `/data`).
Future developments will focus on handling all data in working memory.

---
## ⚙️ Runtime Configuration
The motion-control stack is configured through a single bash run script (`start_motion_control_stack.sh`). 
User-defined configuration parameters are specified at the top of the script and passed into the Docker container as 
environment variables at runtime.

### Configuration Parameters
#### General
| Parameter | Description                                                                                                                    |
|-----------|--------------------------------------------------------------------------------------------------------------------------------|
| `DATA_DIR` | Host directory mounted into container at `/data`. Used for logs, temp files, registration outputs, and motion-monitor results. |

#### Fire-Server
| Parameter | Options | Description |
|-----------|---------|-------------|
| `MOCO_FLAG` | `on`, `off` | Enables or disables prospective motion-correction feedback packets sent to the scanner. |
| `REG_TYPE` | `slice`, `smsgroup`, `volume` | Defines the image grouping used as the registration unit. |

#### Queue-Processor
| Parameter | Options | Description |
|-----------|---------|-------------|
| `FIFO_FLAG` | `on`, `off` | Controls queue behavior. `on` processes all incoming data sequentially. `off` skips older items and processes only the most recent data available. |

#### Motion-Monitor
| Parameter | Units | Description |
|-----------|-------|-------------|
| `HEAD_RADIUS` | mm | Assumed head radius used when converting rotational motion into displacement metrics. |
| `MOTION_THRESH` | mm | Framewise displacement threshold used to flag excessive motion. |

The example run script `start_motion_control_stack.sh` is thus written as:
```bash
DATA_DIR=$(pwd)"/data"
MOCO_FLAG="off"
REG_TYPE="smsgroup"
FIFO_FLAG="off"
HEAD_RADIUS=50
MOTION_THRESH=0.3

docker run --rm -it \
  -u $(id -u):$(id -g) \
  -p 9002:9002 \
  -p 8080:8080 \
  -v $DATA_DIR:/data \
  -e MOCO_FLAG="$MOCO_FLAG" \
  -e FIFO_FLAG="$FIFO_FLAG" \
  -e REG_TYPE="$REG_TYPE" \
  -e HEAD_RADIUS="$HEAD_RADIUS" \
  -e MOTION_THRESH="$MOTION_THRESH" \
  jauger/motion-control-stack:dev all
```

### Starting the Stack
After configuring the desired parameters in the bash run script, start the complete motion-control stack with:
```bash
sh start_motion_control_stack.sh
```

This launches all motion-control sub-process services within a single Docker container:
- **Fire-Server** receives and collates image data from the scanner and optionally sends prospective motion-correction feedback back to the scanner.
- **Queue-Processor** performs image registration and motion estimation on the received image data.
- **Motion-Monitor** computes and streams motion metrics in real time.

During operation, the services automatically detect the end of an imaging sequence (including early termination by the MRI control computer), reset their internal state, and wait for the next sequence to begin.

All services continue running until the container is explicitly stopped.

### Stopping the Stack
#### Option 1: Stop from the Launch Terminal
If the container is running in the foreground terminal, stop all services and shut down the Docker container by pressing `Ctrl+C`.

This triggers the container shutdown handler and gracefully terminates all running services.

#### Option 2: Stop the Docker Container Directly
If the stack is running in the background or from a different terminal session, first identify the running container with:

```bash
docker ps
```

or

```bash
docker container ls
```

Example output:

```text
CONTAINER ID   IMAGE                            STATUS
a1b2c3d4e5f6   jauger/motion-control-stack:dev  Up 15 minutes
```

Then, stop the specific container using its container ID:

```bash
docker stop a1b2c3d4e5f6
```

Docker will send a termination signal to the container and gracefully stop all motion-control services.

#### Verify Shutdown
To verify that the container is no longer running:

```bash
docker ps -a
```

The motion-control-stack container should no longer appear in the list of running containers.

---
## ⚠️ Design Considerations

### Real-Time Processing Constraints
The stack is designed for real-time MRI motion monitoring and prospective motion correction. Depending on the selected 
configuration, a trade-off exists between:
- **Completeness** – Process every image received from the scanner. `FIFO_FLAG=on` processes all incoming image data 
sequentially in the order it was acquired by the scanner.
- **Latency** – Produce the most current motion estimate as quickly as possible. `FIFO_FLAG=off` discards stale queue 
entries and processes only the newest available data.

### Shared-Volume Communication
All services communicate through shared local files stored in the mounted Docker volume (`/data`).
Advantages:
- Simple and transparent data flow
- Easy debugging and inspection of intermediate results

Considerations:
- Multiple services may access files simultaneously
- Partial writes or race conditions must be handled carefully
- File I/O performance can directly impact end-to-end latency

---
## 🔄 Configuration for Prospective Motion Correction
For prospective motion correction (MOCO), low latency is generally preferred over processing every image.

Therefore, the recommended settings are:
```bash
MOCO_FLAG="on"
FIFO_FLAG="off"
```

Workflow:
1. Fire-Server receives image data from the scanner.
2. Queue-Processor prioritizes the most recent image data and computes the current subject position.
3. Fire-Server converts the registration transform into the scanner coordinate system.
4. Fire-Server transmits a motion-correction (MoCo) feedback packet back to the scanner.
5. The scanner receives the MoCo feedback packet and updates the imaging field-of-view for all subsequent acquisitions, 
or until a new MoCo feedback packet is received.

This forms a closed-loop control system that continuously estimates subject motion and updates scan geometry during 
image acquisition. 

⚠️ IMPORTANT NOTE: The Motion-Monitor service is agnostic to any image queue clearing and will calculate framewise displacements between whichever two subsequent alignment transforms are produced, and transmit those motion results. If, in the interest of latency, stale image data has been skipped in favor of more recent image data (`FIFO_FLAG` = `off`), then the motion traces reported by the Motion-Monitor will also have that data omitted. For a more complete trace of a subject's continuous motion pattern during a scan, it is recommended to run a retrospective motion analysis on the un-corrected image data.