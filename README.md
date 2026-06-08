# Motion-Control-Stack

The **motion-control-stack** implements a real-time, multi-process pipeline for image-based motion estimation and 
prospective motion correction (MOCO).  

The system operates as a **closed-loop control pipeline**, continuously ingesting image data, estimating motion, and 
optionally feeding corrections back to the scanner.

---
## 🔄 High-Level Pipeline
```mermaid
flowchart LR

A[MRD Stream] --> B[Fire Server]

B --> B1[Parse MRD]
B1 --> B2[Write .nhdr]
B1 --> B3[Write .raw]
B1 --> B4[Write pointer .txt]

B4 --> C[Queue Processor]
C --> C1[Load .nhdr + .raw]
C1 --> C2[Run sms-mi-reg Registration]

C2 --> D[Transform File (.tfm)]

D --> E[Motion Monitor]
E --> E1[Framewise Displacement]
E1 --> E2[Motion Classification]
E2 --> E3[Live Visualization]

D --> F{MOCO Enabled?}

F -->|Yes| G[Fire Server MOCO Handler]
G --> H[Convert to Scanner Frame]
H --> I[Send Feedback Packet]
I --> A

subgraph S[Shared Data Volume (/data)]
B2
B3
B4
D
end
```

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

### Running the Stack
After configuring the desired parameters in the bash run script, start all services with:
```bash
sh start_motion_control_stack.sh
```

This launches the complete motion-control stack of services:
- **Fire-Server** receives and collates image data from the scanner and optionally sends prospective motion-correction
feedback back to the scanner.
- **Queue-Processor** performs image registration and motion estimation on the received image data.
- **Motion-Monitor** computes and streams motion metrics in real time.

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
## ⚡ Configuration for Prospective Motion Correction
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

⚠️ It is important to note that the Motion-Monitor is agnostic to any image queue clearing and will calculate framewise 
displacements between whichever two subsequent alignment transforms are produced, and transmit those motion results. 
If, in the interest of latency, stale image data has been skipped in favor of more recent image data, the motion traces 
reported by the Motion-Monitor will also have that data omitted. For a more complete trace of a subject's continuous 
motion pattern during a scan, it is recommended to run a retrospective motion analysis on the un-corrected image data.