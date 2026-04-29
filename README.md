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

### 🔥 Fire Server
- Listens for incoming MRD image data stream (port `9002`)
- Compiles streamed image data into:
  - Series metadata (`.json`)
  - NRRD format detached image header (`.nhdr`)
  - Raw byte string of image data (`.raw`)
  - Pointer file denoting acquisition groups (`.txt`)

#### When MOCO is enabled:
- Monitors data directory for new transform files (`.tfm`)
- Converts alignment transforms into scanner's global coordinate frame
- Sends MoCo feedback message to scanner to update imaging field-of-view

---
### 🔁 Queue Processor
- Monitors data directory for new pointer files (`.txt`) from fire-server
- Loads image data (`.nhdr` + `.raw`) referenced in the pointer file
- Executes image registration using `sms-mi-reg`
- Outputs alignment transforms (`.tfm`) that map the reference volume to the new target slice(s)

#### Queue Behavior Modes
**FIFO ON**
- First-in-first-out (FIFO)
- Processes all incoming data sequentially as acquired

**FIFO OFF**
- Last-in-first-out
- Drops stale pointer files to process only the most recent image data received
- Minimizes latency for real-time correction to the latest known position

---
### 📈 Motion Monitor
- Monitors data directory for new alignment transform files (`.tfm`) from queue-processor registration
- Maintains a running ledger of alignment transforms
- Computes **framewise displacement** between sequential transforms
- Classifies motion events based on set motion threshold
- Generates a motion summary figure and livestreams to a web interface (`http://<host IP address>:8080/stream.mjpg`)

---
### ⚡ MOCO (Motion Correction) Feedback
When prospective motion correction is enabled:
1. Queue Processor prioritizes the most recent image data (FIFO off)
2. Registration computes current subject position  
3. Fire Server converts transform into scanner coordinates  
4. A MOCO feedback packet is sent to the scanner (MOCO on)
5. The scanner updates the imaging field-of-view

This creates a **closed-loop control system** that reduces motion artifacts during acquisition.

---
## 📂 Data Flow
All components communicate through a shared filesystem (typically mounted as `/data`).
Future developments will focus on handling all data in working memory.

---
## ⚠️ Design Considerations

### Real-Time Constraints
- MOCO mode must prioritize **low latency** over completeness

### Shared Volume File-Based Communication
- Acts as both data store and inter-process communication layer
- Potential for race conditions (e.g., partial file writes)

### Queue Strategy Tradeoffs
- **FIFO** → first-in-first-out for completeness 
- **Latest-only** → last-in-first-out for lowest latency

---
## ⚙️ Configuration
Runtime behavior is controlled via environment variables:
- `MOCO_FLAG`
- `REG_TYPE`
- `FIFO_FLAG`
- `HEAD_RADIUS`
- `MOTION_THRESH`