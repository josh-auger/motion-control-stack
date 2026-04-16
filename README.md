# Motion-Control-Stack

The **motion-control-stack** implements a real-time, multi-process pipeline for image-based motion estimation and 
prospective motion correction (MOCO).  

The system operates as a **closed-loop control pipeline**, continuously ingesting image data, estimating motion, and 
optionally feeding corrections back to the scanner.

---
## 🔄 High-Level Pipeline
```mermaid
flowchart LR

%% =========================
%% Data Ingestion
%% =========================
A[MRD Stream] --> B[Fire Server]

B --> B1[Parse MRD]
B1 --> B2[Write .nhdr]
B1 --> B3[Write .raw]
B1 --> B4[Write pointer .txt]

%% =========================
%% Processing Pipeline
%% =========================
B4 --> C[Queue Processor]
C --> C1[Load .nhdr + .raw]
C1 --> C2[Run sms-mi-reg Registration]
C2 --> D[Transform File (.tfm)]

%% =========================
%% Motion Estimation
%% =========================
D --> E[Motion Monitor]
E --> E1[Framewise Displacement]
E1 --> E2[Motion Classification]
E2 --> E3[Live Visualization]

%% =========================
%% MOCO Feedback Loop
%% =========================
D --> F{MOCO Enabled?}

F -->|Yes| G[Fire Server MOCO Handler]
G --> H[Convert to Scanner Frame]
H --> I[Send Feedback Packet]
I --> A

%% =========================
%% Shared State
%% =========================
subgraph S[Shared Data Volume (/data)]
B2
B3
B4
D
end

---
## 🧱 Core Components

### 🔥 Fire Server
- Listens for incoming MRD image data streams (port `9002`)
- Converts streamed data into:
  - NRRD header files (`.nhdr`)
  - Raw image data (`.raw`)
  - Pointer files (`.txt`)

#### When MOCO is enabled:
- Monitors for new transform files (`.tfm`)
- Converts transforms into scanner coordinate frame
- Sends real-time feedback packets to the scanner

---
### 🔁 Queue Processor
- Monitors the data directory for new pointer files (`.txt`)
- Loads referenced image data (`.nhdr` + `.raw`)
- Executes image registration using `sms-mi-reg`
- Outputs alignment transforms (`.tfm`)

#### Queue Behavior Modes
**MOCO OFF**
- First-In-First-Out (FIFO)
- Processes all incoming data

**MOCO ON**
- Latest-only processing
- Drops stale pointer files
- Minimizes latency for real-time correction

---
### 📈 Motion Monitor
- Monitors transform files (`.tfm`)
- Maintains a running sequence of transforms
- Computes **framewise displacement**
- Classifies motion events based on thresholds
- Generates a live-updating motion summary figure
- Streams results to a web interface

---
### ⚡ MOCO (Motion Correction) Feedback
When motion correction is enabled:

1. Queue Processor prioritizes the most recent image data  
2. Registration computes current subject position  
3. Fire Server converts transform into scanner coordinates  
4. A MOCO feedback packet is sent to the scanner  
5. The scanner updates the imaging field-of-view  

This creates a **closed-loop control system** that reduces motion artifacts during acquisition.

---
## 📂 Data Flow
All components communicate through a shared filesystem (typically mounted as `/data`).

---
## ⚠️ Design Considerations

### Real-Time Constraints
- MOCO mode prioritizes **low latency** over completeness

### File-Based Communication
- Potential for race conditions (e.g., partial file writes)

### Queue Strategy Tradeoffs
- **FIFO** → first-in-first-out complete processing  
- **Latest-only** → most-recent-position responsive processing  

### Shared Volume
- Acts as both:
  - Data store
  - Inter-process communication layer

---
## ⚙️ Configuration
Runtime behavior is controlled via environment variables:
- `MOCO_FLAG`
- `REG_TYPE`
- `HEAD_RADIUS`
- `MOTION_THRESH`