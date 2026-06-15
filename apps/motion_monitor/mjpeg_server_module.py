# mjpeg_server.py
import threading
import time
import numpy as np
import cv2
from datetime import datetime
from zoneinfo import ZoneInfo
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

BOUNDARY = "frameboundary"
latest_jpeg = None
lock = threading.Lock()

def generate_dummy_frame(width=1600, height=900, text="Stream connected. Awaiting frames..."):
    """
    Create a white background JPEG frame with centered black text.
    """
    img = 255 * np.ones((height, width, 3), dtype=np.uint8)
    font = cv2.FONT_HERSHEY_SIMPLEX
    font_scale = 1.1
    thickness = 2
    text_size, _ = cv2.getTextSize(text, font, font_scale, thickness)
    text_x = (width - text_size[0]) // 2
    text_y = (height + text_size[1]) // 2
    cv2.putText(img,text,(text_x, text_y),font,font_scale,(0, 0, 0),thickness,cv2.LINE_AA,)
    success, jpeg = cv2.imencode(".jpg", img, [int(cv2.IMWRITE_JPEG_QUALITY), 80])
    if not success:
        raise RuntimeError("Failed to generate dummy JPEG frame")
    return jpeg.tobytes()

def initialize_stream_frame():
    global latest_jpeg
    with lock:
        latest_jpeg = generate_dummy_frame()

def update_frame(jpeg_bytes):
    global latest_jpeg
    with lock:
        latest_jpeg = jpeg_bytes

def add_banner_to_frame(img):
    """
    Add banner text to the bottom of an existing image while preserving all plotted data.
    """
    output = img.copy()
    height, width = output.shape[:2]
    banner_height = 70

    # Create overlay for semi-transparent banner
    overlay = output.copy()
    cv2.rectangle(overlay, (0, height - banner_height), (width, height), (40, 40, 40), thickness=-1)
    # Blend overlay
    alpha = 0.70
    cv2.addWeighted(overlay, alpha, output, 1.0 - alpha, 0, output)

    # Timestamp
    timestamp = datetime.now(ZoneInfo("America/New_York")).strftime("%d-%b-%Y %H:%M:%S %Z")

    # Main message
    line1 = f"MOTION MONITOR RESET : {timestamp}"
    line2 = "Awaiting New Sequence..."


    font = cv2.FONT_HERSHEY_SIMPLEX
    scale = 0.80
    thickness = 2

    (w1, h1), _ = cv2.getTextSize(line1, font, scale, thickness)
    (w2, h2), _ = cv2.getTextSize(line2, font, scale, thickness)

    # Position text within the thinner banner
    y1 = height - 38
    y2 = height - 12

    cv2.putText(
        output,
        line1,
        ((width - w1) // 2, y1),
        font,
        scale,
        (255, 255, 255),
        thickness,
        cv2.LINE_AA
    )

    cv2.putText(
        output,
        line2,
        ((width - w2) // 2, y2),
        font,
        scale,
        (220, 220, 220),
        thickness,
        cv2.LINE_AA
    )

    return output

class MJPEGHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        if self.path != "/stream.mjpg":
            self.send_error(404)
            return

        self.send_response(200)
        self.send_header("Cache-Control", "no-cache")
        self.send_header("Pragma", "no-cache")
        self.send_header(
            "Content-Type",
            f"multipart/x-mixed-replace; boundary={BOUNDARY}"
        )
        self.end_headers()

        try:
            while True:
                with lock:
                    frame = latest_jpeg

                if frame is None:
                    time.sleep(0.01)
                    continue

                self.wfile.write(
                    b"--" + BOUNDARY.encode() + b"\r\n"
                    b"Content-Type: image/jpeg\r\n"
                    b"Content-Length: " + str(len(frame)).encode() + b"\r\n\r\n"
                    + frame + b"\r\n"
                )
                time.sleep(1 / 2)  # client-side pacing only

        except (BrokenPipeError, ConnectionResetError):
            pass

def start_server(port):
    initialize_stream_frame()
    server = ThreadingHTTPServer(("0.0.0.0", port), MJPEGHandler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server
