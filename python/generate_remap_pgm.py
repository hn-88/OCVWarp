import numpy as np
from scipy.ndimage import zoom

# --- Parameters ---
map_path = "EP_xyuv_1920.map"
input_w, input_h = 4096, 4096   # input video resolution
out_w, out_h = 3840, 2160       # output remap grid resolution

# --- Read .map file ---
with open(map_path, "r") as f:
    lines = f.readlines()

nx, ny = map(int, lines[1].split())  # e.g. 100 x 60
data = np.array([[float(x) for x in l.split()] for l in lines[2:]])
grid = data.reshape(ny, nx, 5)

# Extract u,v (normalized 0..1)
u = grid[:, :, 2]
v = grid[:, :, 3]

# --- Interpolate to output resolution ---
scale_x = out_w / nx
scale_y = out_h / ny
u_hr = zoom(u, (scale_y, scale_x), order=1)
v_hr = zoom(v, (scale_y, scale_x), order=1)

# --- Convert to pixel coords in 4096x4096 input ---
map_x = (u_hr * (input_w - 1)).astype(np.uint16)
map_y = (v_hr * (input_h - 1)).astype(np.uint16)

# --- Flip Y (so top is curved, bottom flat) ---
map_y = (input_h - 1) - map_y

# --- Save PGM files with direct coordinates ---
def save_pgm_direct(path, array, maxval):
    h, w = array.shape
    with open(path, "wb") as f:
        f.write(bytearray(f"P5\n{w} {h}\n{maxval}\n", "ascii"))
        array.tofile(f)

save_pgm_direct("map_x_direct.pgm", map_x, input_w - 1)
save_pgm_direct("map_y_direct.pgm", map_y, input_h - 1)

print("Saved map_x_direct.pgm and map_y_direct.pgm with direct pixel coordinates")
