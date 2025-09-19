import numpy as np
from scipy.ndimage import zoom

# --- Parameters ---
map_path = "EP_xyuv_1920.map"
input_w, input_h = 4096, 4096   # input video resolution
out_w, out_h = 3840, 2160       # output remap grid resolution

# --- Read .map file ---
with open(map_path, "r") as f:
    lines = f.readlines()

nx, ny = map(int, lines[1].split())  # 100 x 60
data = np.array([[float(x) for x in l.split()] for l in lines[2:]])
grid = data.reshape(ny, nx, 5)

# Extract u,v
u = grid[:, :, 2]  # normalized 0..1
v = grid[:, :, 3]  # normalized 0..1

# --- Interpolate to 3840x2160 ---
scale_x = out_w / nx
scale_y = out_h / ny
u_hr = zoom(u, (scale_y, scale_x), order=1)
v_hr = zoom(v, (scale_y, scale_x), order=1)

# --- Convert to pixel coords in 4096x4096 input ---
map_x = u_hr * (input_w - 1)
map_y = v_hr * (input_h - 1)

# --- Flip Y (so top is curved, bottom flat) ---
map_y = (input_h - 1) - map_y

# --- Scale to uint16 for PGM ---
map_x_u16 = np.clip(map_x / (input_w - 1) * 65535, 0, 65535).astype(np.uint16)
map_y_u16 = np.clip(map_y / (input_h - 1) * 65535, 0, 65535).astype(np.uint16)

# --- Helper to save PGM (big-endian, P5) ---
def save_pgm_u16(path, array):
    h, w = array.shape
    with open(path, "wb") as f:
        f.write(bytearray(f"P5\n{w} {h}\n65535\n", "ascii"))
        array.byteswap().tofile(f)  # ensure big-endian

# --- Save flipped maps ---
save_pgm_u16("map_x_flipped.pgm", map_x_u16)
save_pgm_u16("map_y_flipped.pgm", map_y_u16)

print("Saved map_x_flipped.pgm and map_y_flipped.pgm")
