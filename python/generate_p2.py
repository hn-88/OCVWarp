import numpy as np
from scipy.ndimage import zoom

# --- Parameters ---
map_path = "EP_xyuv_1920.map"
input_w, input_h = 3840, 2160   # source resolution
out_w, out_h = 3840, 2160       # output resolution (same as source)

# --- Read Paul Bourke .map file ---
with open(map_path, "r") as f:
    lines = f.readlines()

nx, ny = map(int, lines[1].split())
data = np.array([[float(x) for x in l.split()] for l in lines[2:]])
grid = data.reshape(ny, nx, 5)

# Extract normalized u,v
#u = grid[:, :, 2]
#v = grid[:, :, 3] results in warping with curvature at bottom when using ffmpeg remap
#
# Extract u,v (normalized 0..1), BUT flip vertically
u = grid[::-1, :, 2]  # reverse row order
v = grid[::-1, :, 3]


# Orientation fix for ffmpeg remap
v = 1 - v

# Flip horizontally (to fix left-right inversion)
# u = 1 - u


# --- Interpolate u,v to desired resolution ---
scale_x = out_w / nx
scale_y = out_h / ny
u_hr = zoom(u, (scale_y, scale_x), order=1)
v_hr = zoom(v, (scale_y, scale_x), order=1)

# --- Convert normalized -> integer source pixel coordinates ---
map_x = np.round(u_hr * (input_w - 1)).astype(np.uint16)
map_y = np.round(v_hr * (input_h - 1)).astype(np.uint16)

# --- Flip Y so top is curved, bottom flat ---
#map_y = (input_h - 1) - map_y

# --- Save as ASCII PGM (P2) with declared maxval=65535 ---
def save_pgm_p2(path, arr):
    h, w = arr.shape
    with open(path, "w") as f:
        f.write(f"P2\n{w} {h}\n65535\n")
        # Write values row by row
        for row in arr:
            f.write(" ".join(map(str, row.tolist())) + "\n")

save_pgm_p2("map_x_directp2.pgm", map_x)
save_pgm_p2("map_y_directp2.pgm", map_y)

print("✅ Saved map_x_directp2.pgm and map_y_directp2.pgm as P2 ASCII with maxval=65535 (direct coords)")
