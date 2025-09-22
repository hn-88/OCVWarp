import numpy as np
from scipy.ndimage import zoom
from PIL import Image

# --- Parameters ---
map_path = "EP_xyuv_1920.map"
input_w, input_h = 3840, 2160   # source resolution
out_w, out_h = 3840, 2160       # output resolution (same as source)

# --- Read Paul Bourke .map file ---
with open(map_path, "r") as f:
    lines = f.readlines()

nx, ny = map(int, lines[1].split())
data = np.array([[float(x) for x in l.split()] for l in lines[2:]])
grid = data.reshape(ny, nx, 5)  # 5 columns: x, y, u, v, weight

# --- Extract normalized u,v ---
u = grid[::-1, :, 2]  # vertical flip
v = 1 - grid[::-1, :, 3]  # flip for ffmpeg remap

# --- Extract weight (fifth column) ---
weight = grid[::-1, :, 4]

# --- Interpolate u,v to desired resolution ---
scale_x = out_w / nx
scale_y = out_h / ny
u_hr = zoom(u, (scale_y, scale_x), order=1)
v_hr = zoom(v, (scale_y, scale_x), order=1)
weight_hr = zoom(weight, (scale_y, scale_x), order=1)

# --- Convert normalized -> integer source pixel coordinates ---
map_x = np.round(u_hr * (input_w - 1)).astype(np.uint16)
map_y = np.round(v_hr * (input_h - 1)).astype(np.uint16)

# --- Save as ASCII PGM (P2) with maxval=65535 ---
def save_pgm_p2(path, arr):
    h, w = arr.shape
    with open(path, "w") as f:
        f.write(f"P2\n{w} {h}\n65535\n")
        for row in arr:
            f.write(" ".join(map(str, row.tolist())) + "\n")

save_pgm_p2("map_x_directp2.pgm", map_x)
save_pgm_p2("map_y_directp2.pgm", map_y)

print("✅ Saved map_x_directp2.pgm and map_y_directp2.pgm as P2 ASCII with maxval=65535")

# --- Save weight as greyscale PNG (0..255) ---
weight_img = (np.clip(weight_hr, 0, 1) * 255).astype(np.uint8)
Image.fromarray(weight_img, mode='L').save("weight_alpha_mask.png")
print("✅ Saved weight_alpha_mask.png as greyscale alpha mask")
