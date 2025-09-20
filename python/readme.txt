We can use ffmpeg's remap filter to warp with the default map, using the following:

ffmpeg -i fisheye_resized.jpg -i map_x.pgm -i map_y.pgm -lavfi remap warped_output.png
# here, fisheye_resized.jpg, map_x.pgm and map_y.pgm must have the same size as the output warped_output.png.

Resize the input to 3840×2160 automatically inside FFmpeg - 
----------------------------------------------------------

ffmpeg -i 4096.png -i map_x_directp2.pgm -i map_y_directp2.pgm -filter_complex "[0:v]scale=3840:2160[scaled]; [scaled][1:v][2:v]remap[out]" -map "[out]" 4096output.png

Add the mask as in Paul Bourke's .map file - 
-------------------------------------------
(which is the "weight" field in the .map file)

ffmpeg -i input.png -i map_x.pgm -i map_y.pgm -i edge_mask.pgm \
-filter_complex "
  [0:v]scale=3840:2160[scaled]; 
  [scaled][1:v][2:v]remap[remapped]; 
  [3:v]format=rgb,scale=3840:2160[mask]; 
  [remapped][mask]blend=all_mode=multiply[out]
" -map "[out]" output.png

(
ChatGPT says,
For large images, it’s often more efficient to broadcast a single-channel grayscale mask across RGB using the geq (general equation) filter instead of fully converting it to RGB. This avoids channel duplication overhead. Here’s how you can do it in your pipeline:

ffmpeg -i 4096.png -i map_x_directp2.pgm -i map_y_directp2.pgm -i weight_alpha_mask.png -filter_complex "
  [0:v]scale=3840:2160[scaled]; 
  [scaled][1:v][2:v]remap[remapped]; 
  [3:v]format=gray,scale=3840:2160[mask_gray]; 
  [remapped][mask_gray]geq=r='r(X,Y)*lum(X,Y)/255':g='g(X,Y)*lum(X,Y)/255':b='b(X,Y)*lum(X,Y)/255'[out]
" -map "[out]" 4096wmasked.jpg

Explanation

[3:v]format=gray,scale=3840:2160[mask_gray]

Ensures the mask is single-channel grayscale at the target resolution.

[remapped][mask_gray]geq=...

r='r(X,Y)*lum(X,Y)/255' → multiply the red channel of the remapped image by the grayscale value of the mask

Similarly for green (g) and blue (b) channels

lum(X,Y) is the luminance of the second input (mask_gray) at that pixel

Dividing by 255 normalizes the mask (since PGM/PNG grayscale is 0–255)

[out] → final RGB output with edges faded according to the mask

)


