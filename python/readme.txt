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
  [3:v]format=gray,scale=3840:2160[mask]; 
  [remapped][mask]alphamerge[out]
" -map "[out]" output.png


