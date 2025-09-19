We can use ffmpeg's remap filter to warp with the default map, using the following:

ffmpeg -i fisheye_resized.jpg -vf "scale=4096:4096" \
  -i map_x_flipped.pgm -i map_y_flipped.pgm \
  -lavfi remap=map_x=map_x_flipped.pgm:map_y=map_y_flipped.pgm \
  warped_output.png

crop/pad the warped output back to 3840×2160 automatically inside FFmpeg - 

ffmpeg -i fisheye_resized.jpg -vf "scale=4096:4096" \
  -i map_x_flipped.pgm -i map_y_flipped.pgm \
  -lavfi "[0:v][1:v][2:v]remap=map_x=map_x_flipped.pgm:map_y=map_y_flipped.pgm,scale=3840:2160" \
  warped_output.jpg
