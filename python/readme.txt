We can use ffmpeg's remap filter to warp with the default map, using the following:

ffmpeg -i fisheye_resized.jpg -i map_x.pgm -i map_y.pgm -lavfi remap warped_output.png
# here, fisheye_resized.jpg, map_x.pgm and map_y.pgm must have the same size as the output warped_output.png.

Resize the input to 3840×2160 automatically inside FFmpeg - 
----------------------------------------------------------

ffmpeg -i 4096.png -i map_x_directp2.pgm -i map_y_directp2.pgm -filter_complex "[0:v]scale=3840:2160[scaled]; [scaled][1:v][2:v]remap[out]" -map "[out]" 4096output.png

Add the mask as in Paul Bourke's .map file - 
-------------------------------------------
(which is the "weight" field in the .map file)

ffmpeg -i 4096.png -i map_x_directp2.pgm -i map_y_directp2.pgm -i weight_alpha_mask.png -filter_complex "
  [0:v]scale=3840:2160[scaled];
  [scaled][1:v][2:v]remap[remapped];
  [3:v]format=gray,scale=3840:2160[mask_gray];
  [mask_gray]colorchannelmixer=rr=1:gg=1:bb=1[mask_rgb];
  [remapped][mask_rgb]blend=all_mode=multiply[out]
" -map "[out]" 4096wmasked.jpg

Adding NVidia codec and copying the audio stream - 
------------------------------------------------

ffmpeg -i input.mp4 -i map_x_directp2.pgm -i map_y_directp2.pgm -i weight_alpha_mask.png -filter_complex "
  [0:v]scale=3840:2160[scaled];
  [scaled][1:v][2:v]remap[remapped];
  [3:v]format=gray,scale=3840:2160,colorchannelmixer=rr=1:gg=1:bb=1[mask_rgb];
  [remapped][mask_rgb]blend=all_mode=multiply[out]
" -map "[out]" -map 0:a -c:v hevc_nvenc -preset p5 -cq 23 -rc vbr -maxrate 15M -bufsize 26M \
-c:a aac -b:a 128k output.mp4

This runs at 13 fps on a desktop running NVidia RTX 1060 as against 2.5 fps using OCVWarp saving to avc1 codec out.

For starting at two minutes and stopping after 10 seconds duration,
ffmpeg -ss 00:02:00 -i input.mp4 \
       -i map_x_directp2.pgm \
       -i map_y_directp2.pgm \
        -i weight_alpha_mask.png \
       -t 10 \
-filter_complex "
 [0:v]scale=3840:2160[scaled];
  [scaled][1:v][2:v]remap[remapped];
  [3:v]format=gray,scale=3840:2160,colorchannelmixer=rr=1:gg=1:bb=1[mask_rgb];
  [remapped][mask_rgb]blend=all_mode=multiply[out]
" -map "[out]" -map 0:a -c:v hevc_nvenc -c:a copy output.mp4


-------------------------------------------------------------------------------------------------------------------------------

ChatGPT also suggested the following below, but they were not very useful.

With optimized static mask - but this also ran at 13 fps
--------------------------

ffmpeg -i input.mp4 -i map_x_directp2.pgm -i map_y_directp2.pgm -loop 1 -i weight_alpha_mask.png -filter_complex "
  [0:v]scale=3840:2160[scaled];
  [scaled][1:v][2:v]remap[remapped];
  [3:v]format=gray,scale=3840:2160,colorchannelmixer=rr=1:gg=1:bb=1[mask_rgb];
  [remapped][mask_rgb]blend=all_mode=multiply[out]
" -map "[out]" -map 0:a -c:v hevc_nvenc -c:a copy output.mp4

With pre-processed mask - "even faster" - 20fps for some parts of the video - but this had the green video output issue.
--------------------------

ffmpeg -loop 1 -i weight_alpha_mask.png \
-vf "format=gray,scale=3840:2160,colorchannelmixer=rr=1:gg=1:bb=1,format=rgb24" \
-t 1 -c:v libx264 -pix_fmt yuv444p mask_preprocessed.mp4

(-t 1 creates a one second video)

ffmpeg -i input.mp4 -i map_x_directp2.pgm -i map_y_directp2.pgm -stream_loop -1 -i mask_preprocessed.mp4 -filter_complex "
  [0:v]scale=3840:2160[scaled];
  [scaled][1:v][2:v]remap[remapped];
  [remapped][3:v]blend=all_mode=multiply[out]
" -map "[out]" -map 0:a -c:v hevc_nvenc -c:a copy output.mp4


Or, use .y4m to avoid having to specify -pix_fmt rgb24 -s 3840x2160 - but this is only 12 fps and has the green channel issue
-------------------------------------------------------------------

ffmpeg -loop 1 -i weight_alpha_mask.png -vf "format=gray,scale=3840:2160,format=yuv444p" -frames:v 1 -pix_fmt yuv444p mask_3840x2160.y4

ffmpeg -i input.mp4 -i map_x_directp2.pgm -i map_y_directp2.pgm \
-stream_loop -1 -i mask_3840x2160.y4m \
-filter_complex "
  [0:v]scale=3840:2160[scaled];
  [scaled][1:v][2:v]remap[remapped];
  [remapped][3:v]blend=all_mode=multiply[out]
" -map "[out]" -map 0:a -c:v hevc_nvenc -c:a copy output.mp4






