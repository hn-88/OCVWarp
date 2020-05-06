if [%Platform%] EQ [x64] 
set OPENCV_DIR C:\tools\OpenCV\Build\x64\vc14
else
if [%Platform%] EQ [Win32] 
set OPENCV_DIR C:\tools\OpenCV\Build\x86\vc14
