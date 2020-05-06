if [%Platform%] EQU [x64] 
set OPENCV_DIR C:\tools\OpenCV\Build\x64\vc14
else
if [%Platform%] EQU [Win32] 
set OPENCV_DIR C:\tools\OpenCV\Build\x86\vc14
