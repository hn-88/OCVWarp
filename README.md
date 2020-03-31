# OCVWarp
Warping images and videos for planetarium fulldome display using OpenCV.

If OpenCV and CMake are installed on your system, instructions for building:

```
cd build
cmake ..
make OCVWarp.bin
make OCVWarpNorth.bin
```
Initial behaviour and parameters are set using OCVWarp.ini in the build folder.

Keyboard commands are
```
ESC, x or X to exit
u, + or =   to increase angley by 1 degree - the angle seen vertically in output
m, - or _   to decrease angley by 1 degree 
U           to increase angley by 10 degrees
M           to decrease angley by 10 degrees
k, ] or }   to increase anglex by 1 degree - the angle seen horizontally in output
h, [ or {   to decrease anglex by 1 degree - the angle seen horizontally in output
K           to increase anglex by 10 degrees
H           to decrease anglex by 10 degrees
```

