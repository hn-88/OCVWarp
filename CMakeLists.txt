cmake_minimum_required(VERSION 3.2)
cmake_policy(SET CMP0037 NEW)
PROJECT(OCVWarp)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules/")

# copy-pasting from https://github.com/hn-88/pan2fulldome/blob/main/CMakeLists.txt
# to build only OCVWarp, nothing else.
find_package(OpenCV REQUIRED)

# #include_directories(~/OpenCVLocal/include/opencv4)
include_directories(${OpenCV_INCLUDE_DIRS})
add_executable(OCVWarp.bin OCVWarp.cpp tinyfiledialogs.c)
target_link_libraries(OCVWarp.bin ${OpenCV_LIBS})
