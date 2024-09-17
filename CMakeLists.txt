cmake_minimum_required(VERSION 3.10)
project(ProcessAudio)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

include_directories(${PROJECT_NAME} "/Users/joeypeterson/Desktop/kissfft-master")
include_directories(${PROJECT_NAME} "/Users/joeypeterson/Desktop/AudioFile")
link_directories(/Users/joeypeterson/Desktop/kissfft-master/build)

# Add your executable
add_executable(sinewave_cpp index.cpp)

# Link the Kiss FFT library
target_link_libraries(sinewave_cpp kissfft-float)