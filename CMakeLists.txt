cmake_minimum_required(VERSION 3.24)
project(main)

find_package(SDL2 REQUIRED)
find_package(GLUT REQUIRED)
find_package(OpenGL REQUIRED)

add_executable(main src/main.cpp src/render.cpp src/simulation.cpp src/common.cpp src/logging.cpp)
include_directories(main ${SDL2_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS})
target_link_libraries(main ${SDL2_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})

add_executable(circle src/generators/circle.cpp)
add_executable(inf src/generators/inf.cpp)
