cmake_minimum_required(VERSION 2.6)
project(dIdVsuperc) 

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

add_subdirectory(src)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR} ${GSLCBLAS_INCLUDE_DIR})
set(LIBS ${LIBS} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})

target_link_libraries(${LIBS})
link_libraries(${LIBS})

add_executable(
	dIdVsuperc 
		src/common.h 
		src/constraints.c
		src/fit.c 
		src/functions.c 
		src/globals.c 
		src/init.c 
		src/main.c 
		src/misc.c 
		src/plot.c 
		src/ui.c
) 




