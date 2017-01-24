cmake_minimum_required(VERSION 3.1)
project(reallocation)

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(SOURCE_PATH "./MatrixReallocation/source")
set(HEADER_PATH "./MatrixReallocation/include")
set(SOURCE_FILES ${SOURCE_PATH}/reallocation.cpp
				 ${SOURCE_PATH}/service.cpp
				 ${SOURCE_PATH}/taskdata.cpp)
set(HEADER_FILES ${HEADER_PATH}/reallocation.h
				 ${HEADER_PATH}/service.h
				 ${HEADER_PATH}/taskdata.h)

include_directories(${HEADER_PATH})

add_library(reallocation SHARED	${SOURCE_FILES}
		   						${HEADER_FILES})

target_compile_features(reallocation PUBLIC cxx_nullptr)
target_compile_options(reallocation PUBLIC -Wall -O3 -mavx -msse4.1)

install(TARGETS reallocation
		LIBRARY
		DESTINATION lib)
install(FILES ${HEADER_PATH}/reallocation.h
		DESTINATION include)

