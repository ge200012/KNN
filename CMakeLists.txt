cmake_minimum_required(VERSION 3.16)
project(gnns)


set(CMAKE_CXX_STANDARD 14)
#自己添加的内容，导入pybind11库
set(CMAKE_CXX_STANDARD_REQUIRED ON)




#set(ENV{PYTHONPATH} "${CMAKE_CURRENT_SOURCE_DIR}")



#include_directories(gnns PRIVATE ${pybind11_INCLUDE_DIR})
#include_directories(gnns PRIVATE ${Python3_INCLUDE_DIRS})
#include_directories(gnns PRIVATE ${PythonLibs_INCLUDE_DIRS})



#set(pybind11_DIR /root/miniconda3/lib/python3.12/site-packages/pybind11/share/cmake/pybind11)
#set(Pyhton3_DIR /root/miniconda3/include/python3.12)
#set(PhthonLibs_DIR /root/miniconda3/lib)

#find_package(pybind11 REQUIRED)

#set(PYTHON_INCLUDE_DIRS /home/KNN /root/miniconda3/include/python3.12)
#set(PYTHON_LIBRARIES /home/KNN /root/miniconda3/lib)

#set(PythonLibs_DIR /root/miniconda3/lib)


set(pybind11_DIR /root/anaconda3/lib/python3.11/site-packages/pybind11/share/cmake/pybind11/)
#message(STATUS "pybind11_DIR set to: ${pybind11_DIR}")
#find_package(Python3 3.12 EXACT REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)
#find_package(PythonLibs REQUIRED)


#find_package(OpenMP REQUIRED)


#set(ENV{PYTHONPATH} "/home/KNN:${CMAKE_CURRENT_SOURCE_DIR}:${CMAKE_CURRENT_BINARY_DIR}")
set(SOURCE_FILES main.cpp constant.h define.h io.h distance.cpp general.h distance.h IndexFilter.cpp IndexFilter.h similarTrajectory.cpp similarTrajectory.h utils.cpp utils.h exactS.cpp exactS.h pss.cpp pss.h experiment.cpp experiment.h mostSimilar.cpp mostSimilar.h)
add_executable(gnns ${SOURCE_FILES})



#target_link_libraries(gnns PRIVATE ${Python3_LIBRARIES} pybind11::pybind11)

#target_link_libraries(gnns PRIVATE pthread dl util m python3.12)
target_link_libraries(gnns PRIVATE pybind11::embed)
#target_link_libraries(gnns PRIVATE OpenMP::OpenMP_CXX)

#自己添加的内容，链接pybind11库和python库
#target_link_libraries(gnns PRIVATE ${Python3_LIBRARIES})


#find_package(PythonLibs 3.12 EXACT REQUIRED)


#include_directories(/root/miniconda3/include/python3.12
# /root/miniconda3/lib/python3.12/site-packages/pybind11/include
# /home/KNN)
#link_directories(/root/miniconda3/lib
# /root/miniconda3/lib/python3.12/site-packages/pybind11/include/pybind11
# /home/KNN)



include_directories(algorithms util io evaluation .)

#target_include_directories(gnns PRIVATE ${Python3_INCLUDE_DIRS})
#target_link_libraries(gnns PRIVATE ${Python3_LIBRARIES})

#set_target_properties(gnns PROPERTIES ENVIRONMENT "PYTHONPATH=${PYTHONPATH}")



find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_compile_options(gnns PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(gnns PRIVATE ${OpenMP_CXX_LIBRARIES})
endif()
















