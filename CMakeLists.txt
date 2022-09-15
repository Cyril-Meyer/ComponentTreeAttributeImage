cmake_minimum_required(VERSION 3.5)

project(ComponentTreeAttributeImage LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(ComponentTreeAttributeImage
    preprocess_nenist.cpp
#    Algorithms/ComponentTree.h
#    Algorithms/ComponentTree.hxx
#    Algorithms/Morphology.h
#    Algorithms/Morphology.hxx
#    Common/FlatSE.h
#    Common/FlatSE.hxx
#    Common/Image.h
#    Common/Image.hxx
#    Common/ImageIO.hxx
#    Common/ImageIterators.h
#    Common/Point.h
#    Common/Types.h
    )

target_include_directories(ComponentTreeAttributeImage PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Common/)
target_include_directories(ComponentTreeAttributeImage PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Algorithms/)
