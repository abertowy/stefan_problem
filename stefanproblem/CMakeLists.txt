cmake_minimum_required(VERSION 3.5)

project(stefanProblemLevelSet VERSION 0.1 LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt5 REQUIRED COMPONENTS Widgets Charts)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets Charts)

set(PROJECT_SOURCES
        main.cpp
        stefanProblem.cpp
        stefanProblem.h
)

add_executable(stefanProblemLevelSet
    ${PROJECT_SOURCES}
)

target_link_libraries(stefanProblemLevelSet PUBLIC
    Qt${QT_VERSION_MAJOR}::Widgets
    Qt${QT_VERSION_MAJOR}::Charts)

set_target_properties(stefanProblemLevelSet PROPERTIES
    WIN32_EXECUTABLE TRUE
)
