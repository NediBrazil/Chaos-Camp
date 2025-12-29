cmake_minimum_required(VERSION 3.16)

project(DirectXQt LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

find_package(Qt6 REQUIRED COMPONENTS Widgets Gui)

set(SOURCES
    main.cpp
    MainWindow.cpp
    ViewportWidget.cpp
    DXRenderer.cpp
)

set(HEADERS
    MainWindow.h
    ViewportWidget.h
    DXRenderer.h
)

include_directories(${CMAKE_SOURCE_DIR}/external/d3dx12)

qt_add_executable(DirectXQt
    ${SOURCES}
    ${HEADERS}
    mainwindow.ui
)

target_link_libraries(DirectXQt
    PRIVATE
        Qt6::Widgets
        Qt6::Gui
        d3d12.lib
        dxgi.lib
        d3dcompiler.lib
        dxguid.lib
)

qt_finalize_executable(DirectXQt)

if(MSVC)
    set_target_properties(DirectXQt PROPERTIES
        WIN32_EXECUTABLE YES
    )
endif()
