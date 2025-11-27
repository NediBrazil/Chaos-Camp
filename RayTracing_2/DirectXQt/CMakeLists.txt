cmake_minimum_required(VERSION 3.16)

project(DirectXQt LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable Qt MOC, RCC, UIC
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# Qt
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

qt_add_executable(DirectXQt
    MANUAL_FINALIZATION
    ${SOURCES}
    ${HEADERS}
    mainwindow.ui
)

target_link_libraries(DirectXQt
    PRIVATE
        Qt6::Widgets
        dxgi
        d3d12
)

qt_finalize_executable(DirectXQt)

if(MSVC)
    set_target_properties(DirectXQt PROPERTIES
        WIN32_EXECUTABLE YES
    )
endif()
