cmake_minimum_required(VERSION 3.15)
project(LoRaCpp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB LORA_SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_executable(lora ${LORA_SRC})

if (MSVC)
    target_compile_options(lora PRIVATE /utf-8 /O2)
else()
    target_compile_options(lora PRIVATE -O2)
endif()
