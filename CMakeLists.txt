cmake_minimum_required(VERSION 3.10)
project(SEALExample)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 查找SEAL库
find_package(SEAL 4.1 REQUIRED)


add_executable(External main_external.cpp)

# 链接SEAL库

target_link_libraries(External SEAL::seal)

