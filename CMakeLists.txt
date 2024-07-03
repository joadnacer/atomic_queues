cmake_minimum_required(VERSION 3.14)
project(AtomicQueues LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(atomic_queues INTERFACE)
target_include_directories(atomic_queues INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.4.0
)

FetchContent_MakeAvailable(Catch2)

add_executable(atomic_queue_test src/atomic_queue_test.cpp)
target_link_libraries(atomic_queue_test atomic_queues Catch2::Catch2WithMain)

add_custom_target(run_tests
    COMMAND atomic_queue_test
    DEPENDS atomic_queue_test
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

add_executable(spsc_bench src/spsc_bench.cpp)
target_link_libraries(spsc_bench atomic_queues)
target_compile_options(spsc_bench PRIVATE -O3)

add_executable(mpmc_bench src/mpmc_bench.cpp)
target_link_libraries(mpmc_bench atomic_queues)
target_compile_options(mpmc_bench PRIVATE -O3)

add_custom_target(run_spsc_bench
    COMMAND spsc_bench
    DEPENDS spsc_bench
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

add_custom_target(run_mpmc_bench
    COMMAND mpmc_bench 1 2 4 8 16
    DEPENDS mpmc_bench
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)