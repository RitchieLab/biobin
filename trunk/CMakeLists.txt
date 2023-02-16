cmake_minimum_required(VERSION 3.17.3)

project(biobin)

add_definitions(-DPACKAGE_STRING="biobin 2.3.1")
add_definitions(-DPACKAGE_BUGREPORT="Software_RitchieLab@uphs.upenn.edu")


set_property(GLOBAL PROPERTY USE_FOLDERS ON)

add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=0)

find_package(Boost REQUIRED COMPONENTS filesystem program_options system iostreams thread regex)
include_directories(${Boost_INCLUDE_DIR})

set(GSL_DIR /appl/gsl-2.5)
set(GSL_INCLUDE_DIR ${GSL_DIR}/include)
set(GSL_LIBRARY ${GSL_DIR}/lib/libgsl.so)
set(GSL_CBLAS_LIBRARY ${GSL_DIR}/lib/libgslcblas.so)
find_package(GSL REQUIRED)
include_directories(${GSL_DIR} ${GSL_INCLUDE_DIR})

find_package(SQLite3 REQUIRED)


include_directories(src)

file(GLOB_RECURSE src_knowledge src/knowledge/*.cpp src/knowledge/*.h src/knowledge/*.hpp src/knowledge/liftover/*.cpp src/knowledge/liftover/*.h src/knowledge/detail/*.hpp)

add_library(knowledge ${src_knowledge})

target_link_libraries(knowledge ${Boost_LIBRARIES})

file(GLOB_RECURSE src_biobin src/biobin/*.cpp src/biobin/*.h src/biobin/util/*.cpp src/biobin/util/*.h src/biobin/tests/*.cpp src/biobin/tests/*.h src/biobin/tests/detail/*.cpp src/biobin/tests/detail/*.h)

add_executable(${PROJECT_NAME} ${src_biobin})

target_link_libraries(${PROJECT_NAME} PUBLIC knowledge ${GSL_LIBRARIES} sqlite3 ${Boost_LIBRARIES}               
)


