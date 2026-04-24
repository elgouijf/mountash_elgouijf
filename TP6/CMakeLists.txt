cmake_minimum_required(VERSION 3.20)

project(TP6 CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -march=native")

include_directories(${CMAKE_SOURCE_DIR}/include)

add_compile_definitions(PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# --- GoogleTest ---
include(FetchContent) # Pour pouvoir récupérer une dépendance directement à partir de l'internet

# Récupère la bibliothèque GoogleTest depuis GitHub
FetchContent_Declare(
  googletest # nom de la dépendance
  GIT_REPOSITORY https://github.com/google/googletest.git # url
  GIT_TAG v1.17.0 # version
)

FetchContent_MakeAvailable(googletest) # Récupère la dépendance et la rend disponible pour être utilisée dans le projet

# Active l'infrastructure de tests CTest
enable_testing()

# Sous-projets, comme ça on ne lance qu'un seul CMake pour tout le projet depuis la racine
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(doc)
add_subdirectory(demo)