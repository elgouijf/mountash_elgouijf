cmake_minimum_required(VERSION 3.20)
project(TP6 CXX)

# =====================
# Standard C++
# =====================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# =====================
# OpenMP
# Nécessaire pour utiliser #pragma omp parallel for dans le code C++.
# Les flags de compilation/link sont propagés aux targets qui linkent OpenMP::OpenMP_CXX.
# =====================
find_package(OpenMP REQUIRED)

# =====================
# Flags de compilation
# -O3           : optimisations maximales
# -march=native : instructions SIMD spécifiques au CPU (AVX2, SSE4, etc.)
# -funroll-loops: dépliage des petites boucles (utile pour les boucles x/y/z)
# -ffast-math   : optimisations flottantes agressives (réordonnancement, approximations)
# -flto         : Link Time Optimization (optimisation inter-fichiers)
# -DNDEBUG      : désactive les assert() en mode release
# =====================
add_compile_options(
    -Wall
    -O3
    -march=native
    -funroll-loops
    -ffast-math
    -flto
    -DNDEBUG
)

# =====================
# Includes globaux
# =====================
include_directories(${CMAKE_SOURCE_DIR}/include)

# Rend PROJECT_SOURCE_DIR accessible dans le code C++ (utilisé dans trouver_script_python)
add_compile_definitions(PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# =====================
# GoogleTest
# =====================
include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.17.0
)

FetchContent_MakeAvailable(googletest)

# Active l'infrastructure de tests CTest
enable_testing()

# =====================
# Sous-projets
# =====================
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(doc)
add_subdirectory(demo)