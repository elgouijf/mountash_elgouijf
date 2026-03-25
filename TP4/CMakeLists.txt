## Indique qu'elle version de cmake minimale est utilisée.
cmake_minimum_required(VERSION 3.16.3)

## Specifie le nom du projet et le langage de programmation du projet
project(TP4 CXX)

## Force le standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## Indique les options de compilations pour le projet
## -Wall : affiche les warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -march=native")

## Définir la localisation des entêtes.
## include_directories sera propagé à l'ensemble du projet
include_directories(${CMAKE_SOURCE_DIR}/include)

add_compile_definitions(PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

## Parcours les sous répertoires contenant les définitions (.cxx)
## On commence par créer une bibliothèque
add_subdirectory(src)

## On créée les tests pour la bibliothèque.
## Ces tests sont unitaires ou fonctionnels
add_subdirectory(test)


