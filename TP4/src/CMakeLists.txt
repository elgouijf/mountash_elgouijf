# Créer une bibliothèque à partir d'un ensemble de fichiers de définitions (.cxx) sans main.
# Ces fichiers ne contiennent pas de fonction main
add_library(particules
    vecteur.cxx
    particule.cxx
    cellule.cxx
    univers.cxx
    io.cxx
)

# Créer l'exécutable principal de la simulation
add_executable(simu
    main.cxx
)

target_compile_definitions(simu PRIVATE PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
# Lier l'exécutable principal avec la bibliothèque
target_link_libraries(simu particules)

