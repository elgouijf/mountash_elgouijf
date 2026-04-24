# Dossier de tests
# Créer un executable à partir de fichiers sources
# add_executable(test_plot_ur test_plot_ur.cxx)
# add_executable(test_3d test_3d.cxx)

# target_compile_definitions(test_plot_ur PRIVATE PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
# target_compile_definitions(test_3d PRIVATE PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")


# target_link_libraries(test_3d particules)


# Exécutable unique de tests unitaires
add_executable(unit_tests
    test_vecteur.cxx
    test_particule.cxx
    test_cellule.cxx
    test_univers.cxx
)

# On lie les tests avec :
# - la bibliothèque du projet
# - GoogleTest
target_link_libraries(unit_tests
    particules
    GTest::gtest_main
)



include(GoogleTest)
gtest_discover_tests(unit_tests)