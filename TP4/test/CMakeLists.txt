# Dossier de tests
# Créer un executable à partir de fichiers sources
add_executable(test_plot_ur test_plot_ur.cxx)
add_executable(test_3d test_3d.cxx)

target_compile_definitions(test_plot_ur PRIVATE PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
target_compile_definitions(test_3d PRIVATE PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")


target_link_libraries(test_3d particules)
# Fait le lien avec la bibliothèque créée précédemment
# target_link_libraries(test_point particules)
# target_link_libraries(test_univers particules)
# target_link_libraries(test_interactions particules)