# Créer un executable à partir de fichiers sources
#add_executable(test_1  test_1.cxx)
set(CMAKE_BUILD_TYPE Release)
add_executable(test_point test_point.cxx)
add_executable(test_univers test_univers.cxx)
add_executable(test_interactions_1 test_interactions_1.cxx)
add_executable(test_interactions_2 test_interactions_2.cxx)
add_executable(test_interactions_3 test_interactions_3.cxx)
add_executable(test_interactions_4 test_interactions_4.cxx)
add_executable(test_interactions_5 test_interactions_5.cxx)

# Fait le lien avec la bibliothèque créée précedemment
# Les entêtes sont automatiquement détectés.
#target_link_libraries(test_1 particules)

target_link_libraries(test_point particules)
target_link_libraries(test_univers particules)
target_link_libraries(test_interactions particules)
target_link_libraries(test_interactions_2 particules)
target_link_libraries(test_interactions_3 particules)
target_link_libraries(test_interactions_4 particules)
target_link_libraries(test_interactions_5 particules)