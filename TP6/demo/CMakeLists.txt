# Simulation principale
add_executable(simu_tp4 simu_tp4.cxx)
target_link_libraries(simu_tp4 particules)

# Programme 3D
add_executable(test_3d test_3d.cxx)
target_link_libraries(test_3d particules)

# Programme plot UR
add_executable(test_plot_ur test_plot_ur.cxx)
target_link_libraries(test_plot_ur particules)