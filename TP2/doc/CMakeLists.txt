find_package(Doxygen REQUIRED)

set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/doc/Doxyfile)
set(DOXYGEN_OUT ${CMAKE_BINARY_DIR}/docs)

add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_IN}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
    COMMENT "Generation de la documentation Doxygen"
    VERBATIM
)