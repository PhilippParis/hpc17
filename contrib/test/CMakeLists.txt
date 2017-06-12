add_definitions(-std=c99)

add_executable(test_coll_implem
${CMAKE_SOURCE_DIR}/contrib/test/test_helpers/buffer_handling.c
${CMAKE_SOURCE_DIR}/contrib/test/test_helpers/test_utils.c
${CMAKE_SOURCE_DIR}/contrib/test/test_helpers/mpi_colls_check.c
${CMAKE_SOURCE_DIR}/contrib/test/verify_coll_implem.c
${CMAKE_SOURCE_DIR}/contrib/collectives/gather_scatter_implem.c
)
TARGET_LINK_LIBRARIES(test_coll_implem ${COMMON_LIBRARIES} )
SET_TARGET_PROPERTIES(test_coll_implem PROPERTIES COMPILE_FLAGS "${MY_COMPILE_FLAGS}"
                                                  C_STANDARD 99)




