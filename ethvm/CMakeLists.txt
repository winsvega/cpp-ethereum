aux_source_directory(. SRC_LIST)
include_directories(..)

add_executable(ethvm ${SRC_LIST})

find_package(Eth)

target_link_libraries(ethvm ${Eth_ETHASH_LIBRARIES})
target_link_libraries(ethvm ${Eth_ETHASHSEAL_LIBRARIES})

if (APPLE)
	install(TARGETS ethvm DESTINATION bin)
else()
	eth_install_executable(ethvm)
endif()
