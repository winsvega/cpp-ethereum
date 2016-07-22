aux_source_directory(. SRC_LIST)
include_directories(..)

file(GLOB HEADERS "*.h")

add_executable(ethkey ${SRC_LIST} ${HEADERS})

find_package(Eth)

target_link_libraries(ethkey ${Eth_ETHCORE_LIBRARIES})

if (APPLE)
	install(TARGETS ethkey DESTINATION bin)
else()
	eth_install_executable(ethkey)
endif()
