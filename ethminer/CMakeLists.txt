aux_source_directory(. SRC_LIST)
include_directories(..)

file(GLOB HEADERS "*.h")

eth_simple_add_executable(ethminer ${SRC_LIST} ${HEADERS})

eth_use(ethminer REQUIRED JsonRpc::Client)

jsonrpcstub_client_create(ethminer minerfarm.json FarmClient ${CMAKE_CURRENT_SOURCE_DIR} FarmClient)
jsonrpcstub_client_create(ethminer minerphonehome.json PhoneHome ${CMAKE_CURRENT_SOURCE_DIR} PhoneHome)

find_package(Eth)

target_link_libraries(ethminer ${Eth_ETHASHSEAL_LIBRARIES})

if (APPLE)
	install(TARGETS ethminer DESTINATION bin)
else()
	eth_install_executable(ethminer)
endif()
