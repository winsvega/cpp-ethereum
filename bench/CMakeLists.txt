aux_source_directory(. SRC_LIST)
include_directories(..)
include_directories(../utils)

add_executable(bench ${SRC_LIST})

find_package(Boost QUIET REQUIRED COMPONENTS date_time)
find_package(Dev)

target_link_libraries(bench ${Dev_DEVCORE_LIBRARIES})
target_link_libraries(bench ${Dev_DEVCRYPTO_LIBRARIES})
if (DEFINED MSVC)
	target_link_libraries(bench ${Boost_DATE_TIME_LIBRARIES})
endif()

if (APPLE)
	install(TARGETS bench DESTINATION bin)
else()
	eth_install_executable(bench)
endif()

