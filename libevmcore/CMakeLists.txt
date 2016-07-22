aux_source_directory(. SRC_LIST)
include_directories(..)

file(GLOB HEADERS "*.h")

add_library(evmcore ${SRC_LIST} ${HEADERS})

install(TARGETS evmcore RUNTIME DESTINATION bin ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include/evmcore)

