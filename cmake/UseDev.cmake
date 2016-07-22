function(eth_apply TARGET REQUIRED SUBMODULE)

	if (${SUBMODULE} STREQUAL "base")
		if (UNIX AND NOT APPLE)
			target_link_libraries(${TARGET} pthread)
		endif()
	endif()

endfunction()
