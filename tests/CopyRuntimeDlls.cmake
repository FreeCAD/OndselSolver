# TARGET_RUNTIME_DLLS can be empty with a static Google Test installation.
if(TEST_RUNTIME_DLLS)
    file(COPY ${TEST_RUNTIME_DLLS} DESTINATION "${TEST_RUNTIME_DESTINATION}")
endif()
