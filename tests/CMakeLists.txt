include(GoogleTest)
find_package(GTest QUIET)
if(NOT GTest_FOUND)
    include(FetchContent)
    FetchContent_Declare(
      googletest
      URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()

if(MSVC)
    add_compile_options(/wd4251)

    option(
        gtest_force_shared_crt
        "Use shared (DLL) run-time lib even when Google Test is built as static lib."
        ON)
    option(gtest_disable_pthreads "Disable uses of pthreads in gtest." ON)

    set(Google_Tests_LIBS
        oldnames.lib
        debug msvcrtd.lib
        debug msvcprtd.lib
        optimized msvcrt.lib
        optimized msvcprt.lib
    )

    # Universal C runtime introduced in VS 2015 (cl version 19)
    if(NOT(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "19"))
        list(APPEND Google_Tests_LIBS
            debug vcruntimed.lib
            debug ucrtd.lib
            debug concrtd.lib
            optimized vcruntime.lib
            optimized ucrt.lib
            optimized concrt.lib
        )
    endif()
endif()

if(WIN32)
    add_definitions(-D_USE_MATH_DEFINES)
endif(WIN32)


# Add test executables here
add_executable(test_run)
target_include_directories(
	test_run PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../OndselSolver
)
target_sources(test_run
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/test.cpp
)
if(GTest_FOUND)
    target_link_libraries(test_run
        GTest::gtest_main
        OndselSolver
    )
else()
    target_link_libraries(test_run
        gtest_main
        ${Google_Tests_LIBS}
        OndselSolver
    )

    include(GoogleTest)
endif()
# discovers tests by asking the compiled test executable to enumerate its tests
set(CMAKE_GTEST_DISCOVER_TESTS_DISCOVERY_MODE PRE_TEST)

gtest_discover_tests(test_run)

add_executable(test_dynamics TestDynamics.cpp)
target_include_directories(test_dynamics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../OndselSolver)
target_link_libraries(test_dynamics PRIVATE GTest::gtest_main OndselSolver)
gtest_discover_tests(test_dynamics DISCOVERY_MODE PRE_TEST PROPERTIES TIMEOUT 120 LABELS dynamics)

# With a DLL-based Google Test installation, allow both test executables to run
# without requiring users to modify PATH. No-op for a static installation.
if(WIN32 AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.21)
    foreach(test_target IN ITEMS test_run test_dynamics)
        add_custom_command(TARGET ${test_target} POST_BUILD
            COMMAND ${CMAKE_COMMAND}
                "-DTEST_RUNTIME_DLLS=$<TARGET_RUNTIME_DLLS:${test_target}>"
                "-DTEST_RUNTIME_DESTINATION=$<TARGET_FILE_DIR:${test_target}>"
                -P "${CMAKE_CURRENT_SOURCE_DIR}/CopyRuntimeDlls.cmake"
            VERBATIM)
    endforeach()
endif()
