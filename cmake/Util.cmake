include(FetchContent)

# Default to disconnected updates so re-configure does not try to reach the network.
if(NOT DEFINED FETCHCONTENT_UPDATES_DISCONNECTED)
    set(FETCHCONTENT_UPDATES_DISCONNECTED ON CACHE BOOL
        "Disable FetchContent update steps during configure (enable to force dependency updates)")
endif()

function(find_or_fetch DEPENDENCY)
    string(TOUPPER ${DEPENDENCY} _DEP_UPPER)
    # Cache-variable local-source override: -DLIBNEO_PATH=<dir>
    if(DEFINED ${_DEP_UPPER}_PATH AND EXISTS "${${_DEP_UPPER}_PATH}/CMakeLists.txt")
        set(${DEPENDENCY}_SOURCE_DIR "${${_DEP_UPPER}_PATH}")
        message(STATUS "Using ${DEPENDENCY} from ${${_DEP_UPPER}_PATH} (${_DEP_UPPER}_PATH)")
    else()
        set(REPO_URL https://github.com/itpplasma/${DEPENDENCY}.git)
        # Cache-variable ref override: -DLIBNEO_REF=<branch|tag|sha>
        if(DEFINED ${_DEP_UPPER}_REF AND NOT "${${_DEP_UPPER}_REF}" STREQUAL "")
            set(REMOTE_BRANCH "${${_DEP_UPPER}_REF}")
        else()
            get_branch_or_main(${REPO_URL} REMOTE_BRANCH)
        endif()
        message(STATUS "Using ${DEPENDENCY} ref ${REMOTE_BRANCH} from ${REPO_URL}")

        FetchContent_Declare(
            ${DEPENDENCY}
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
            GIT_REPOSITORY ${REPO_URL}
            GIT_TAG ${REMOTE_BRANCH}
        )
        FetchContent_GetProperties(${DEPENDENCY})
        if(NOT ${DEPENDENCY}_POPULATED)
            if(DEFINED ${DEPENDENCY}_SOURCE_DIR AND EXISTS "${${DEPENDENCY}_SOURCE_DIR}/CMakeLists.txt")
                message(STATUS "Reusing existing ${DEPENDENCY} sources at ${${DEPENDENCY}_SOURCE_DIR}")
            else()
                FetchContent_Populate(${DEPENDENCY})
            endif()
        endif()
    endif()

    if(NOT DEFINED ${DEPENDENCY}_SOURCE_DIR OR NOT EXISTS "${${DEPENDENCY}_SOURCE_DIR}/CMakeLists.txt")
        message(FATAL_ERROR "Unable to locate sources for ${DEPENDENCY}: expected ${${DEPENDENCY}_SOURCE_DIR}")
    endif()

    # Disable tests for fetched dependencies to avoid CTest registering
    # tests that won't be built (since we use EXCLUDE_FROM_ALL)
    # Support both old (LIBNEO_ENABLE_TESTS) and new (LIBNEO_BUILD_TESTING) variable names
    set(LIBNEO_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(LIBNEO_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
    set(LIBNEO_ENABLE_GOLDEN_TESTS OFF CACHE BOOL "" FORCE)

    add_subdirectory(${${DEPENDENCY}_SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/${DEPENDENCY}
        EXCLUDE_FROM_ALL
    )
endfunction()


function(get_branch_or_main REPO_URL REMOTE_BRANCH)
    execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    execute_process(
        COMMAND git ls-remote --heads ${REPO_URL} ${BRANCH}
        OUTPUT_VARIABLE BRANCH_EXISTS
        ERROR_VARIABLE GIT_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(BRANCH_EXISTS)
        set(${REMOTE_BRANCH} ${BRANCH} PARENT_SCOPE)
    else()
        set(${REMOTE_BRANCH} "main" PARENT_SCOPE)
    endif()
endfunction()
