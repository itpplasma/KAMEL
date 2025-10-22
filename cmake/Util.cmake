include(FetchContent)

# Default to disconnected updates so re-configure does not try to reach the network.
if(NOT DEFINED FETCHCONTENT_UPDATES_DISCONNECTED)
    set(FETCHCONTENT_UPDATES_DISCONNECTED ON CACHE BOOL
        "Disable FetchContent update steps during configure (enable to force dependency updates)" )
endif()

function(find_or_fetch DEPENDENCY)
    if(DEFINED ENV{CODE} AND EXISTS $ENV{CODE}/${DEPENDENCY})
        set(${DEPENDENCY}_SOURCE_DIR $ENV{CODE}/${DEPENDENCY})
        message(STATUS "Using ${DEPENDENCY} in $ENV{CODE}/${DEPENDENCY}")
    else()
        set(REPO_URL https://github.com/itpplasma/${DEPENDENCY}.git)
        get_branch_or_main(${REPO_URL} REMOTE_BRANCH)
        message(STATUS "Using ${DEPENDENCY} branch ${REMOTE_BRANCH} from ${REPO_URL}")

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
                FetchContent_GetProperties(${DEPENDENCY})
            endif()
        endif()
    endif()

    if(NOT DEFINED ${DEPENDENCY}_SOURCE_DIR OR NOT EXISTS "${${DEPENDENCY}_SOURCE_DIR}/CMakeLists.txt")
        message(FATAL_ERROR "Unable to locate sources for ${DEPENDENCY}: expected ${${DEPENDENCY}_SOURCE_DIR}")
    endif()

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
