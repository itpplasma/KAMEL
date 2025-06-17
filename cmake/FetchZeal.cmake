include(FetchContent)

find_program(GFORTRAN_EXECUTABLE gfortran REQUIRED
    DOC "Plain gfortran for building Zeal"
)

FetchContent_Declare(zeal
    URL https://elsevier.digitalcommonsdata.com/public-files/datasets/yc9vv7rwyj/files/be55b15b-6d9c-4f0e-b5f7-09b300cfb806/file_downloaded
    URL_HASH MD5=5450326083c6114b026fd765a9733928
)

# Download & extract Zeal sources, making them available
FetchContent_MakeAvailable(zeal)

# Now `zeal_SOURCE_DIR` is set and sources can be compiled by subprojects
