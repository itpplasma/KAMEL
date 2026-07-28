add_executable(test_kilca_lorentz_torque test_lorentz_torque.cpp)
target_link_libraries(test_kilca_lorentz_torque PRIVATE kilca_lib)

add_test(
    NAME kilca_lorentz_torque
    COMMAND test_kilca_lorentz_torque
)
