add_executable(static_target
  static_target.cpp
)
target_link_libraries(static_target core solver_eigen)
set_target_properties(static_target PROPERTIES OUTPUT_NAME static_target)

add_executable(constant_velocity_target
  constant_velocity_target.cpp
)

set_target_properties(constant_velocity_target PROPERTIES OUTPUT_NAME constant_velocity_target)

target_link_libraries(constant_velocity_target core solver_eigen solver_pcg)
