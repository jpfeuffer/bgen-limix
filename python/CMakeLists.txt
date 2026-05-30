find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_core STABLE_ABI NB_STATIC ${CMAKE_CURRENT_SOURCE_DIR}/bindings.cpp)
target_link_libraries(_core PRIVATE BGEN::bgen)
target_include_directories(_core PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(_core PRIVATE ${Python_NumPy_INCLUDE_DIRS})

install(TARGETS _core LIBRARY DESTINATION cbgen)
