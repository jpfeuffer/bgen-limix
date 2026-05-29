function(bgen_add_test name)
  set(target test_${name})
  add_executable(${target} src/${name}.c)
  target_link_libraries(${target} PRIVATE BGEN::bgen)
  target_compile_options(${target} PRIVATE ${BGEN_C_WARNING_FLAGS})
  target_compile_definitions(
    ${target} PRIVATE $<$<BOOL:${WIN32}>:_CRT_SECURE_NO_WARNINGS>)
  target_compile_definitions(
    ${target} PRIVATE $<$<BOOL:${WIN32}>:_CRT_NONSTDC_NO_DEPRECATE>)

  set(datadir "${CMAKE_CURRENT_BINARY_DIR}/data/")
  target_compile_definitions(${target} PUBLIC "TEST_DATADIR=\"${datadir}\"")

  add_test(
    NAME ${name}
    COMMAND ${target} -E environment
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

  add_custom_command(
    TARGET ${target}
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
            "${CMAKE_CURRENT_BINARY_DIR}/${name}.tmp")

endfunction()

function(bgen_copy filename)
  file(COPY data/${filename} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/data")
endfunction()

bgen_add_test(v12_zlib_layout2)
bgen_add_test(zero_len_chrom_id)
bgen_add_test(pbgen_example)
bgen_add_test(metadata_example)
bgen_add_test(haplotype)
bgen_add_test(complex)
bgen_add_test(index)
bgen_add_test(open_genotype)
bgen_add_test(assert_interface_2)
bgen_add_test(assert_interface_3)
bgen_add_test(assert_interface_3_complex)
bgen_add_test(wrong_bgen_files)
bgen_add_test(wrong_metadata_files)
bgen_add_test(variant_position_overflow)
bgen_add_test(one_million)
bgen_add_test(create_metafile)

bgen_copy(example.matrix)
bgen_copy(example.1bits.bgen)
bgen_copy(example.1bits.bgen.metafile)
bgen_copy(example.14bits.bgen)
bgen_copy(example.14bits.bgen.metafile)
bgen_copy(example.14bits.bgen.metafile.truncated)
bgen_copy(example.14bits.bgen.truncated)
bgen_copy(example.32bits.bgen)
bgen_copy(example.32bits.bgen.metafile)
bgen_copy(example.v11.bgen)

bgen_copy(pbgen_example.pgen)
bgen_copy(zero_len_chrom_id.bgen)
bgen_copy(wrong.metafile)
bgen_copy(haplotypes.bgen)
bgen_copy(complex.23bits.bgen)

bgen_copy(random.bgen)
bgen_copy(random.bgen.metafile)

bgen_copy(roundtrip1.bgen)

# --- S3 integration test (only when S3 support is compiled in) ---
if(BGEN_ENABLE_S3)
  add_executable(test_s3_open src/s3_open.c)
  target_link_libraries(test_s3_open PRIVATE BGEN::bgen)
  target_compile_options(test_s3_open PRIVATE
    $<$<COMPILE_LANGUAGE:C>:${BGEN_C_WARNING_FLAGS}>)
  target_compile_definitions(test_s3_open PUBLIC
    "TEST_DATADIR=\"${CMAKE_CURRENT_BINARY_DIR}/data/\"")

  add_test(NAME s3_open COMMAND test_s3_open)
  set_tests_properties(s3_open PROPERTIES
    ENVIRONMENT "BGEN_TEST_S3_URI=$ENV{BGEN_TEST_S3_URI}")

  # Credential-chain unit test (no network required)
  add_executable(test_credentials_file src/credentials_file.c)
  target_link_libraries(test_credentials_file PRIVATE BGEN::bgen)
  target_compile_options(test_credentials_file PRIVATE
    $<$<COMPILE_LANGUAGE:C>:${BGEN_C_WARNING_FLAGS}>)
  target_compile_definitions(test_credentials_file PRIVATE BGEN_S3_SUPPORT)
  add_test(NAME credentials_file COMMAND test_credentials_file)
endif()
