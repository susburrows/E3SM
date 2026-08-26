if(NOT DEFINED PROGRAM OR NOT DEFINED EXPECTED OR NOT DEFINED ENCODED_ARGS)
  message(FATAL_ERROR "expected-failure test is missing required inputs")
endif()

string(REPLACE "|" ";" ARGS "${ENCODED_ARGS}")
execute_process(
  COMMAND "${PROGRAM}" ${ARGS}
  RESULT_VARIABLE RESULT
  OUTPUT_VARIABLE STDOUT
  ERROR_VARIABLE STDERR)

if(RESULT EQUAL 0)
  message(FATAL_ERROR
          "command unexpectedly succeeded\nstdout:\n${STDOUT}\nstderr:\n${STDERR}")
endif()

string(CONCAT COMBINED "${STDOUT}" "\n" "${STDERR}")
string(FIND "${COMBINED}" "${EXPECTED}" MATCH_INDEX)
if(MATCH_INDEX EQUAL -1)
  message(FATAL_ERROR
          "command failed without expected text '${EXPECTED}'"
          "\nstdout:\n${STDOUT}\nstderr:\n${STDERR}")
endif()
