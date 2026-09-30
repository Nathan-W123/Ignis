# Regenerate Version.hpp at build time, not just at configure time.
#
# `git describe` used to run once, when the build directory was configured, so
# a directory configured at one commit went on stamping every later binary --
# and every result those binaries wrote -- with that commit.  A result that
# names the wrong version is worse than one that names none.
#
# configure_file() only rewrites its output when the text changes, so this
# costs a recompile of the few files that include the header only when the
# describe string actually moves, not on every build.
set(IGNIS_GIT_DESCRIBE "unknown")
if(GIT_EXECUTABLE AND EXISTS "${ROOT}/.git")
  execute_process(COMMAND ${GIT_EXECUTABLE} describe --always --dirty
                  WORKING_DIRECTORY ${ROOT}
                  OUTPUT_VARIABLE IGNIS_GIT_DESCRIBE
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_QUIET)
  if(NOT IGNIS_GIT_DESCRIBE)
    set(IGNIS_GIT_DESCRIBE "unknown")
  endif()
endif()
configure_file(${SRC} ${DST} @ONLY)
