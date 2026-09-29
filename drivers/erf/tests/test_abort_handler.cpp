// ---------------------------------------------------------------------------
// Auto-installs a fatal handler that terminates via std::_Exit with a non-zero
// code instead of std::abort()'s SIGABRT. Linked ONLY into the expected-abort
// (WILL_FAIL) test executables: CTest's WILL_FAIL inverts a non-zero *exit*
// status reliably across platforms, whereas a signal death ("Subprocess
// aborted") is not treated as an ordinary failure and so is not inverted.
//
// The installer runs from a global constructor (.init_array) before main, so it
// works whether the program's main is C++ or Fortran (NoahmpIO_abort in Fortran
// routes through this same handler via NoahmpIO_fatal_c).
// ---------------------------------------------------------------------------
#include <NoahmpFatal.H>
#include <cstdlib>

namespace {
    void exit_on_fatal(const char*) { std::_Exit(7); }  // non-zero, no signal
    const bool g_installed = []() {
        NoahmpIO_set_fatal_handler(exit_on_fatal);
        return true;
    }();
}
