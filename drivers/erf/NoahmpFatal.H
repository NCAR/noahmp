#ifndef NOAHMP_FATAL_H
#define NOAHMP_FATAL_H

// Fatal-error abort, abstracted away from the parallel runtime. Noah-MP must not
// depend on MPI/AMReX, so it calls NoahmpIO_fatal() and the host (ERF) installs a
// handler that propagates the failure across ranks (e.g. amrex::Abort()). With no
// handler installed, NoahmpIO_fatal() prints the message and falls back to std::abort().
//
// Usage (ERF, once during initialization, before any Noah-MP work):
//     NoahmpIO_set_fatal_handler([](const char* msg){ amrex::Abort(msg); });

using FatalHandler = void (*)(const char* msg);

// Install the host application's abort handler. Not thread-safe; call once during
// initialization before any Noah-MP work begins.
void NoahmpIO_set_fatal_handler(FatalHandler handler);

// Terminate after a fatal error. Routes through the installed handler (responsible
// for cross-rank propagation) and never returns. `msg` may be null.
[[noreturn]] void NoahmpIO_fatal(const char* msg);

#endif  // NOAHMP_FATAL_H
