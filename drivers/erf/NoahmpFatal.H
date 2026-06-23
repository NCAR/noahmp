#ifndef NOAHMP_FATAL_H
#define NOAHMP_FATAL_H

// ---------------------------------------------------------------------------
// Fatal-error abort, abstracted away from the parallel runtime.
//
// Noah-MP must not depend on MPI/AMReX directly, so it never calls MPI_Abort.
// Instead it calls NoahmpIO_fatal(), and the embedding application (ERF) installs
// a handler that propagates the failure across all ranks -- e.g. amrex::Abort(),
// which calls MPI_Abort under the hood. The parallel-environment dependency thus
// lives entirely in the host, not in Noah-MP.
//
// If no handler is installed, NoahmpIO_fatal() prints the message and falls back to
// std::abort() (single-process termination -- correct for a serial build, and a
// safe last resort otherwise).
//
// Usage (ERF, once during initialization, before any Noah-MP work):
//     NoahmpIO_set_fatal_handler([](const char* msg){ amrex::Abort(msg); });
// ---------------------------------------------------------------------------

using FatalHandler = void (*)(const char* msg);

// Install the host application's abort handler. Not thread-safe; call once during
// initialization before any Noah-MP work begins.
void NoahmpIO_set_fatal_handler(FatalHandler handler);

// Terminate after a fatal error. Routes through the installed handler (which is
// responsible for cross-rank propagation) and never returns. `msg` may be null.
[[noreturn]] void NoahmpIO_fatal(const char* msg);

#endif  // NOAHMP_FATAL_H
