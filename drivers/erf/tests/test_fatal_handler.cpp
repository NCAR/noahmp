// ---------------------------------------------------------------------------
// test_fatal_handler -- NoahmpIO_fatal routing (Tier 1).
//
// Noah-MP carries no MPI/AMReX dependency: fatal errors go through
// NoahmpIO_fatal(), which calls a host-installed handler (ERF wires
// amrex::Abort -> MPI_Abort) and, if none is installed, falls back to
// std::abort(). Because NoahmpIO_fatal is [[noreturn]], each path is exercised
// in a forked child and the outcome inspected from the parent:
//   * handler installed  -> child exits via the handler (_Exit(123))
//   * no handler         -> child dies on SIGABRT (serial fallback)
// ---------------------------------------------------------------------------
#include <NoahmpFatal.H>
#include <unistd.h>
#include <sys/wait.h>
#include <cstdlib>
#include "test_util.H"

static void handler_exit123(const char*) { _Exit(123); }

// Run `body` in a child; return its raw waitpid status to the parent.
template <typename F>
static int run_in_child(F body) {
    pid_t pid = fork();
    if (pid == 0) {          // child
        body();
        _Exit(0);            // body was supposed to never return
    }
    int status = 0;
    waitpid(pid, &status, 0);
    return status;
}

int main() {
    // 1. Installed handler is dispatched (and short-circuits std::abort).
    int st = run_in_child([]{
        NoahmpIO_set_fatal_handler(handler_exit123);
        NoahmpIO_fatal("boom");
    });
    CHECK(WIFEXITED(st));
    CHECK_EQ(WEXITSTATUS(st), 123);

    // 2. No handler -> std::abort() fallback (SIGABRT), not a normal exit.
    int st2 = run_in_child([]{
        NoahmpIO_fatal("boom-no-handler");
    });
    CHECK(WIFSIGNALED(st2));
    CHECK_EQ(WTERMSIG(st2), SIGABRT);

    TEST_SUMMARY("test_fatal_handler");
}
