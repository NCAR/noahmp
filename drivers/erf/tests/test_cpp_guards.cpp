// ---------------------------------------------------------------------------
// test_cpp_guards -- NoahmpIO_vector sizing guards (Tier 1, WILL_FAIL).
//
// NoahmpIO_vector must be sized exactly once and never with size 0 (the Fortran
// mirror is allocated once; relocation would dangle the self-referential
// pointers). Each violation prints a specific diagnostic and calls NoahmpIO_fatal
// (routed to a clean non-zero _Exit by the linked test_abort_handler.cpp).
//
// Detection is by OUTPUT (PASS/FAIL_REGULAR_EXPRESSION in tests/CMakeLists.txt),
// not exit code: the specific guard message must appear, and the "expected abort
// did NOT occur" fall-through below must not. This is stricter than WILL_FAIL,
// which would treat ANY non-zero exit (e.g. a mistyped selector) as a pass.
//
// Selector argv[1]: double | zero.
// ---------------------------------------------------------------------------
#include <NoahmpIO.H>
#include <cstring>
#include <cstdio>

int main(int argc, char** argv) {
    const char* which = (argc > 1) ? argv[1] : "";

    if (std::strcmp(which, "double") == 0) {
        NoahmpIO_vector vec;
        vec.resize(1, 0);
        vec.resize(1, 0);        // second resize -> abort
    } else if (std::strcmp(which, "zero") == 0) {
        NoahmpIO_vector vec;
        vec.resize(0, 0);        // zero size -> abort
    } else {
        std::fprintf(stderr, "test_cpp_guards: unknown selector '%s'\n", which);
        return 2;
    }

    // Unreachable if the guard fired as expected.
    std::fprintf(stderr, "test_cpp_guards: expected abort did NOT occur for '%s'\n",
                 which);
    return 0;
}
