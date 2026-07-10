// ---------------------------------------------------------------------------
// test_cpp_guards -- NoahmpIO_vector sizing guards (Tier 1, WILL_FAIL).
//
// NoahmpIO_vector must be sized exactly once and never with size 0 (the Fortran
// mirror is allocated once; relocation would dangle the self-referential
// pointers). Each violation calls NoahmpIO_fatal -> std::abort (no handler),
// terminating non-zero. Registered WILL_FAIL, so the abort = pass; a normal
// return (exit 0) would flag a guard that stopped firing.
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
