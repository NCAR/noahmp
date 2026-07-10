// ---------------------------------------------------------------------------
// test_abi_guard -- the ABI/precision guards (Tier 1).
//
// The C++ NoahmpIO_type_fi and the Fortran bind(C) mirror must agree on layout
// (flat array of pointers) and on the coupling real precision. Layout/count are
// pinned at compile time by static_asserts in NoahmpIO.H; precision is the one
// hazard only a run-time value check can catch (both structs are all pointers,
// so their byte size is precision-invariant). This test re-checks all three at
// run time against the built library.
// ---------------------------------------------------------------------------
#include <NoahmpIO.H>
#include <NoahmpArray.H>
#include <type_traits>
#include "test_util.H"

// Fortran-side probe: sizeof(real(c_kind_noahmp)) as the Fortran compiler sees it.
extern "C" std::size_t NoahmpRealSize_fi();

int main() {
    // Precision agreement across the two separately-invoked compilers.
    const std::size_t f_real = NoahmpRealSize_fi();
    CHECK(f_real == 4 || f_real == 8);
    CHECK_EQ(f_real, sizeof(noahmp_real));

    // The flat-pointer invariant that the header's static_assert encodes, mirrored
    // at run time for good measure.
    CHECK(std::is_standard_layout<NoahmpIO_type_fi>::value);
    CHECK_EQ(sizeof(NoahmpIO_type_fi),
             NOAHMP_IO_FI_NUM_MEMBERS * sizeof(void*));
    CHECK(NOAHMP_IO_FI_NUM_MEMBERS > 0);

    // Runs the one-time precision guard; must return (not abort) when consistent.
    NoahmpIO_AssertAbi();

    TEST_SUMMARY("test_abi_guard");
}
