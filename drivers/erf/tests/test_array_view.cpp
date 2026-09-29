// ---------------------------------------------------------------------------
// test_array_view -- NoahmpArray{1,2,3}D index math (Tier 1, header-only).
//
// The NoahmpArray views are how C++ reads/writes Fortran-owned arrays with no
// copy, so their Fortran-column-major addressing with arbitrary (non-1) lower
// bounds is a correctness core of the whole coupling. This test lays out a
// plain buffer, wraps it in each view, and checks that element (i,j,k) maps to
// the Fortran linear offset -- independently of the struct/library machinery.
// ---------------------------------------------------------------------------
#include <NoahmpArray.H>
#include <vector>
#include "test_util.H"

int main() {
    // ---- 1-D: data[i - begin] -------------------------------------------
    {
        const int b = 5, e = 9;                 // extent 5
        std::vector<noahmp_real> buf(e - b + 1);
        NoahmpArray1D<noahmp_real> a(buf.data(), b, e);
        for (int i = b; i <= e; ++i) a(i) = noahmp_real(100 + i);
        for (int i = b; i <= e; ++i) {
            CHECK_EQ(buf[i - b], noahmp_real(100 + i));   // physical layout
            CHECK_EQ(a(i), noahmp_real(100 + i));         // view round-trip
        }
    }

    // ---- 2-D: data[(i-b0) + n0*(j-b1)] (column-major) -------------------
    {
        const std::array<int,2> b{2, 3}, e{5, 7};
        const int n0 = e[0] - b[0] + 1;          // 4
        std::vector<noahmp_real> buf((e[0]-b[0]+1) * (e[1]-b[1]+1));
        NoahmpArray2D<noahmp_real> a(buf.data(), b, e);
        for (int j = b[1]; j <= e[1]; ++j)
            for (int i = b[0]; i <= e[0]; ++i)
                a(i, j) = noahmp_real(1000 * j + i);
        for (int j = b[1]; j <= e[1]; ++j)
            for (int i = b[0]; i <= e[0]; ++i) {
                const std::size_t off = (i - b[0]) + n0 * (j - b[1]);
                CHECK_EQ(buf[off], noahmp_real(1000 * j + i));
                CHECK_EQ(a(i, j), noahmp_real(1000 * j + i));
            }
        // const view reads the same storage
        const NoahmpArray2D<noahmp_real> ca(buf.data(), b, e);
        CHECK_EQ(ca(e[0], e[1]), noahmp_real(1000 * e[1] + e[0]));
    }

    // ---- 3-D: data[(i-b0) + n0*((k-b1) + n1*(j-b2))] -------------------
    {
        const std::array<int,3> b{2, 1, 3}, e{4, 3, 5};
        const int n0 = e[0] - b[0] + 1;          // 3
        const int n1 = e[1] - b[1] + 1;          // 3
        std::vector<noahmp_real> buf(std::size_t(n0) * n1 * (e[2]-b[2]+1));
        NoahmpArray3D<noahmp_real> a(buf.data(), b, e);
        for (int j = b[2]; j <= e[2]; ++j)
            for (int k = b[1]; k <= e[1]; ++k)
                for (int i = b[0]; i <= e[0]; ++i)
                    a(i, k, j) = noahmp_real(100*j + 10*k + i);
        for (int j = b[2]; j <= e[2]; ++j)
            for (int k = b[1]; k <= e[1]; ++k)
                for (int i = b[0]; i <= e[0]; ++i) {
                    const std::size_t off =
                        (i - b[0]) + n0 * ((k - b[1]) + std::size_t(n1) * (j - b[2]));
                    CHECK_EQ(buf[off], noahmp_real(100*j + 10*k + i));
                    CHECK_EQ(a(i, k, j), noahmp_real(100*j + 10*k + i));
                }
    }

    TEST_SUMMARY("test_array_view");
}
