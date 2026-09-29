// ---------------------------------------------------------------------------
// test_sharing -- cross-language zero-copy proof (Tier 1, the crown jewel).
//
// The whole point of the coupling is that ERF (C++) and Noah-MP (Fortran) work
// on the SAME array bytes without copying across the language boundary each
// step. This test proves both directions on real coupled arrays:
//
//   C++ -> Fortran : C++ writes vec[0].XLAT / .SMOIS through NoahmpArray views;
//                    Fortran helpers read NoahmpIO_vect(0)%NoahmpIO(0)%XLAT/SMOIS
//                    (the module-global storage) and return the values.
//   Fortran -> C++ : a Fortran helper writes ...%HFX; C++ reads vec[0].HFX.
//
// SMOIS also exercises 3-D column-major (i, layer, j) index agreement. If the
// two sides did not alias the same memory (e.g. a hidden copy crept in), these
// reads would return stale/garbage values and the checks would fail.
// ---------------------------------------------------------------------------
#include <NoahmpIO.H>
#include "test_util.H"

extern "C" {
    void        noahmp_test_prep_block(int level, int blkid);
    noahmp_real noahmp_test_read_xlat (int level, int blkid, int i, int j);
    noahmp_real noahmp_test_read_smois(int level, int blkid, int i, int k, int j);
    void        noahmp_test_write_hfx (int level, int blkid, int i, int j, noahmp_real v);
}

static constexpr int XS = 1, XE = 3, YS = 1, YE = 2;
static constexpr int KMS = 1, KME = 2, NSOIL = 4, NSNOW = 3, NUMRAD = 2;

int main() {
    NoahmpIO_vector vec;
    vec.resize(1, 0);
    NoahmpIO_type& nio = vec[0];
    nio.xstart = XS; nio.xend = XE; nio.ystart = YS; nio.yend = YE;
    nio.kms = KMS; nio.kme = KME;
    nio.nsoil = NSOIL; nio.nsnow = NSNOW; nio.numrad = NUMRAD; nio.rank = 0;
    nio.ScalarInitDefault();
    noahmp_test_prep_block(0, 0);
    nio.VarInitDefault();

    // --- C++ -> Fortran, 2-D (XLAT) --------------------------------------
    for (int j = YS; j <= YE; ++j)
        for (int i = XS; i <= XE; ++i)
            nio.XLAT(i, j) = noahmp_real(1000 * j + i);
    for (int j = YS; j <= YE; ++j)
        for (int i = XS; i <= XE; ++i)
            CHECK_EQ(noahmp_test_read_xlat(0, 0, i, j), noahmp_real(1000 * j + i));

    // --- C++ -> Fortran, 3-D column-major (SMOIS[i, layer, j]) -----------
    for (int j = YS; j <= YE; ++j)
        for (int k = 1; k <= NSOIL; ++k)
            for (int i = XS; i <= XE; ++i)
                nio.SMOIS(i, k, j) = noahmp_real(100 * j + 10 * k + i);
    for (int j = YS; j <= YE; ++j)
        for (int k = 1; k <= NSOIL; ++k)
            for (int i = XS; i <= XE; ++i)
                CHECK_EQ(noahmp_test_read_smois(0, 0, i, k, j),
                         noahmp_real(100 * j + 10 * k + i));

    // --- Fortran -> C++, 2-D (HFX) --------------------------------------
    for (int j = YS; j <= YE; ++j)
        for (int i = XS; i <= XE; ++i)
            noahmp_test_write_hfx(0, 0, i, j, noahmp_real(7 * i + 13 * j));
    for (int j = YS; j <= YE; ++j)
        for (int i = XS; i <= XE; ++i)
            CHECK_EQ(nio.HFX(i, j), noahmp_real(7 * i + 13 * j));

    TEST_SUMMARY("test_sharing");
}
