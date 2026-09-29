// ---------------------------------------------------------------------------
// test_lifecycle -- NoahmpIO_type / NoahmpIO_vector lifecycle (Tier 1).
//
// Drives the documented cold-init order through the public C++ API on a tiny
// block: resize(once) -> set dims -> ScalarInitDefault -> VarInitDefault. The
// key assertion is that every NoahmpArray view comes back with bounds equal to
// the C++-owned dimension scalars: that can only happen if the Fortran allocate
// read those scalars through the wired pointers, so it proves scalar wiring +
// array-view construction end to end (no external data files needed).
// ---------------------------------------------------------------------------
#include <NoahmpIO.H>
#include "test_util.H"

// Fortran test helper: set the IOPT_* branch scalars ReadNamelist would normally
// set, so VarInitDefault's allocation branches are well-defined.
extern "C" void noahmp_test_prep_block(int level, int blkid);

// Small representative geometry (deliberately non-square, lower bounds != 0).
static constexpr int XS = 1, XE = 4, YS = 1, YE = 3;
static constexpr int KMS = 1, KME = 3, NSOIL = 4, NSNOW = 3, NUMRAD = 2;

static void set_dims(NoahmpIO_type& nio) {
    nio.xstart = XS; nio.xend = XE; nio.ystart = YS; nio.yend = YE;
    nio.kms = KMS;   nio.kme = KME;
    nio.nsoil = NSOIL; nio.nsnow = NSNOW; nio.numrad = NUMRAD;
    nio.rank = 0;
}

int main() {
    NoahmpIO_vector vec;
    CHECK(vec.empty());
    vec.resize(1, /*level=*/0);
    CHECK_EQ(vec.size(), std::size_t(1));
    CHECK(!vec.empty());

    NoahmpIO_type& nio = vec[0];
    set_dims(nio);
    nio.ScalarInitDefault();      // wire scalar pointers back to these C++ members
    noahmp_test_prep_block(0, 0); // safe IOPT_* (stand-in for ReadNamelist)
    nio.VarInitDefault();         // allocate Fortran arrays + build C++ views

    // 2-D views: bounds must equal {xstart,ystart}..{xend,yend}.
    CHECK_EQ(nio.XLAT.begin[0], XS);  CHECK_EQ(nio.XLAT.end[0], XE);
    CHECK_EQ(nio.XLAT.begin[1], YS);  CHECK_EQ(nio.XLAT.end[1], YE);
    CHECK_EQ(nio.HFX.begin[0],  XS);  CHECK_EQ(nio.HFX.end[0],  XE);
    CHECK_EQ(nio.HFX.begin[1],  YS);  CHECK_EQ(nio.HFX.end[1],  YE);

    // 3-D soil view SMOIS[xstart:xend, 1:nsoil, ystart:yend].
    CHECK_EQ(nio.SMOIS.begin[0], XS);    CHECK_EQ(nio.SMOIS.end[0], XE);
    CHECK_EQ(nio.SMOIS.begin[1], 1);     CHECK_EQ(nio.SMOIS.end[1], NSOIL);
    CHECK_EQ(nio.SMOIS.begin[2], YS);    CHECK_EQ(nio.SMOIS.end[2], YE);

    // 3-D atmos view U_PHY[xstart:xend, kms:kme, ystart:yend] (i,k,j order).
    CHECK_EQ(nio.U_PHY.begin[0], XS);    CHECK_EQ(nio.U_PHY.end[0], XE);
    CHECK_EQ(nio.U_PHY.begin[1], KMS);   CHECK_EQ(nio.U_PHY.end[1], KME);
    CHECK_EQ(nio.U_PHY.begin[2], YS);    CHECK_EQ(nio.U_PHY.end[2], YE);

    // 3-D radiation view ALBSFCDIRXY[xstart:xend, 1:numrad, ystart:yend].
    CHECK_EQ(nio.ALBSFCDIRXY.begin[1], 1);
    CHECK_EQ(nio.ALBSFCDIRXY.end[1],   NUMRAD);

    // Views point at real storage: a write/read round-trip must not fault.
    nio.XLAT(XS, YS)          = noahmp_real(0.5);
    nio.SMOIS(XE, NSOIL, YE)  = noahmp_real(0.25);
    CHECK_EQ(nio.XLAT(XS, YS),         noahmp_real(0.5));
    CHECK_EQ(nio.SMOIS(XE, NSOIL, YE), noahmp_real(0.25));

    TEST_SUMMARY("test_lifecycle");
}
