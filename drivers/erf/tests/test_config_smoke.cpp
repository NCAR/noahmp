// ---------------------------------------------------------------------------
// test_config_smoke -- config parsers through the C++ API (Tier 2).
//
// Drives the real Fortran configuration readers via the public C++ interface,
// using only lightweight, already-shipped inputs (no NetCDF):
//   * ReadNamelist() parses the staged namelist.erf; we assert the coupling
//     scalars it populates (nsoil/nsnow/ZLVL) show up on the C++ side -- which
//     also re-confirms scalar sharing through a real code path.
//   * ReadTable() parses the shipped NoahmpTable.TBL; a Fortran probe returns
//     two well-known general-parameter table scalars and we assert they were
//     filled with their expected physical values (not the -9999 dummy).
//
// CMake stages namelist.erf and NoahmpTable.TBL into the working directory.
// ---------------------------------------------------------------------------
#include <NoahmpIO.H>
#include "test_util.H"

extern "C" noahmp_real noahmp_test_read_table(int level, int blkid, int which);

int main() {
    NoahmpIO_vector vec;
    vec.resize(1, 0);
    NoahmpIO_type& nio = vec[0];
    nio.rank = 0;
    nio.ScalarInitDefault();

    // ---- namelist (NOAHLSM_OFFLINE) -------------------------------------
    nio.ReadNamelist();
    CHECK_EQ(nio.nsoil, 4);          // NSOIL = 4 in the fixture
    CHECK_EQ(nio.nsnow, 4);          // NSNOW = 4
    CHECK_CLOSE(nio.ZLVL, 10.0, 1e-3);

    // ---- table (NoahmpTable.TBL general parameters) --------------------
    nio.ReadTable();
    const noahmp_real zbot  = noahmp_test_read_table(0, 0, 0);  // ZBOT_TABLE
    const noahmp_real csoil = noahmp_test_read_table(0, 0, 1);  // CSOIL_TABLE
    CHECK_CLOSE(zbot,  -8.0,   1e-3);      // ZBOT_DATA  = -8.0
    CHECK_CLOSE(csoil, 2.0e6,  1.0);       // CSOIL_DATA = 2.00E+6
    CHECK(zbot  != noahmp_real(-9999.0));  // not the undefined_real dummy
    CHECK(csoil != noahmp_real(-9999.0));

    TEST_SUMMARY("test_config_smoke");
}
