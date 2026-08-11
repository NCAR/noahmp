// ---------------------------------------------------------------------------
// test_io_driver_cpp -- full Noah-MP step through the PUBLIC C++ API, in ERF's
// call order.
//
// This mirrors the per-block sequence ERF uses in
//   ../../../Source/LandSurfaceModel/Noah-MP/ERF_NOAHMP_Init.cpp   (cold init)
//   ../../../Source/LandSurfaceModel/Noah-MP/ERF_NOAHMP_Advance.cpp (per step)
// driving NoahmpIO_type / NoahmpIO_vector exactly as the coupler does:
//
//   resize(1,lev) -> blkid/level -> ScalarInitDefault -> rank/comm
//     -> ReadNamelist -> ReadLandHeader -> set domain/tile/memory bounds
//     -> VarInitDefault -> ReadTable -> ReadLandMain -> InitMain -> WriteLand(0)
//   then per step: stage forcing -> itimestep -> DriverMain()
//
// The C++ mirror does not expose the Fortran-only inputs (the wrfinput path in
// erf_setup_file_*, and YR/JULIAN), and ERF supplies those from the namelist /
// its clock; here small bind(C) test shims stand in. The same shims own the MPI
// world and generate the synthetic wrfinput, so this file needs no MPI/NetCDF
// headers. Generic like the Fortran driver test: synthesizes its own fixture, or
// reads any real file via NOAHMP_TEST_WRFINPUT.
//
// Beyond the finite/physical surface-state checks, it also asserts (a) the driver
// actually produced output -- TSK must differ from the cold-init snapshot, or a
// no-op driver would pass every plausibility check -- and (b) C++<->Fortran index
// agreement: an asymmetric pattern written through the C++ views is read back
// through Fortran accessors at the same (i,j)/(i,k,j), catching a transposed or
// offset index map that is invisible when every access goes through one side.
// ---------------------------------------------------------------------------
#include <NoahmpIO.H>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include "test_util.H"

// Test-support shims (Fortran, bind(C)) from test_io_support.F90.
extern "C" {
    int  noahmp_test_mpi_init_c();
    void noahmp_test_mpi_finalize_c();
    void noahmp_test_make_wrfinput_c(const char* path, int plen, int nx, int ny, int nsoil);
    void noahmp_test_set_setup_file_c(int level, int blkid, const char* path, int plen);
    void noahmp_test_set_time_c(int level, int blkid, int yr, noahmp_real julian);
    void noahmp_test_query_dims_c(const char* path, int plen, int* nx, int* ny, int* nsoil);
    // Read Fortran-owned storage directly, to cross-check C++ view index mapping.
    noahmp_real noahmp_test_read_hfx_c (int level, int blkid, int i, int j);
    noahmp_real noahmp_test_read_tslb_c(int level, int blkid, int i, int k, int j);
}

static bool finite_view2d(const NoahmpArray2D<noahmp_real>& a) {
    for (int j = a.begin[1]; j <= a.end[1]; ++j)
        for (int i = a.begin[0]; i <= a.end[0]; ++i)
            if (!std::isfinite(double(a(i, j)))) return false;
    return true;
}

static bool finite_view3d(const NoahmpArray3D<noahmp_real>& a) {
    // axes are (i, layer, j); operator() takes them positionally.
    for (int j = a.begin[2]; j <= a.end[2]; ++j)
        for (int k = a.begin[1]; k <= a.end[1]; ++k)
            for (int i = a.begin[0]; i <= a.end[0]; ++i)
                if (!std::isfinite(double(a(i, k, j)))) return false;
    return true;
}

int main() {
    const int comm = noahmp_test_mpi_init_c();  // Fortran MPI_COMM_WORLD handle

    // ---- pick fixture: synthetic (default) or any real WPS file via env ----
    const char* env = std::getenv("NOAHMP_TEST_WRFINPUT");
    const bool synthetic = (env == nullptr || env[0] == '\0');
    std::string path;
    int nx = 4, ny = 4, nsoil = 4;
    if (synthetic) {
        path = "test_wrfinput_cpp.nc";
        noahmp_test_make_wrfinput_c(path.c_str(), int(path.size()), nx, ny, nsoil);
    } else {
        path = env;
        noahmp_test_query_dims_c(path.c_str(), int(path.size()), &nx, &ny, &nsoil);
        std::fprintf(stderr, "Using external wrfinput: %s  nx=%d ny=%d nsoil=%d\n",
                     path.c_str(), nx, ny, nsoil);
    }

    // ===================== ERF_NOAHMP_Init.cpp order =====================
    NoahmpIO_vector vec;
    vec.resize(1, /*level=*/0);
    NoahmpIO_type& nio = vec[0];

    nio.blkid = 0;
    nio.level = 0;
    nio.ScalarInitDefault();
    nio.rank = 0;
    nio.comm = comm;

    nio.ReadNamelist();
    // ERF asserts the namelist NSOIL matches the fab sizing; here it must match
    // the fixture's soil dimension.
    CHECK_EQ(nio.nsoil, nsoil);

    // The wrfinput path is Fortran-only (namelist-driven in ERF); point it here.
    noahmp_test_set_setup_file_c(0, 0, path.c_str(), int(path.size()));
    nio.ReadLandHeader();

    // Domain / tile / memory bounds span the whole block (kds..kde = 1..2).
    nio.xstart = 0; nio.xend = nx-1; nio.ystart = 0; nio.yend = ny-1;
    nio.ids = 0; nio.ide = nx-1; nio.jds = 0; nio.jde = ny-1; nio.kds = 1; nio.kde = 2;
    nio.its = 0; nio.ite = nx-1; nio.jts = 0; nio.jte = ny-1; nio.kts = 1; nio.kte = 2;
    nio.ims = 0; nio.ime = nx-1; nio.jms = 0; nio.jme = ny-1; nio.kms = 1; nio.kme = 2;

    // Fortran-only time scalars (ERF fills these from its clock).
    noahmp_test_set_time_c(0, 0, 2023, noahmp_real(229.0));

    nio.VarInitDefault();
    nio.ReadTable();
    nio.ReadLandMain();
    nio.InitMain();
    nio.WriteLand(0);   // initial land plotfile (tag 0), as ERF does

    // Sanity: views must span the domain (guards against a vacuous check below).
    CHECK_EQ(nio.TSK.end[0] - nio.TSK.begin[0] + 1, nx);
    CHECK_EQ(nio.TSK.end[1] - nio.TSK.begin[1] + 1, ny);

    // ===================== ERF_NOAHMP_Advance.cpp order ==================
    // Stage benign atmospheric forcing through the C++ views (lowest level; the
    // driver mirrors level 1 -> level 2 internally).
    for (int j = nio.ystart; j <= nio.yend; ++j) {
        for (int i = nio.xstart; i <= nio.xend; ++i) {
            nio.T_PHY(i, 1, j)   = noahmp_real(290.0);
            nio.QV_CURR(i, 1, j) = noahmp_real(0.006);
            nio.U_PHY(i, 1, j)   = noahmp_real(3.0);
            nio.V_PHY(i, 1, j)   = noahmp_real(1.0);
            nio.P8W(i, 1, j)     = noahmp_real(1.0e5);
            nio.SWDOWN(i, j)     = noahmp_real(400.0);
            nio.GLW(i, j)        = noahmp_real(340.0);
            nio.COSZEN(i, j)     = noahmp_real(0.6);
            nio.RAINBL(i, j)     = noahmp_real(0.0);
            nio.SR(i, j)         = noahmp_real(0.0);
            nio.MP_RAINNC(i, j)  = noahmp_real(0.0);
            nio.MP_SNOW(i, j)    = noahmp_real(0.0);
            nio.MP_GRAUP(i, j)   = noahmp_real(0.0);
            nio.MP_HAIL(i, j)    = noahmp_real(0.0);
        }
    }

    // Snapshot TSK to prove the driver produces output: every plausibility check
    // below already holds at cold-init (TSK read from the file at ~290 K, HFX/TSLB
    // finite), so a driver that silently did nothing would pass them all. The
    // post-loop delta catches that.
    std::vector<noahmp_real> tsk0;
    for (int j = nio.TSK.begin[1]; j <= nio.TSK.end[1]; ++j)
        for (int i = nio.TSK.begin[0]; i <= nio.TSK.end[0]; ++i)
            tsk0.push_back(nio.TSK(i, j));

    // Advance two steps: itimestep==1 (initial-guess branch) then a normal step.
    for (int it = 1; it <= 2; ++it) {
        nio.itimestep = it;
        nio.DriverMain();

        CHECK(finite_view2d(nio.TSK));    // no NaN/Inf in surface temperature
        CHECK(finite_view2d(nio.HFX));    // no NaN/Inf in sensible heat flux
        CHECK(finite_view3d(nio.TSLB));   // no NaN/Inf in soil temperature
        bool tsk_ok = true;
        for (int j = nio.TSK.begin[1]; j <= nio.TSK.end[1]; ++j)
            for (int i = nio.TSK.begin[0]; i <= nio.TSK.end[0]; ++i)
                if (!(nio.TSK(i, j) > noahmp_real(150.0) && nio.TSK(i, j) < noahmp_real(400.0)))
                    tsk_ok = false;
        CHECK(tsk_ok);                    // TSK stays physically plausible [K]
        std::fprintf(stderr, "  DriverMain iteration %d completed\n", it);
    }

    // Not a no-op: the driver must have written surface temperature. A driver that
    // returned without touching TSK would leave it bit-identical to the cold-init
    // read and still satisfy every finiteness/range check above.
    bool tsk_changed = false;
    {
        std::size_t idx = 0;
        for (int j = nio.TSK.begin[1]; j <= nio.TSK.end[1]; ++j)
            for (int i = nio.TSK.begin[0]; i <= nio.TSK.end[0]; ++i, ++idx)
                if (nio.TSK(i, j) != tsk0[idx]) tsk_changed = true;
    }
    CHECK(tsk_changed);                   // driver updated TSK (not a no-op)

    // -------- cross-language index agreement (C++ (i,j) == Fortran (i,j)) --------
    // Everything above touches the arrays exclusively through the C++ views, so a
    // transposed or offset index map would be self-consistent and invisible (and
    // the synthetic domain is square, nx==ny, so a transpose would not even change
    // the shape). Prove the C++ view and the Fortran module-global block address
    // the SAME element by writing an ASYMMETRIC pattern through the C++ view and
    // reading it back through a Fortran accessor at the same indices. Done last, on
    // output fields, so it does not perturb the physics above.
    bool hfx_alias_ok = true;             // 2-D
    for (int j = nio.ystart; j <= nio.yend; ++j)
        for (int i = nio.xstart; i <= nio.xend; ++i)
            nio.HFX(i, j) = noahmp_real(1000 * i + j);        // 1000*i+j != 1000*j+i
    for (int j = nio.ystart; j <= nio.yend; ++j)
        for (int i = nio.xstart; i <= nio.xend; ++i)
            if (noahmp_test_read_hfx_c(0, 0, i, j) != noahmp_real(1000 * i + j))
                hfx_alias_ok = false;
    CHECK(hfx_alias_ok);                  // C++ HFX(i,j) aliases Fortran HFX(i,j)

    bool tslb_alias_ok = true;            // 3-D column-major (i, layer, j)
    for (int j = nio.ystart; j <= nio.yend; ++j)
        for (int k = 1; k <= nsoil; ++k)
            for (int i = nio.xstart; i <= nio.xend; ++i)
                nio.TSLB(i, k, j) = noahmp_real(100 * i + 10 * k + j);
    for (int j = nio.ystart; j <= nio.yend; ++j)
        for (int k = 1; k <= nsoil; ++k)
            for (int i = nio.xstart; i <= nio.xend; ++i)
                if (noahmp_test_read_tslb_c(0, 0, i, k, j) !=
                    noahmp_real(100 * i + 10 * k + j))
                    tslb_alias_ok = false;
    CHECK(tslb_alias_ok);                 // C++ TSLB(i,k,j) aliases Fortran TSLB(i,k,j)

    noahmp_test_mpi_finalize_c();
    TEST_SUMMARY("test_io_driver_cpp");
}
