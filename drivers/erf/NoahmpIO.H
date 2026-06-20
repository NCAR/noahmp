#ifndef NOAHMP_IO_H
#define NOAHMP_IO_H

#include <NoahArray.H>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <cstddef>
#include <type_traits>

extern "C" {
    struct NoahmpIO_type_fi {

        int *ids, *ide, *jds, *jde, *kds, *kde;
        int *ims, *ime, *jms, *jme, *kms, *kme;
        int *its, *ite, *jts, *jte, *kts, *kte;
        int *xstart, *xend, *ystart, *yend;
        int *nsoil, *nsnow;
        int *itimestep, *ntime;
        int *rank, *blkid, *level;
        int *comm;

        noahmp_real* DTBL = nullptr;       // Timestep
        noahmp_real* ZLVL = nullptr;       // Reference height [m]

        noahmp_real* XLAT = nullptr;       // latitude [rad]
        noahmp_real* WSLAKEXY = nullptr;   // lake water storage [mm]
        noahmp_real* U_PHY = nullptr;
        noahmp_real* T_PHY = nullptr;
        noahmp_real* V_PHY = nullptr;
        noahmp_real* QV_CURR = nullptr;
        noahmp_real* HFX = nullptr;
        noahmp_real* LH = nullptr;
        noahmp_real* SWDOWN = nullptr;
        noahmp_real* GLW = nullptr;
        noahmp_real* TSK = nullptr;
        noahmp_real* EMISS = nullptr;
        noahmp_real* ALBSFCDIRXY = nullptr;
        noahmp_real* ALBSFCDIFXY = nullptr;
        noahmp_real* COSZEN = nullptr;
        noahmp_real* P8W = nullptr;
        noahmp_real* TAU_EW = nullptr;
        noahmp_real* TAU_NS = nullptr;

        NoahmpIO_type_fi(        
          int* ids, int* ide, int* jds, int* jde, int* kds, int* kde,
          int* ims, int* ime, int* jms, int* jme, int* kms, int* kme,
          int* its, int* ite, int* jts, int* jte, int* kts, int* kte,
          int* xstart, int* xend, int* ystart, int* yend, int* itimestep,
          int* ntime, int* nsoil, int* nsnow, int* rank,
          int* blkid, int* level, int *comm, noahmp_real* DTBL, noahmp_real* ZLVL
        ):
          ids(ids), ide(ide), jds(jds), jde(jde), kds(kds), kde(kde),
          ims(ims), ime(ime), jms(jms), jme(jme), kms(kms), kme(kme),
          its(its), ite(ite), jts(jts), jte(jte), kts(kts), kte(kte),
          xstart(xstart), xend(xend), ystart(ystart), yend(yend), itimestep(itimestep),
          ntime(ntime), nsoil(nsoil), nsnow(nsnow), rank(rank), blkid(blkid), 
          level(level), comm(comm), DTBL(DTBL), ZLVL(ZLVL) {}
    };
}

// ---------------------------------------------------------------------------
// ABI guards
//
// NoahmpIO_type_fi mirrors the bind(C) Fortran derived type of the same name
// in NoahmpIO_fi.F90. Every member on BOTH sides is an opaque handle (a C
// pointer / type(C_PTR)), so a correct C<->Fortran mapping requires exactly
// the same NUMBER of members, in the same ORDER, with no padding. No single
// check covers all three; the guards are layered, and each catches what the
// others cannot:
//
//   * count + packing  -- static_assert below (compile time). Standard-layout
//     and sizeof == N*sizeof(void*) pins the member count and rules out padding.
//   * total byte size  -- NoahmpIOTypeFiSize_fi() (run time). Confirms the
//     Fortran compiler sees the same struct size; catches a member added/removed
//     on only one side.
//   * element precision -- NoahmpRealSize_fi() (run time). The struct is all
//     pointers, so its size is identical for float vs double; this is the only
//     way to detect a noahmp_real vs c_kind_noahmp mismatch (DOUBLE_PREC applied
//     to only one compiler).
//   * member ORDER     -- NoahmpIOCheckMemberOrder_fi() (run time). A pure
//     reordering has the same count, size, and precision, yet silently crosses
//     the whole pointer map; only a positional handshake catches it.
//
// The run-time guards are all run once via NoahmpIO_AssertAbi(). If you add or
// reorder a variable, update NoahmpIO.H, NoahmpIO_fi.F90, NOAHMP_IO_FI_NUM_MEMBERS,
// and the canonical order list in NoahmpIOCheckMemberOrder_fi together.
// ---------------------------------------------------------------------------
inline constexpr std::size_t NOAHMP_IO_FI_NUM_MEMBERS = 50;

static_assert(std::is_standard_layout<NoahmpIO_type_fi>::value,
              "NoahmpIO_type_fi must be standard-layout to share an ABI with the "
              "Fortran bind(C) NoahmpIO_type_fi.");
static_assert(sizeof(NoahmpIO_type_fi) == NOAHMP_IO_FI_NUM_MEMBERS * sizeof(void*),
              "NoahmpIO_type_fi layout drifted from the Fortran bind(C) type: "
              "unexpected member count or padding. Update NoahmpIO.H, "
              "NoahmpIO_fi.F90, and NOAHMP_IO_FI_NUM_MEMBERS together.");

// Returns sizeof(NoahmpIO_type_fi) as seen by the Fortran compiler.
extern "C" std::size_t NoahmpIOTypeFiSize_fi();

// Returns sizeof(real(c_kind_noahmp)) as seen by the Fortran compiler.
extern "C" std::size_t NoahmpRealSize_fi();

// Positional member-order handshake: see NoahmpIOCheckMemberOrder_fi in
// NoahmpIO_fi.F90. Returns the canonical index of the first mismatched member,
// or -1 if every member agrees.
extern "C" int NoahmpIOCheckMemberOrder_fi(NoahmpIO_type_fi* probe);

// Runs all run-time ABI guards (size, precision, member order) exactly once,
// aborting with a diagnostic on any mismatch. Cheap to call repeatedly. Invoked
// from both NoahmpIO_type construction (single-instance path) and
// NoahmpIO_vector::resize (block-array path) so neither path can skip the check.
void NoahmpIO_AssertAbi();

class NoahmpIO_type {

    public:
        int ids, ide, jds, jde, kds, kde;
        int ims, ime, jms, jme, kms, kme;
        int its, ite, jts, jte, kts, kte;
        int xstart, xend, ystart, yend;
        int nsoil, nsnow;
        int itimestep, ntime;
        int rank, blkid, level;
        int comm;

        noahmp_real DTBL;
        noahmp_real ZLVL = -9999.0;    // Reference height [m]; set via C++ API or falls back to namelist

        NoahArray2D<noahmp_real> XLAT;
        NoahArray2D<noahmp_real> WSLAKEXY;
        NoahArray3D<noahmp_real> U_PHY;
        NoahArray3D<noahmp_real> T_PHY;
        NoahArray3D<noahmp_real> V_PHY;
        NoahArray3D<noahmp_real> QV_CURR;
        NoahArray2D<noahmp_real> HFX;
        NoahArray2D<noahmp_real> LH;
        NoahArray2D<noahmp_real> SWDOWN;
        NoahArray2D<noahmp_real> GLW;
        NoahArray2D<noahmp_real> TSK;
        NoahArray2D<noahmp_real> EMISS;
        NoahArray3D<noahmp_real> ALBSFCDIRXY;
        NoahArray3D<noahmp_real> ALBSFCDIFXY;
        NoahArray2D<noahmp_real> COSZEN;
        NoahArray3D<noahmp_real> P8W;
        NoahArray2D<noahmp_real> TAU_EW, TAU_NS;

        NoahmpIO_type():
          fptr(&ids, &ide, &jds, &jde, &kds, &kde,
               &ims, &ime, &jms, &jme, &kms, &kme,
               &its, &ite, &jts, &jte, &kts, &kte,
               &xstart, &xend, &ystart, &yend, &itimestep,
               &ntime, &nsoil, &nsnow, &rank, &blkid,
               &level, &comm, &DTBL, &ZLVL) { NoahmpIO_AssertAbi(); }

        // fptr stores the addresses of THIS object's own scalar members, so the
        // struct is self-referential. A default memberwise copy/move would leave
        // the new object's fptr pointing at the SOURCE's members -- dangling once
        // the source dies, and silent memory corruption thereafter. Copying is
        // therefore forbidden. Move cannot simply be deleted because
        // std::vector<NoahmpIO_type>::resize requires the type be move-
        // constructible; the move ctor below re-points fptr at the destination.
        NoahmpIO_type(const NoahmpIO_type&)            = delete;
        NoahmpIO_type& operator=(const NoahmpIO_type&) = delete;
        NoahmpIO_type& operator=(NoahmpIO_type&&)      = delete;

        NoahmpIO_type(NoahmpIO_type&& o) noexcept:
          ids(o.ids), ide(o.ide), jds(o.jds), jde(o.jde), kds(o.kds), kde(o.kde),
          ims(o.ims), ime(o.ime), jms(o.jms), jme(o.jme), kms(o.kms), kme(o.kme),
          its(o.its), ite(o.ite), jts(o.jts), jte(o.jte), kts(o.kts), kte(o.kte),
          xstart(o.xstart), xend(o.xend), ystart(o.ystart), yend(o.yend),
          nsoil(o.nsoil), nsnow(o.nsnow),
          itimestep(o.itimestep), ntime(o.ntime),
          rank(o.rank), blkid(o.blkid), level(o.level),
          comm(o.comm),
          DTBL(o.DTBL), ZLVL(o.ZLVL),
          XLAT(o.XLAT), WSLAKEXY(o.WSLAKEXY),
          U_PHY(o.U_PHY), T_PHY(o.T_PHY), V_PHY(o.V_PHY), QV_CURR(o.QV_CURR),
          HFX(o.HFX), LH(o.LH), SWDOWN(o.SWDOWN), GLW(o.GLW), TSK(o.TSK), EMISS(o.EMISS),
          ALBSFCDIRXY(o.ALBSFCDIRXY), ALBSFCDIFXY(o.ALBSFCDIFXY),
          COSZEN(o.COSZEN), P8W(o.P8W),
          TAU_EW(o.TAU_EW), TAU_NS(o.TAU_NS),
          fptr(o.fptr)
        {
          // fptr(o.fptr) copied all handles: the array handles (set later by the
          // Fortran side via C_LOC into Fortran-owned memory) are correct as-is;
          // the self-referential scalar handles must be re-pointed from the
          // source's members to THIS object's members.
          fptr.ids=&ids; fptr.ide=&ide; fptr.jds=&jds; fptr.jde=&jde; fptr.kds=&kds; fptr.kde=&kde;
          fptr.ims=&ims; fptr.ime=&ime; fptr.jms=&jms; fptr.jme=&jme; fptr.kms=&kms; fptr.kme=&kme;
          fptr.its=&its; fptr.ite=&ite; fptr.jts=&jts; fptr.jte=&jte; fptr.kts=&kts; fptr.kte=&kte;
          fptr.xstart=&xstart; fptr.xend=&xend; fptr.ystart=&ystart; fptr.yend=&yend;
          fptr.nsoil=&nsoil; fptr.nsnow=&nsnow;
          fptr.itimestep=&itimestep; fptr.ntime=&ntime;
          fptr.rank=&rank; fptr.blkid=&blkid; fptr.level=&level;
          fptr.comm=&comm;
          fptr.DTBL=&DTBL; fptr.ZLVL=&ZLVL;
        }

        void ScalarInitDefault();
        void VarInitDefault();
        void InitMain();
        void ReadNamelist();
        void ReadTable();
        void ReadLandHeader();
        void ReadLandMain();
        void DriverMain();
        void WriteLand(int filenum);
        void WriteRestart(const std::string& dir);
        void ReadRestart(const std::string& dir);

    private:
        NoahmpIO_type_fi fptr;
};

// ---------------------------------------------------------------------------
// NoahmpIO_vector
//
// Mirrors a Fortran allocatable array (NoahmpIO_vect(level)%NoahmpIO) that is
// sized exactly ONCE via NoahmpIOTypeVectInit_fi. Two invariants make most of
// std::vector's mutating API unsafe here:
//   1. The C++ vector size and the Fortran allocatable extent must stay equal;
//      only resize(size, level) keeps them in sync.
//   2. Each element's fptr stores pointers into its own members, so any
//      operation that relocates elements (grow/reserve/insert/copy) leaves
//      dangling C pointers that later corrupt memory.
//
// Inheritance is therefore PRIVATE: the mutating half of std::vector's API is
// not part of this type's interface and cannot be reached -- not even through a
// std::vector<NoahmpIO_type>& upcast, which is ill-formed for an inaccessible
// base. Only the size-once entry point (resize(size, level)) and a read/element
// -access subset are re-exported. Any unsupported operation (raw resize,
// push_back, insert, erase, swap, clear, reserve, assign, ...) is a COMPILE
// error rather than a silent relocation -- stronger than a run-time guard.
// ---------------------------------------------------------------------------
class NoahmpIO_vector : private std::vector<NoahmpIO_type> {

    using base = std::vector<NoahmpIO_type>;
    bool initialized_ = false;

public:
    NoahmpIO_vector() = default;

    // Copying/moving the whole container would relocate elements and break the
    // self-referential fptr pointers; disallow at compile time.
    NoahmpIO_vector(const NoahmpIO_vector&)            = delete;
    NoahmpIO_vector& operator=(const NoahmpIO_vector&) = delete;
    NoahmpIO_vector(NoahmpIO_vector&&)                 = delete;
    NoahmpIO_vector& operator=(NoahmpIO_vector&&)      = delete;

    // The one supported sizing entry point; errors if called more than once.
    void resize(size_t size, size_t level);

    // Read-only / element-access subset that is safe once the array is sized.
    // (Element access and iteration do not relocate elements.)
    using base::operator[];
    using base::at;
    using base::front;
    using base::back;
    using base::data;
    using base::begin;
    using base::end;
    using base::cbegin;
    using base::cend;
    using base::size;
    using base::empty;
};

#endif // NOAHMP_IO_H
