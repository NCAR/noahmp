#ifndef NOAHMP_IO_H
#define NOAHMP_IO_H

/*
 * #include <NoahmpArray.H>
 *
 * class NoampIO {
 *
 *   private:
 *       NoampIO_struct fstruct;
 *
 *   public:
 *        NoahmpArray2D<double> XLAT;
 *        NoahmpArray3D<double> DZ8W;
 *        NoahmpArray1D<double> DZS;
 *
 *
 *        NoahmpIO();
 *        Init(ids, ide, jds, jde, kds, kde) {
 *          this->fstruct.ids = ids
 *          this->fstruct.jds = jds
 *
 *          XLAT.allocate(fstruct.ids, fstruct.jds);
 *          this->fstruct.XLAT = XLAT.data;
 *
 *          DZS.allocate();
 *
 *          NoahmpIOVarInitDefault(&this->fstruct);
 *        }
 *
 *        Advance() {
 *
 *        }
 * }
 *
 */

extern "C" {
    struct NoahmpIO_struct {

        int ids, ide, jds, jde, kds, kde;
        int ims, ime, jms, jme, kms, kme;
        int its, ite, jts, jte, kts, kte;

        int xstart, xend, ystart, yend;

        double* XLAT = nullptr;       // latitude [rad]
        double* WSLAKEXY = nullptr;   // lake water storage [mm]

    };

    // Function to expose the struct to Fortran
    void NoahmpIOVarInitDefault(NoahmpIO_struct* ptr);
}

#endif // NOAHMP_IO_H
