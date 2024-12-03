#ifndef NOAHMP_IO_H
#define NOAHMP_IO_H

#include <NoahArray.H>

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
    void NoahmpInitMain(NoahmpIO_struct* ptr);
}

struct NoahmpIO_arrays {
        NoahArray2D<double> XLAT;
        NoahArray2D<double> WSLAKEXY;
   };

extern void NoahmpIOVarInit(NoahmpIO_struct* sptr, NoahmpIO_arrays* aptr);

#endif // NOAHMP_IO_H
