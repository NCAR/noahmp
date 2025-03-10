#ifndef NOAHMP_IO_H
#define NOAHMP_IO_H

#include <NoahArray.H>

extern "C" {
    struct NoahmpIO_type_fi {

        int *ids, *ide, *jds, *jde, *kds, *kde;
        int *ims, *ime, *jms, *jme, *kms, *kme;
        int *its, *ite, *jts, *jte, *kts, *kte;

        int *xstart, *xend, *ystart, *yend;

        int *nsoil, *nsnow;

        double* XLAT = nullptr;       // latitude [rad]
        double* WSLAKEXY = nullptr;   // lake water storage [mm]

        NoahmpIO_type_fi(        
          int* ids, int* ide, int* jds, int* jde, int* kds, int* kde,
          int* ims, int* ime, int* jms, int* jme, int* kms, int* kme,
          int* its, int* ite, int* jts, int* jte, int* kts, int* kte,
          int* xstart, int* xend, int* ystart, int* yend, int *nsoil,
          int* nsnow
        ):
          ids(ids), ide(ide), jds(jds), jde(jde), kds(kds), kde(kde),
          ims(ims), ime(ime), jms(jms), jme(jme), kms(kms), kme(kme),
          its(its), ite(ite), jts(jts), jte(jte), kts(kts), kte(kte),
          xstart(xstart), xend(xend), ystart(ystart), yend(yend), nsoil(nsoil),
          nsnow(nsnow) {}
    };

    // Function to expose the struct to Fortran
    void NoahmpIOVarInitDefault_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpInitMain_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpReadNamelist_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpReadTable_fi(NoahmpIO_type_fi* noahmpio);
}

struct NoahmpIO_type {

        NoahmpIO_type_fi fptr;

        int ids, ide, jds, jde, kds, kde;
        int ims, ime, jms, jme, kms, kme;
        int its, ite, jts, jte, kts, kte;

        int xstart, xend, ystart, yend;

        int nsoil, nsnow;

        NoahArray2D<double> XLAT;
        NoahArray2D<double> WSLAKEXY;

        NoahmpIO_type():
          fptr(&ids, &ide, &jds, &jde, &kds, &kde,
               &ims, &ime, &jms, &jme, &kms, &kme,
               &its, &ite, &jts, &jte, &kts, &kte,
               &xstart, &xend, &ystart, &yend, &nsoil, &nsnow) {}

   };

extern void NoahmpIOVarInitDefault(NoahmpIO_type* noahmpio);
extern void NoahmpInitMain(NoahmpIO_type* noahmpio);
extern void NoahmpReadNamelist(NoahmpIO_type* noahmpio);
extern void NoahmpReadTable(NoahmpIO_type* noahmpio);

#endif // NOAHMP_IO_H
