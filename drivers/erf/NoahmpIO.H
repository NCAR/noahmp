#ifndef NOAHMP_IO_H
#define NOAHMP_IO_H

#include <NoahArray.H>
#include <vector>

extern "C" {
    struct NoahmpIO_type_fi {

        int *ids, *ide, *jds, *jde, *kds, *kde;
        int *ims, *ime, *jms, *jme, *kms, *kme;
        int *its, *ite, *jts, *jte, *kts, *kte;
        int *xstart, *xend, *ystart, *yend;
        int *nsoil, *nsnow;
        int *itimestep, *ntime;
        char *llanduse;
        int *rank, *blkid;

        double* XLAT = nullptr;       // latitude [rad]
        double* WSLAKEXY = nullptr;   // lake water storage [mm]
        double* U_PHY = nullptr;
        double* T_PHY = nullptr;
        double* V_PHY = nullptr;
        double* QV_CURR = nullptr;
        double* SHBXY = nullptr;
        double* EVBXY = nullptr;

        NoahmpIO_type_fi(        
          int* ids, int* ide, int* jds, int* jde, int* kds, int* kde,
          int* ims, int* ime, int* jms, int* jme, int* kms, int* kme,
          int* its, int* ite, int* jts, int* jte, int* kts, int* kte,
          int* xstart, int* xend, int* ystart, int* yend, int* itimestep,
          int* ntime, int *nsoil, int* nsnow, char* llanduse, int* rank, int* blkid
        ):
          ids(ids), ide(ide), jds(jds), jde(jde), kds(kds), kde(kde),
          ims(ims), ime(ime), jms(jms), jme(jme), kms(kms), kme(kme),
          its(its), ite(ite), jts(jts), jte(jte), kts(kts), kte(kte),
          xstart(xstart), xend(xend), ystart(ystart), yend(yend), itimestep(itimestep),
          ntime(ntime), nsoil(nsoil), nsnow(nsnow), llanduse(llanduse), rank(rank), blkid(blkid) {}
    };
}

class NoahmpIO_type {

    public:
        int ids, ide, jds, jde, kds, kde;
        int ims, ime, jms, jme, kms, kme;
        int its, ite, jts, jte, kts, kte;
        int xstart, xend, ystart, yend;
        int nsoil, nsnow;
        int itimestep, ntime;
        char llanduse[257] = "";
        int rank, blkid;

        NoahArray2D<double> XLAT;
        NoahArray2D<double> WSLAKEXY;
        NoahArray3D<double> U_PHY;
        NoahArray3D<double> T_PHY;
        NoahArray3D<double> V_PHY;
        NoahArray3D<double> QV_CURR;
        NoahArray2D<double> SHBXY;
        NoahArray2D<double> EVBXY;

        NoahmpIO_type():
          fptr(&ids, &ide, &jds, &jde, &kds, &kde,
               &ims, &ime, &jms, &jme, &kms, &kme,
               &its, &ite, &jts, &jte, &kts, &kte,
               &xstart, &xend, &ystart, &yend, &itimestep,
               &ntime, &nsoil, &nsnow, &llanduse[0], &rank, &blkid) {}

        void ScalarInitDefault();
        void VarInitDefault();
        void InitMain();
        void ReadNamelist();
        void ReadTable();
        void ReadLandHeader();
        void ReadLandMain();
        void DriverMain();

    private:
        NoahmpIO_type_fi fptr;
};

class NoahmpIO_vector : public std::vector<NoahmpIO_type> {
public:
    void resize(size_t size);

};

#endif // NOAHMP_IO_H
