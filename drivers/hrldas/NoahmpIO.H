#ifndef NOAHMP_IO_H
#define NOAHMP_IO_H

extern "C" {
    struct NoahmpIO_struct {
        int ids, ide, jds, jde, kds, kde;
        int ims, ime, jms, jme, kms, kme;
        int its, ite, jts, jte, kts, kte;
    };

    // Function to expose the struct to Fortran
    void NoahmpIOVarInitDefault(NoahmpIO_struct** ptr);
}

#endif // NOAHMP_IO_H
