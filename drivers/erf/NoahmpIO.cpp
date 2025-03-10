#include <NoahmpIO.H>
#include <NoahArray.H>

void NoahmpIOVarInitDefault(NoahmpIO_type* noahmpio) {

      NoahmpIO_type_fi* fptr = &noahmpio->fptr;
      NoahmpIOVarInitDefault_fi(fptr);

      noahmpio->XLAT = NoahArray2D<double>(fptr->XLAT, 
                                          {noahmpio->xstart, noahmpio->ystart}, 
                                          {noahmpio->xend, noahmpio->yend});

      noahmpio->WSLAKEXY = NoahArray2D<double>(fptr->WSLAKEXY, 
                                              {noahmpio->xstart, noahmpio->ystart}, 
                                              {noahmpio->xend, noahmpio->yend});

};

void NoahmpInitMain(NoahmpIO_type* noahmpio) {
     NoahmpIO_type_fi* fptr = &noahmpio->fptr;
     NoahmpInitMain_fi(fptr);
};

void NoahmpReadNamelist(NoahmpIO_type* noahmpio) {
     NoahmpIO_type_fi* fptr = &noahmpio->fptr;
     NoahmpReadNamelist_fi(fptr);
};

void NoahmpReadTable(NoahmpIO_type* noahmpio) {
     NoahmpIO_type_fi* fptr = &noahmpio->fptr;
     NoahmpReadTable_fi(fptr);
};
