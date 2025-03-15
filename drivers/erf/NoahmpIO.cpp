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

      noahmpio->T_PHY = NoahArray3D<double>(fptr->T_PHY, 
                                           {noahmpio->xstart, noahmpio->kds, noahmpio->ystart}, 
                                           {noahmpio->xend, noahmpio->kde, noahmpio->yend});

      noahmpio->U_PHY = NoahArray3D<double>(fptr->U_PHY, 
                                           {noahmpio->xstart, noahmpio->kds, noahmpio->ystart}, 
                                           {noahmpio->xend, noahmpio->kde, noahmpio->yend});

      noahmpio->V_PHY = NoahArray3D<double>(fptr->V_PHY, 
                                           {noahmpio->xstart, noahmpio->kds, noahmpio->ystart}, 
                                           {noahmpio->xend, noahmpio->kde, noahmpio->yend});

      noahmpio->QV_CURR = NoahArray3D<double>(fptr->QV_CURR, 
                                             {noahmpio->xstart, noahmpio->kds, noahmpio->ystart}, 
                                             {noahmpio->xend, noahmpio->kde, noahmpio->yend});

      noahmpio->SHBXY = NoahArray2D<double>(fptr->SHBXY, 
                                            {noahmpio->xstart, noahmpio->ystart}, 
                                            {noahmpio->xend, noahmpio->yend});

      noahmpio->EVBXY = NoahArray2D<double>(fptr->EVBXY, 
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

void NoahmpReadHDRInfo(NoahmpIO_type* noahmpio) {
     NoahmpIO_type_fi* fptr = &noahmpio->fptr;
     NoahmpReadHDRInfo_fi(fptr);
};

void NoahmpReadLandMain(NoahmpIO_type* noahmpio) {
     NoahmpIO_type_fi* fptr = &noahmpio->fptr;
     NoahmpReadLandMain_fi(fptr);
};

void NoahmpDriverMain(NoahmpIO_type* noahmpio) {
     NoahmpIO_type_fi* fptr = &noahmpio->fptr;
     NoahmpDriverMain_fi(fptr);
};
