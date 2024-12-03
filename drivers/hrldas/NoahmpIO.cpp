#include <NoahmpIO.H>
#include <NoahArray.H>

void NoahmpIOVarInit(NoahmpIO_struct* sptr, NoahmpIO_arrays* aptr) {
      NoahmpIOVarInitDefault(sptr); 

      aptr->XLAT = NoahArray2D<double>(sptr->XLAT, {sptr->xstart, sptr->ystart}, {sptr->xend, sptr->yend});
      aptr->WSLAKEXY = NoahArray2D<double>(sptr->WSLAKEXY, {sptr->xstart, sptr->ystart}, {sptr->xend, sptr->yend});
};
