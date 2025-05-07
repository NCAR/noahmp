#include <NoahmpIO.H>
#include <NoahArray.H>

extern "C" {
    void NoahmpIOScalarInitDefault_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpIOVarInitDefault_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpInitMain_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpReadNamelist_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpReadTable_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpReadLandHeader_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpReadLandMain_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpDriverMain_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpIOTypeVectInit_fi(int* NBlocks);
}

void NoahmpIO_type::ScalarInitDefault() {
     NoahmpIOScalarInitDefault_fi(&fptr);
};


void NoahmpIO_type::InitMain() {
     NoahmpInitMain_fi(&fptr);
};

void NoahmpIO_type::ReadNamelist() {
     NoahmpReadNamelist_fi(&fptr);
};

void NoahmpIO_type::ReadTable() {
     NoahmpReadTable_fi(&fptr);
};

void NoahmpIO_type::ReadLandHeader() {
     NoahmpReadLandHeader_fi(&fptr);
};

void NoahmpIO_type::ReadLandMain() {
     NoahmpReadLandMain_fi(&fptr);
};

void NoahmpIO_type::DriverMain() {
     NoahmpDriverMain_fi(&fptr);
};


void NoahmpIO_type::VarInitDefault() {

      NoahmpIOVarInitDefault_fi(&fptr);

      XLAT     = NoahArray2D<double>(fptr.XLAT,     {xstart,ystart}, {xend,yend});
      WSLAKEXY = NoahArray2D<double>(fptr.WSLAKEXY, {xstart,ystart}, {xend,yend});

      T_PHY   = NoahArray3D<double>(fptr.T_PHY,   {xstart,kds,ystart}, {xend,kde,yend});
      U_PHY   = NoahArray3D<double>(fptr.U_PHY,   {xstart,kds,ystart}, {xend,kde,yend});
      V_PHY   = NoahArray3D<double>(fptr.V_PHY,   {xstart,kds,ystart}, {xend,kde,yend});
      QV_CURR = NoahArray3D<double>(fptr.QV_CURR, {xstart,kds,ystart}, {xend,kde,yend});

      SHBXY = NoahArray2D<double>(fptr.SHBXY, {xstart,ystart}, {xend,yend});
      EVBXY = NoahArray2D<double>(fptr.EVBXY, {xstart,ystart}, {xend,yend});

};

void NoahmpIO_vector::resize(size_t size) {
     std::vector<NoahmpIO_type>::resize(size);
     int _size = size;
     NoahmpIOTypeVectInit_fi(&_size);
};
