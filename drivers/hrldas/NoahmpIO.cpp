#include<NoahmpIO.H>

extern "C" {
    // Function to expose the struct to Fortran
    void NoahmpIOVarInitDefault(NoahmpIO_struct* ptr);
}
