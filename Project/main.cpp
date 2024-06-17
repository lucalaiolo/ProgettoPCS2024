#include "Fractures.hpp"
#include "Eigen/Eigen"
#include "UCDUtilities.hpp"
#include <iostream>
using namespace std;
using namespace DFNLibrary;
int main(int argc, char *argv[])
{
    const double tol1D_default = 1e-12;
    const double tol2D_default = 1e-12;

    const double tol1D_user = atof(argv[1]);
    const double tol2D_user = atof(argv[2]);

    const double tol1D = max(tol1D_default, tol1D_user);
    const double tol2D = max(tol2D_default,tol2D_user);

    Fractures F;
    Traces T;
    importFractureList("./DFN/FR200_data.txt", F);
    computeTraces(F,T,tol1D);
    exportParaview("./DFN200.inp",F);
    computePolygonalMesh(F,T,tol1D, tol2D);
    exportFractureMesh("./mesh.inp",F,4);
    if(!printFractureMesh("./Cell0Ds.csv", "./Cell1Ds.csv", "./Cell2Ds.csv", F)) {
        cerr << "Something went wrong!" << endl;
        return 1;
    }
    return 0;
}
