#include "Fractures.hpp"
#include "Eigen/Eigen"
#include "UCDUtilities.hpp"
#include <iostream>
using namespace std;
using namespace DFNLibrary;
int main()
{
    const double tol = 1e-12;
    Fractures F;
    Traces T;
    importFractureList("./DFN/FR200_data.txt", F);
    computeTraces(F,T,tol);
    exportParaview("./DFN200.inp",F);
    computePolygonalMesh(F,T,tol);
    exportFractureMesh("./mesh.inp",F,4);
    if(!printFractureMesh("./Cell0Ds.csv", "./Cell1Ds.csv", "./Cell2Ds.csv", F)) {
        cerr << "Something went wrong!" << endl;
        return 1;
    }

    return 0;
}
