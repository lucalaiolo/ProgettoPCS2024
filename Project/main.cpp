#include "Fractures.hpp"
#include "Eigen/Eigen"
#include "UCDUtilities.hpp"
#include <iostream>
using namespace std;
using namespace DFNLibrary;
int main()
{
    Fractures F;
    Traces T;
    importFractureList("./DFN/FR200_data.txt", F, T);

    exportParaview("./DFN200.inp",F);

    return 0;
}
