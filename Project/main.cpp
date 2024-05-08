#include "Fractures.hpp"
#include "Eigen/Eigen"
#include "UCDUtilities.hpp"
#include <iostream>
using namespace std;
using namespace FractureLibrary;
int main()
{
    Fractures F;
    importFractureList("./DFN/FR3_data.txt", F);
    computeTraces(F);
    return 0;
}
