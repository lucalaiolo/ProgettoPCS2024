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
    importFractureList("./DFN/FR10_data.txt", F, T);
    return 0;
}
