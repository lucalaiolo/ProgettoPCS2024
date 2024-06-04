#include "Fractures.hpp"
#include "Eigen/Eigen"
#include "UCDUtilities.hpp"
#include <iostream>
using namespace std;
using namespace DFNLibrary;
int main()
{
    const double tol = 1e-13;
    Fractures F;
    Traces T;
    importFractureList("./DFN/FR200_data.txt", F);
    computeTraces(F,T,tol);
    /*
    Vector3d P0 = {0.0,0.0,0.0};
    Vector3d P1 = {1.0,0.0,0.0};
    Vector3d P2 = {1.0,1.0,0.0};
    Vector3d P3 = {0.0,1.0,0.0};

    Vector3d Q0 = {0.5,0.0,0.0};
    Vector3d Q1 = {0.5,1.0,0.0};

    Vector3d Q2 = {1.0,0.5,0.0};
    Vector3d Q3 = {0.0,0.5,0.0};

    Vector3d Q4 = {0.0,0.0,0.0};
    Vector3d Q5 = {1.0,0.5,0.0};

    Vector3d Q6 = {1.0,0.5,0.0};
    Vector3d Q7 = {0.75,1.0,0.0};

    F.FractVertices.resize(1);
    MatrixXd frac = MatrixXd::Zero(3,4);
    frac.col(0) = P0;
    frac.col(1) = P1;
    frac.col(2) = P2;
    frac.col(3) = P3;
    F.FractVertices[0] = frac;

    map<bool,vector<unsigned int>> tracce;
    tracce.insert({false,{0,1,2,3}});
    F.listTraces.push_back(tracce);

    MatrixXd t1 = MatrixXd::Zero(3,2);
    MatrixXd t2 = MatrixXd::Zero(3,2);
    MatrixXd t3 = MatrixXd::Zero(3,2);
    MatrixXd t4 = MatrixXd::Zero(3,2);

    t1.col(0) = Q0;
    t1.col(1) = Q1;

    t2.col(0) = Q2;
    t2.col(1) = Q3;

    t3.col(0) = Q4;
    t3.col(1) = Q5;

    t4.col(0) = Q6;
    t4.col(1) = Q7;

    T.TraceCoordinates.push_back(t1);
    T.TraceCoordinates.push_back(t2);
    T.TraceCoordinates.push_back(t3);
    T.TraceCoordinates.push_back(t4);
    */
    exportParaview("./DFN200.inp",F);
    computePolygonalMesh(F,T,tol);
    exportFractureMesh("./mesh.inp",F,4);
    return 0;
}
