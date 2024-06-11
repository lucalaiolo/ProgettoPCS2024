#pragma once

#include <gtest/gtest.h>
#include "Fractures.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include "UCDUtilities.hpp"
#include "Polygons.hpp"
using namespace Eigen;
using namespace std;
using namespace DFNLibrary;
using namespace GeometryLibrary;
TEST(DFNTest, Mesh) {
    const double tol = 1e-13;
    Fractures F;
    F.FractMesh.resize(1);
    Traces T;
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

    Vector3d Q8 = {0.25,0.625,0.0};
    Vector3d Q9 = {0.25,0.80,0.0};

    Vector3d Q10 = {0.0,0.0,0.0};
    Vector3d Q11 = {0.25,0.4,0.0};

    Vector3d Q12 = {0.5,0.75,0.0};
    Vector3d Q13 = {0.25,1.0,0.0};


    Vector3d normal = {0.0,0.0,1.0};
    double d = 0.0;
    Plane frac_plane(normal,d);
    F.FractPlanes.push_back(frac_plane);
    F.FractVertices.resize(1);
    MatrixXd frac = MatrixXd::Zero(3,4);
    frac.col(0) = P0;
    frac.col(1) = P1;
    frac.col(2) = P2;
    frac.col(3) = P3;
    F.FractVertices[0] = frac;

    map<bool,vector<unsigned int>> tracce;
    tracce.insert({false,{0,1,2,3}});
    tracce.insert({true, {4,5,6}});
    F.listTraces.push_back(tracce);

    MatrixXd t1 = MatrixXd::Zero(3,2);
    MatrixXd t2 = MatrixXd::Zero(3,2);
    MatrixXd t3 = MatrixXd::Zero(3,2);
    MatrixXd t4 = MatrixXd::Zero(3,2);
    MatrixXd t5 = MatrixXd::Zero(3,2);
    MatrixXd t6 = MatrixXd::Zero(3,2);
    MatrixXd t7 = MatrixXd::Zero(3,2);
    t1.col(0) = Q0;
    t1.col(1) = Q1;

    t2.col(0) = Q2;
    t2.col(1) = Q3;

    t3.col(0) = Q4;
    t3.col(1) = Q5;

    t4.col(0) = Q6;
    t4.col(1) = Q7;

    t5.col(0) = Q8;
    t5.col(1) = Q9;

    t6.col(0) = Q10;
    t6.col(1) = Q11;

    t7.col(0) = Q12;
    t7.col(1) = Q13;

    T.TraceCoordinates.push_back(t1);
    T.TraceCoordinates.push_back(t2);
    T.TraceCoordinates.push_back(t3);
    T.TraceCoordinates.push_back(t4);
    T.TraceCoordinates.push_back(t5);
    T.TraceCoordinates.push_back(t6);
    T.TraceCoordinates.push_back(t7);
    computePolygonalMesh(F,T,tol);
    exportFractureMesh("./mesh_test.inp",F,0);

    int ciao = 0;
    EXPECT_EQ(0,ciao);
}


