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
    const Vector3d P0 = {0.0,0.0,0.0};
    const Vector3d P1 = {1.0,0.0,0.0};
    const Vector3d P2 = {1.0,1.0,0.0};
    const Vector3d P3 = {0.0,1.0,0.0};

    const Vector3d Q0 = {0.5,0.0,0.0};
    const Vector3d Q1 = {0.5,1.0,0.0};

    const Vector3d Q2 = {1.0,0.5,0.0};
    const Vector3d Q3 = {0.0,0.5,0.0};


    const Vector3d Q4 = {0.0,0.0,0.0};
    const Vector3d Q5 = {1.0,0.5,0.0};

    const Vector3d Q6 = {1.0,0.5,0.0};
    const Vector3d Q7 = {0.75,1.0,0.0};

    const Vector3d Q8 = {0.25,0.625,0.0};
    const Vector3d Q9 = {0.25,0.80,0.0};

    const Vector3d Q10 = {0.0,0.0,0.0};
    const Vector3d Q11 = {0.25,0.4,0.0};

    const Vector3d Q12 = {0.5,0.75,0.0};
    const Vector3d Q13 = {0.25,1.0,0.0};


    const Vector3d normal = {0.0,0.0,1.0};
    const double d = 0.0;
    const Plane frac_plane(normal,d);
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
    computePolygonalMesh(F,T,tol,tol);
    exportFractureMesh("./mesh_test.inp",F,0);
}
TEST(DFNTest, ImportAndCompute){
    Fractures F;
    Traces T;
    const double tol = 1e-13;
    importFractureList("./DFN/FR3_data.txt", F);

    EXPECT_DOUBLE_EQ(0.0, F.FractVertices[0](0,0));
    EXPECT_DOUBLE_EQ(0.0, F.FractVertices[0](1,0));
    EXPECT_DOUBLE_EQ(0.0, F.FractVertices[0](2,0));

    EXPECT_DOUBLE_EQ(0.5, F.FractMeanPoint[0][0]);
    EXPECT_DOUBLE_EQ(1.0, F.FractPlanes[0].Normal[2]);
    EXPECT_DOUBLE_EQ(0.0, F.FractPlanes[0].d);

    computeTraces(F, T, tol);

    EXPECT_DOUBLE_EQ(1.0, T.TraceLength[0]);

    EXPECT_DOUBLE_EQ(0.8, T.TraceCoordinates[0](0,0));
    EXPECT_DOUBLE_EQ(1.0, T.TraceCoordinates[0](1,0));
    EXPECT_DOUBLE_EQ(0.0, T.TraceCoordinates[0](2,0));

    EXPECT_DOUBLE_EQ(0.8, T.TraceCoordinates[0](0,1));
    EXPECT_DOUBLE_EQ(0.0, T.TraceCoordinates[0](1,1));
    EXPECT_DOUBLE_EQ(0.0, T.TraceCoordinates[0](2,1));

    EXPECT_EQ(0, T.TraceIDFractures[0][0]);
    EXPECT_EQ(1, T.TraceIDFractures[0][1]);

    EXPECT_EQ(1, F.listTraces[0][true][0]);
    EXPECT_EQ(0, F.listTraces[0][false][0]);
}

TEST(DFNTest, TraceOnEdge) {
    Fractures F;
    F.FractMeanPoint.resize(2);
    F.FractMeanPoint[0]={0.5,0.5,0};
    F.FractMeanPoint[1]={1,0.5,-0.5};
    F.FractPlanes.resize(2);
    F.FractPlanes[0].Normal = {0,0,1};
    F.FractPlanes[1].Normal = {1,0,0};
    F.FractPlanes[0].d=0;
    F.FractPlanes[1].d = 1;
    F.FractVertices.resize(2);
    F.FractVertices[0] = MatrixXd::Zero(3,4);
    F.FractVertices[0].col(0) << 0.0,0.0,0.0;
    F.FractVertices[0].col(1) << 1.0,0.0,0.0;
    F.FractVertices[0].col(2) << 1.0,1.0,0.0;
    F.FractVertices[0].col(3) << 0.0,1.0,0.0;
    F.FractVertices[1]=MatrixXd::Zero(3,4);
    F.FractVertices[1].col(0) << 1.0,0.0,0.0;
    F.FractVertices[1].col(1) << 1.0,0.0,-1.0;
    F.FractVertices[1].col(2) << 1.0,0.5,-1.0;
    F.FractVertices[1].col(3) << 1.0,0.5,0.0;
    F.listTraces.resize(2);
    Traces T;
    computeTraces(F,T,1e-12);

    EXPECT_EQ(0,F.listTraces[0][true][0]);
    EXPECT_EQ(0,F.listTraces[1][false][0]);
    Vector3d P0 = {1,0,0};
    Vector3d P1 = {1,0.5,0};
    EXPECT_EQ(P0,T.TraceCoordinates[0].col(0));
    EXPECT_EQ(P1,T.TraceCoordinates[0].col(1));

}


TEST(DFNTest, Findcase){
    const double tol = 1e-13;
    vector<double> Betas0 = {3.0, 1.0, 3.0, 1.0};
    EXPECT_EQ(0, findCase(Betas0, tol));

    vector<double> Betas1 = {1.0, 3.0, 4.0, 6.0};
    EXPECT_EQ(1, findCase(Betas1, tol));

    vector<double> Betas2 = {2.0, 4.0, 1.0, 5.0};
    EXPECT_EQ(2, findCase(Betas2, tol));

    vector<double> Betas3 = {1.0, 3.0, 2.0, 4.0};
    EXPECT_EQ(3, findCase(Betas3, tol));

    vector<double> Betas4 = {2.0, 4.0, 1.0, 3.0};
    EXPECT_EQ(4, findCase(Betas4, tol));

    vector<double> Betas5 = {4.0, 6.0, 1.0, 3.0};
    EXPECT_EQ(5, findCase(Betas5, tol));

    vector<double> Betas6 = {1.0, 5.0, 2.0, 4.0};
    EXPECT_EQ(6, findCase(Betas6, tol));

    vector<double> Betas_1 = {1.0, 3.0, 1.0, 4.0};
    EXPECT_EQ(-1, findCase(Betas_1, tol));

    vector<double> Betas_2 = {1.0, 4.0, 1.0, 3.0};
    EXPECT_EQ(-2, findCase(Betas_2, tol));

    vector<double> Betas_3 = {1.0, 4.0, 2.0, 4.0};
    EXPECT_EQ(-3, findCase(Betas_3, tol));

    vector<double> Betas_4 = {2.0, 4.0, 1.0, 4.0};
    EXPECT_EQ(-4, findCase(Betas_4, tol));

}

TEST(DFNTest, IntersectionPossibility){
    // Intersection
    const double tol = 1e-13;
    const MatrixXd Matrix_1 = (Eigen::MatrixXd(3, 4)<< 0.0, 2.0, 2.0, 0.0,
                               0.0, 0.0, 2.0, 2.0,
                               0.0, 0.0, 0.0, 0.0).finished();
    const MatrixXd Matrix_2 = (Eigen::MatrixXd(3, 4)<< 1.0, 1.0, 1.0, 1.0,
                               0.0, 0.0, 2.0, 2.0,
                               1.0, -1.0, -1.0, 1.0).finished();
    const Vector3d  meanPoint1 = {1.0, 1.0, 0.0};
    const Vector3d  meanPoint2 = {1.0, 1.0, 0.0};
    bool check1 = checkIntersectionPossibility(Matrix_1, Matrix_2, meanPoint1, meanPoint2, tol);
    EXPECT_TRUE(check1);

    // No intersection
    MatrixXd Matrix_3 = (Eigen::MatrixXd(3, 4)<< 0.0, 2.0, 2.0, 0.0,
                         3.0, 3.0, 5.0, 5.0,
                         0.0, 0.0, 0.0, 0.0).finished();
    MatrixXd Matrix_4 = (Eigen::MatrixXd(3, 4)<< 1.0, 1.0, 1.0, 1.0,
                         0.0, 0.0, 2.0, 2.0,
                         1.0, -1.0, -1.0, 1.0).finished();
    const Vector3d  meanPoint3 = {1.0, 4.0, 0.0};
    const Vector3d  meanPoint4 = {1.0, 1.0, 0.0};
    bool check2 = checkIntersectionPossibility(Matrix_3, Matrix_4, meanPoint3, meanPoint4, tol);
    EXPECT_FALSE(check2);
}

TEST(GeometryTest, SquaredDistance){
    const Vector3d Point1 = {2.0, 0.0, 0.0};
    const Vector3d Point2 = {0.0, 0.0, 0.0};
    double d = computeSquaredDistancePoints(Point1, Point2);
    EXPECT_DOUBLE_EQ(4.0, d);
}



