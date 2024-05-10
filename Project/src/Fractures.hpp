#pragma once

#include<iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace FractureLibrary {

struct Traces {
    vector<array<unsigned int,2>> TraceIDFractures;
    vector<MatrixXd> TraceCoordinates;
};

struct Fractures {
    vector<MatrixXd> FractVertices;
    vector<map<bool,vector<unsigned int>>> listTraces;
    // listTraces[i] accesses the traces of the i-th fracture;
    // listTraces[i][TRUE] access tips traces for i-th fracture (if there are any);
    // listTraces[i][FALSE] access non-tips traces for i-th fracture (if there are any);
    // listTraces[i][True][j] accesses the id of the j-th trace
};

void importFractureList(const string& filepath, Fractures& FractureList);
void computeTraces(Fractures &FractureList, Traces & TracesList);
double computeSquaredDistancePoints(const Vector3d Point1, const Vector3d Point2);
bool checkIntersectionPossibility(const MatrixXd& Fract1_vertices, const MatrixXd& Fract2_vertices);
int findCase(const vector<Vector3d>& A_B_C_D, const vector<double>& betas);
}

