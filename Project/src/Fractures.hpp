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
    vector<vector<unsigned int>> tracesList;
};

void importFractureList(const string& filepath, Fractures& FractureList);
void computeTraces(const Fractures& FractureList);
double computeSquaredDistancePoints(const Vector3d Point1, const Vector3d Point2);
}

