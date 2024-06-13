#include "Fractures.hpp"
#include "Polygons.hpp"
#include "UCDUtilities.hpp"
#include<iostream>
#include "Eigen/Eigen"
#include <fstream>
#include <iomanip>
using namespace std;
using namespace Eigen;
namespace GeometryLibrary {
//*********************************************************
double computeSquaredDistancePoints(const Vector3d &Point1, const Vector3d &Point2) {
    double dist = 0.0;
    for(unsigned int i=0;i<3;i++) {
        dist += (Point1(i)-Point2(i))*(Point1(i)-Point2(i));
    }
    return dist;
}
//*********************************************************
    vector<vector<vector<unsigned int>>> TriangulatePolygons(const vector<vector<unsigned int>>& Cell2DCoordinates)
    {
        const unsigned int numPolygons = Cell2DCoordinates.size();
        vector<vector<vector<unsigned int>>> triangleList(numPolygons);

        for(unsigned int p = 0; p < numPolygons; p++)
        {
            const unsigned int numPolygonVertices = Cell2DCoordinates[p].size();

            for (unsigned int v = 0; v < numPolygonVertices; v++)
            {
                const unsigned int nextVertex = Cell2DCoordinates[p][(v + 1) % numPolygonVertices];
                const unsigned int nextNextVertex = Cell2DCoordinates[p][(v + 2) % numPolygonVertices];

                if ((v + 2) % numPolygonVertices == 0)
                    break;

                vector<unsigned int> triangle_vertices = {Cell2DCoordinates[p][0], nextVertex, nextNextVertex};

                triangleList[p].push_back(triangle_vertices);
            }
        }
        return triangleList;
    }
//*********************************************************
    vector<double> computePolygonsArea(const vector<vector<unsigned int>>& listVertices, const vector<Vector3d>& VerticesCoordinates) {

        vector<vector<vector<unsigned int>>> triangleList = TriangulatePolygons(listVertices);

        const unsigned int numPolygons = listVertices.size();
        vector<double> area(numPolygons,0.0);
        for(unsigned int p = 0; p < numPolygons; p++) {
            for(unsigned int t = 0; t < triangleList[p].size(); t++) {
                Matrix3d points = Matrix3d::Zero();
                points.col(0) = VerticesCoordinates[triangleList[p][t][0]];
                points.col(1) = VerticesCoordinates[triangleList[p][t][1]];
                points.col(2) = VerticesCoordinates[triangleList[p][t][2]];

                Triangle triangle(points);

                area[p] += triangle.computeArea();
            }
        }
        return area;

    }
//*********************************************************
    void FractMeshGedimInterface(const vector<vector<unsigned int>>& Cell2DVertices, vector<vector<unsigned int>>& triangles, VectorXi& materials)
    {
        const unsigned int numPolygons = Cell2DVertices.size();
        vector<vector<vector<unsigned int>>> triangleList = TriangulatePolygons(Cell2DVertices);

        unsigned int numTotalTriangles = 0;
        for(unsigned int p = 0; p < numPolygons; p++)
            numTotalTriangles += triangleList[p].size();

        triangles.reserve(numTotalTriangles);
        materials = VectorXi::Zero(numTotalTriangles);

        unsigned int count = 0;
        for(unsigned int p = 0; p < numPolygons; p++)
        {
            for(unsigned int t = 0; t < triangleList[p].size(); t++)
            {
                triangles.push_back(triangleList[p][t]);
                materials(count) = p;
                count++;
            }
        }
    }
//*********************************************************
    bool pointInsidePolygon(double tol, const Vector3d& P, const vector<unsigned int>& VerticesIDs, const vector<Vector3d>& VerticesCoordinates, Vector3d& PlaneNormal) {
        /// CASE WHERE P LIES ON ONE OF THE EDGES IS NOT SUPPORTED
        Vector3d meanPoint = {0.0,0.0,0.0};
        const unsigned int num_vertices = VerticesIDs.size();
        for(unsigned int i=0;i<num_vertices;i++) {
            meanPoint += VerticesCoordinates[VerticesIDs[i]];
        }
        meanPoint = meanPoint/num_vertices;
        Vector3d n = (VerticesCoordinates[VerticesIDs[1]]-VerticesCoordinates[VerticesIDs[0]]).cross(PlaneNormal);
        if((meanPoint-VerticesCoordinates[VerticesIDs[0]]).dot(n)<-tol) {
            PlaneNormal = -PlaneNormal;
        }
        for(unsigned int i=0; i<num_vertices;i++) {
            const Vector3d P0 = VerticesCoordinates[VerticesIDs[i]];
            const Vector3d P1 = VerticesCoordinates[VerticesIDs[(i+1)%num_vertices]];
            n = (P1-P0).cross(PlaneNormal);
            if((P-P0).dot(n)<-tol) {
                return false;
            }
        }
        return true;
    }
//*********************************************************
}
