#pragma once
#include<iostream>
#include "Eigen/Eigen"
using namespace std;
using namespace Eigen;

namespace GeometryLibrary {
//*********************************************************
///
/// \brief The Triangle class
///
struct Triangle {

    Matrix3d Vertices;
    Triangle(const MatrixXd& Vertices): Vertices(Vertices) {}
    ///
    /// \brief computeArea
    /// \return area of the triangle
    ///
    double computeArea() {
        Vector3d AB = Vertices.col(1)-Vertices.col(0);
        Vector3d AC = Vertices.col(2)-Vertices.col(0);
        double area = (AB.cross(AC)).norm();
        return 0.5*area;
    }
};
//*********************************************************
///
/// \brief The Plane class
/// d is the double such that, for every point P belonging to the plane, Normal.dot(P)=d
struct Plane {
    Vector3d Normal;
    double d;
    Plane() = default;
    Plane(const Vector3d& newNormal, const double new_d) {
        Normal = newNormal;
        d = new_d;
    }
};
//*********************************************************
///
/// \brief The PolygonalMesh class
///
struct PolygonalMesh {

    unsigned int NumberCell0D = 0;
    vector<unsigned int> Cell0DId = {};
    vector<Vector3d> Cell0DCoordinates = {};

    unsigned int NumberCell1D = 0;
    vector<unsigned int> Cell1DId = {};
    vector<Vector2i> Cell1DVertices = {};

    unsigned int NumberCell2D = 0;
    vector<unsigned int> Cell2DId = {};
    vector<vector<unsigned int>> Cell2DVertices = {};
    vector<vector<unsigned int>> Cell2DEdges = {};

};
//*********************************************************
///
/// \brief TriangulatePolygons: triangulates a variable number of polygons
/// \param Cell2DCoordinates: given l=num_polygons=Cell2DCoordinates.size(), Cell2DCoordinates[l] is a vector<unsigned int> that contains the ids of the 0D cells that form that polygon
/// \return triangulation of every polygon
///
vector<vector<vector<unsigned int>>> TriangulatePolygons(const vector<vector<unsigned int>> &Cell2DCoordinates);
//*********************************************************
///
/// \brief computePolygonsArea: computes area of a variable number of polygons
/// \param listVertices: given l=num_polygons=listVertices.size(), listVertices[l] is a vector<unsigned int> that contains the ids of the 0D cells that form that polygon
/// \param VerticesCoordinates: if i is the id of a certain 0D cell, VerticesCoordinates[i] will be the Vector3d containing its coordinates
/// \return area of every polygon
///
vector<double> computePolygonsArea(const vector<vector<unsigned int>>& listVertices, const vector<Vector3d>& VerticesCoordinates);
//*********************************************************
///
/// \brief FractMeshGedimInterface
/// \param Cell2DVertices: given l=num_polygons=Cell2DVertices.size(), Cell2DVertices[l] is a vector<unsigned int> that contains the ids of the 0D cells that form that polygon
/// \param triangles: will contain, for every triangle in the mesh, the ids of the vertices that form that triangle
/// \param materials
///
void FractMeshGedimInterface(const vector<vector<unsigned int>>& Cell2DVertices, vector<vector<unsigned int>>& triangles, VectorXi& materials);
//*********************************************************
///
/// \brief pointInsidePolygon: checks whether a given points is inside or outside a certain convex polygon
/// !!! FUNCTION DOES NOT SUPPORT CASES WHERE P BELONGS TO ONE OF THE EDGES !!!
/// \param tol: tolerance for comparisons between doubles
/// \param P: point
/// \param VerticesIDs: ids of the vertices of the polygon
/// \param VerticesCoordinates: for every id i, VerticesCoordinates[i] contains the coordinates of the vertex with id i
/// \return true if P is inside the polygon, false otherwise
///
bool pointInsidePolygon(double tol, const Vector3d& P, const vector<unsigned int>& VerticesIDs, const vector<Vector3d>& VerticesCoordinates, Vector3d &PlaneNormal);
//*********************************************************
///
/// \brief computeSquaredDistancePoints: computes the squared distance between two given points in R^3
/// \param Point1
/// \param Point2
/// \return squared distance between Point1 and Point2
///
double computeSquaredDistancePoints(const Vector3d& Point1, const Vector3d& Point2);
//*********************************************************
}
