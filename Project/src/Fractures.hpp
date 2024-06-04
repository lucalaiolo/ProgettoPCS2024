#pragma once
#include<iostream>
#include "Eigen/Eigen"
#include "Polygons.hpp"
using namespace std;
using namespace Eigen;
using namespace GeometryLibrary;
namespace DFNLibrary {
//*********************************************************
///
/// \brief The Traces struct
/// For every trace, given its id i, TraceIDFractures[i] is an array<unsigned int,2> whose elements correspond to the ids
/// of the fractures whose intersection is the trace with id i.
/// TraceCoordinates[i] is a (3x2) matrix of doubles where TraceCoordinates[i].col(0) are the coordinates of the first end point of the trace
/// with id i and TraceCoordinates[i].col(1) are the coordinates of the second end point of the same trace.
/// TraceLength[i] is the length of the trace with id i.
struct Traces {

    vector<array<unsigned int,2>> TraceIDFractures = {};
    vector<MatrixXd> TraceCoordinates = {};
    vector<double> TraceLength = {};

};
//*********************************************************
///
/// \brief The Fractures class
/// For every fracture, given its id i, FractVertices[i] is a (3xn_i) matrx where FractVertices[i].col(j) are the coordinates of the j-th
/// vertex of the fracture.
/// listTraces[i][true] is a vector<unsigned int> that contains the ids of the traces belonging to the fracture with id i such that one
/// of their end points do not belong to the edges of the fracture considered.
///  listTraces[i][false] is a vector<unsigned int> that contains the ids of the traces belonging to the fracture with id i such that both
/// of their end points belong to the edges of the fracture considered.
/// FractPlanes[i] contains, for every fracture with id i, the plane in which that fracture lies
/// FractMeanPoint[i] contains the coordinates of the mean point of the fracture with id i
/// FractMesh[i] is a PolygonalMesh object that corresponds to the polygonal mesh of the fracture with id i computed using the function computePolygonalMesh
///
//*********************************************************
struct Fractures {

    vector<MatrixXd> FractVertices = {};
    vector<map<bool,vector<unsigned int>>> listTraces = {};

    vector<Plane> FractPlanes = {};
    vector<Vector3d> FractMeanPoint = {};
    vector<PolygonalMesh> FractMesh = {};

};
//*********************************************************
///
/// \brief importFractureList: Import the list of fractures
/// @param filepath: the filepath of the input file
/// @param FractureList: a Fractures struct
/// \return the result of the reading, true if success, false otherwise
///
bool importFractureList(const string& filepath, Fractures& FractureList);
//*********************************************************
///
/// \brief computeTraces: calculates the intersections between fractures
/// \param FractureList: a Fractures struct
/// \param TracesList: a Traces struct
/// \param tol: tolerance for comparisons between doubles
///
void computeTraces(Fractures &FractureList, Traces &TracesList, const double &tol);
//*********************************************************
///
/// \brief computeSquaredDistancePoints: computes the squared distance between two given points in R^3
/// \param Point1
/// \param Point2
/// \return squared distance between Point1 and Point2
///
inline double computeSquaredDistancePoints(const Vector3d& Point1, const Vector3d& Point2);
//*********************************************************
///
/// \brief checkIntersectionPossibility: check if there is any possibility that two fractures will intersect.
/// Using the coordinates of the vertices of the fractures, we calculate the mean points of the two fractures.
/// If we define rho_i as the maximum distance between the vertices of the i-th fracture and its mean point,
/// If rho_i+rho_j < dist(mean_point_i,mean_point_j), then i-th fracture and j-th fracture will not intersect.
/// It is a sufficient condition, not necessary.
/// \param Fract1_vertices: (3xn) matrix of doubles. n is the number of vertices of the first fracture considered.
/// Thus, Fract1_vertices.col(i) is a Vector3d whose components are the coordinates of the i-th point of the fracture.
/// \param Fract2_vertices: (3xm) matrix of doubles. m is the number of vertices of the second fracture considered.
/// Thus, Fract2_vertices.col(j) is a Vector3d whose components are the coordinates of the j-th point of the fracture.
/// \param meanPoint1: coordinates of the mean point of the first fracture considered
/// \param meanPoint2: coordinates of the mean point of the second fracture considered
/// \param tol: tolerance for comparisons between doubles
/// \return true if intersection is possible, false otherwise
///
bool checkIntersectionPossibility(const MatrixXd& Fract1_vertices, const MatrixXd& Fract2_vertices,
                                  const Vector3d &meanPoint1, const Vector3d &meanPoint2, const double &tol);
//*********************************************************
///
/// \brief findCase: given two fractures, find their relative position (case)
/// \param betas: let r = {x: x = P + beta*t} the line of intersection between the plane in which the first fracture lies and the
/// plane in which the second fracture lies. Let A_B_C_D be the vector<Vector3d> of size 4 where A_B_C_D[0], A_B_C_D[1] are the the points
/// of intersection between r and the edges of the first fracture and A_B_C_D[2], A_B_C_D[3] are the points of intersection between r and
/// the edges of the second fracture. betas[i] are the values such that P + betas[i]t = A_B_C_D[i].
/// \param tol: tolerance for comparisons between doubles
/// \return integer number representing the case
///
int findCase(const vector<double>& betas, const double &tol);
//*********************************************************
///
/// \brief executeCase: according to the integer returned from the findCase function, calculate the trace
/// \param count: unsigned integer representing the id of the candidate trace
/// \param i: id of first fracture
/// \param j: id of second fracture
/// \param pos1: unsigned integer index such that A_B_C_D[pos1] will be the first end point of the trace
/// \param pos2: unsigned integer index such that A_B_C_D[pos2] will be the second end point of the trace
/// \param A_B_C_D: let r = {x: x = P + beta*t} the line of intersection between the plane in which the i-th fracture lies and the
/// plane in which the j-th fracture lies. A_B_C_D is the vector<Vector3d> of size 4 where A_B_C_D[0], A_B_C_D[1] are the the points
/// of intersection between r and the edges of the i-th fracture and A_B_C_D[2], A_B_C_D[3] are the points of intersection between r and
/// the edges of the j-th fracture.
/// \param tips: tips[0] is true if one of the end points of the trace do not belong to the edges of the i-th fracture, false otherwise.
/// tips[1]  is true if one of the end points of the trace do not belong to the edges of the j-th fracture, false otherwise.
/// \param FractureList: a Fractures struct
/// \param TracesList: a Traces struct
///
inline void executeCase(unsigned int& count, const unsigned int& i, const unsigned int& j, const unsigned int& pos1, const unsigned int& pos2,
                        const vector<Vector3d>& A_B_C_D,const array<bool,2>& tips, Fractures& FractureList, Traces& TracesList);
//*********************************************************
///
/// \brief printTraces: exports the list of traces computed
/// \param outputFileName: name of file in which the list of traces will be printed
/// \param TracesList: a Traces struct
///
void exportTraces(const string& outputFileName, const Traces &TracesList);
//*********************************************************
///
/// \brief printFractures: exports the list of fractures along with information about their traces
/// \param outputFileName: name of file in which the list of fractures will be printed
/// \param FractureList: a Fractures struct
/// \param TracesList: a Traces struct
///
void exportFractures(const string& outputFileName, Fractures& FractureList, const Traces &TracesList);
//*********************************************************
///
/// \brief computeTracesSquaredLength: for each trace, computes its squared length
/// \param TracesList: a Traces struct
///
void computeTracesSquaredLength(Traces& TracesList);
//*********************************************************
///
/// \brief exportParaview: exports a .inp file readable by paraview. It's possible to visualize the fractures using that file
/// \param outputFileName: name of output file
///
void exportParaview(const string& outputFileName, const Fractures& FractureList);
//*********************************************************
///
/// \brief computePolygonalMesh: computes, for every fracture, the polygonal mesh obtained by cutting that fracture with its traces
/// \param FractureList:  a Fractures struct
/// \param TracesList: a Traces struct
/// \param tol: tolerance for comparisons between doubles
///
void computePolygonalMesh(Fractures& FractureList, const Traces& TracesList, const double& tol);
//*********************************************************
///
/// \brief cut: cut one polygon of a polygonal mesh using a certain segment (trace, in our case)
/// \param mesh: a PolygonalMesh struct
/// \param solVec: intersection between a particular polygon and a segment
/// \param edges_ids_sol: ids of the edges that intersect the segment
/// \param tempVec: variable that tells us if the intersection coincides with one of the current vertices or not
/// \param iter: vector<unsigned int> of size 2. when computing the intersection between the edges of the polygon and the segment, we access every edge using a for cycle.
/// iter[0] indicates the iteration at which we found the intersection solVec[0]. same for iter[1]
/// \param l: id of the polygon considered
/// \param tol: tolerance for comparisons between doubles
///
void cut(PolygonalMesh& mesh, vector<Vector3d>& solVec, const vector<unsigned int>& edges_ids_sol, vector<int>& tempVec, vector<unsigned int>& iter, const unsigned int &l, const double &tol);
//*********************************************************
///
/// \brief exportFractureMesh: produces an ouput file readable by paraview. using that, it's possible to visualize the polygonal mesh of a certain fracture obtained using computePolygonalMesh
/// \param outputFileName
/// \param FractureList: a Fractures struct
/// \param fractureID: id of the fracture that is considered
///
void exportFractureMesh(const string& outputFileName, const Fractures& FractureList, const unsigned int& fractureID);
//*********************************************************
///
/// \brief printFractureMesh: exports to 3 .csv files the 0D, 1D and 2D cells of the polygonal mesh computed using computePolygonalMesh
/// \param FractureList: a Fractures struct
/// \param outputFileName1: name of file that will contain the information about the 0D cells
/// \param outputFileName2: name of file that will contain the information about the 1D cells
/// \param outputFileName3: name of file that will contain the information about the 2D cells
/// \return true if the operation was successful, false otherwise
///
bool printFractureMesh(const string &outputFileName1, const string &outputFileName2, const string &outputFileName3, const Fractures &FractureList);
//*********************************************************
}


