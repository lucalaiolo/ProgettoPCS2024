#include "Fractures.hpp"
#include<iostream>
#include "Eigen/Eigen"
#include <fstream>
#include <iomanip>
using namespace std;
using namespace Eigen;

namespace FractureLibrary {

    void computeTraces(Fractures &FractureList) {
        const unsigned int numFractures = FractureList.FractVertices.size();
        for(unsigned int i=0;i<numFractures;i++) {
            for(unsigned int j=i+1;j<numFractures;j++){ // we select two fractures without looking for intersections twice
                cout << "Fracture " << i<< endl;
                cout << "Fracture " << j<< endl << endl;
                // we must check if it's possible that two fractures can intersect
                const MatrixXd Fract1_vertices = FractureList.FractVertices[i];
                const MatrixXd Fract2_vertices = FractureList.FractVertices[j];

                const unsigned int Fract1_numVertices = Fract1_vertices.cols();
                const unsigned int Fract2_numVertices = Fract2_vertices.cols();

                if(!checkIntersectionPossibility(Fract1_vertices,Fract2_vertices)) {
                    continue;
                }
                // else, compute traces  
                // Step 1: compute plane in which the i-th fracture lies
                const Vector3d P1_P0 =  (Fract1_vertices.col(1)-Fract1_vertices.col(0));
                const Vector3d P2_P0 = (Fract1_vertices.col(2)-Fract1_vertices.col(0));
                const Vector3d n1 = (P1_P0).cross(P2_P0);
                const double d1 = n1.dot(Fract1_vertices.col(0));
                // Step 2: compute plane in which the j-th fracture lies
                const Vector3d Q1_Q0 =  (Fract2_vertices.col(1)-Fract2_vertices.col(0));
                const Vector3d Q2_Q0 = (Fract2_vertices.col(2)-Fract2_vertices.col(0));
                const Vector3d n2 = (Q1_Q0).cross(Q2_Q0);
                const double d2 = n2.dot(Fract2_vertices.col(0));

                // find line of intersection between the two planes
                const Vector3d t = n1.cross(n2);
                if(t.dot(t) < 2.2e-16) {
                    cout << "Planes are parallel."<< endl;
                    continue;
                }
                Matrix3d M1;
                M1 << n1(0), n1(1), n1(2),
                     n2(0), n2(1), n2(2),
                     t(0), t(1), t(2);
                Vector3d b;
                b << d1,d2,0.0;
                Vector3d P = M1.fullPivLu().solve(b);

                // we have now succesfully calculated the line of intersection between the planes in which fracture i and j lie

                // We now find points of intersection between given line and edges of fracture i
                vector<Vector3d> A_B_C_D = {};
                vector<double> betas;

                MatrixXd M2 = MatrixXd::Zero(3,2);
                for(unsigned int k=0;k<Fract1_numVertices;k++) {
                    Vector3d E1 = Fract1_vertices.col(k);
                    Vector3d E2 = Fract1_vertices.col((k+1)%Fract1_numVertices);
                    M2.col(0) = E1-E2;
                    M2.col(1) = -t;
                    Vector3d c = P - E2;
                    FullPivLU<Eigen::MatrixXd> lu_decomp(M2);
                    if(lu_decomp.rank() != 2) {
                        continue;
                    }
                    Vector2d alpha_beta = lu_decomp.solve(c);
                    if(fabs(alpha_beta(0)) > 1) {
                        continue; // intersection does not belong to the segment
                    }
                    betas.push_back(alpha_beta(1));
                    A_B_C_D.push_back(alpha_beta(0)*E1+(1-alpha_beta(0))*E2);
                }
                if(A_B_C_D.size() != 2) {
                    continue; // no intersection
                }
                // repeat for fracture j
                for(unsigned int k=0;k<Fract2_numVertices;k++) {
                    Vector3d E1 = Fract2_vertices.col(k);
                    Vector3d E2 = Fract2_vertices.col((k+1)%Fract2_numVertices);
                    M2.col(0) = E1-E2;
                    M2.col(1) = -t;
                    Vector3d c = P - E2;
                    FullPivLU<Eigen::MatrixXd> lu_decomp(M2);
                    if(lu_decomp.rank() != 2) {
                        continue;
                    }
                    Vector2d alpha_beta = lu_decomp.solve(c);
                    if(fabs(alpha_beta(0)) > 1) {
                        continue; // intersection does not belong to the segment
                    }
                    betas.push_back(alpha_beta(1));
                    A_B_C_D.push_back(alpha_beta(0)*E1+(1-alpha_beta(0))*E2);
                }
                if(A_B_C_D.size() != 4) {
                    cout << "No intersection" << endl; // no intersection
                    continue;
                }

                // it's possible that there is a trace
                MatrixXd tracePoints = MatrixXd::Zero(3,2);
                for(unsigned int k=0;k<4;k++) {
                    cout << fixed << showpoint << setprecision(30) << betas[k] << endl;
                }
                bool boh = (betas[0]==betas[2]);
                bool boh2 = (betas[1]==betas[3]);
                cout << boh << endl;
                cout << boh2 << endl;
                unsigned int n = findCase(A_B_C_D,betas);


                break;
            }
            break;
        }

    }

    void importFractureList(const string &filepath, Fractures &FractureList) {
        ifstream inputFile(filepath);
        if(inputFile.fail()) {
            throw runtime_error("Something went wrong.");
        }
        string line;
        getline(inputFile,line); // skip header
        getline(inputFile,line);
        unsigned int numFractures;
        istringstream convertN(line);
        convertN >> numFractures;
        convertN.clear(); // clear istringstream object
        FractureList.FractVertices.resize(numFractures);
        for(unsigned int i=0;i<numFractures;i++) {
            unsigned int numVertices;
            getline(inputFile,line); // skip line
            getline(inputFile,line,';');
            getline(inputFile,line);
            convertN.str(line);
            convertN >> numVertices;
            convertN.clear(); // clear istringstream object
            FractureList.FractVertices[i] = MatrixXd::Zero(3,numVertices);
            getline(inputFile,line); // skip line
            for(unsigned int j=0; j<3; j++) {
                getline(inputFile,line);
                replace(line.begin(),line.end(),';',' ');
                convertN.str(line);
                for(unsigned int k=0;k<numVertices;k++) {
                    convertN >> FractureList.FractVertices[i](j,k);
                }
                convertN.clear(); // clear istringstream object
            }
        }
        inputFile.close();
    }

    double computeSquaredDistancePoints(const Vector3d Point1, const Vector3d Point2) {
        double dist=0.0;
        for(unsigned int i=0;i<3;i++) {
            dist += (Point1(i)-Point2(i))*(Point1(i)-Point2(i));
        }
        return dist;
    }

    bool checkIntersectionPossibility(const MatrixXd &Fract1_vertices, const MatrixXd &Fract2_vertices) {
        const unsigned int Fract1_numVertices = Fract1_vertices.cols();
        const unsigned int Fract2_numVertices = Fract2_vertices.cols();

        Vector3d Fract1_meanPoint;
        Vector3d Fract2_meanPoint;
        for(unsigned int k=0;k<Fract1_numVertices;k++) {
            Fract1_meanPoint +=  Fract1_vertices.col(k);
        }
        for(unsigned int k=0;k<Fract2_numVertices;k++) {
            Fract2_meanPoint += Fract2_vertices.col(k);
        }
        Fract1_meanPoint = Fract1_meanPoint/Fract1_numVertices;
        Fract2_meanPoint = Fract2_meanPoint/Fract2_numVertices;

        // We find the maximum distance between Fract[i]_meanPoint and the vertices of the fracture i

        double squared_rho1 = computeSquaredDistancePoints(Fract1_meanPoint, Fract1_vertices.col(0));
        double squared_rho2 = computeSquaredDistancePoints(Fract2_meanPoint,Fract2_vertices.col(0));
        double temp1  = 0.0;
        double temp2 = 0.0;
        for(unsigned int k=1;k<Fract1_numVertices;k++) {
            temp1 = computeSquaredDistancePoints(Fract1_meanPoint,Fract1_vertices.col(k));
            if(temp1 > squared_rho1) {
                squared_rho1 = temp1;
            }
        }
        for(unsigned int k=1;k<Fract2_numVertices;k++) {
            temp2 = computeSquaredDistancePoints(Fract2_meanPoint,Fract2_vertices.col(k));
            if(temp2 > squared_rho2) {
                squared_rho2 = temp2;
            }
        }
        const double rho1 = sqrt(squared_rho1);
        const double rho2 = sqrt(squared_rho2);
        const double meanPointsDistance = sqrt(computeSquaredDistancePoints(Fract1_meanPoint,Fract2_meanPoint));
        if(rho1+rho2 < meanPointsDistance) {
            return false;
        }
        return true;
    }

    unsigned int findCase(const vector<Vector3d> &A_B_C_D, const vector<double> &betas) {
        return 0;
    }
}