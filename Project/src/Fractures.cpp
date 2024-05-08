#include "Fractures.hpp"
#include<iostream>
#include "Eigen/Eigen"
#include <fstream>

using namespace std;
using namespace Eigen;

namespace FractureLibrary {

    void computeTraces(const Fractures &FractureList) {
        const unsigned int numFractures = FractureList.FractVertices.size();
        for(unsigned int i=0;i<numFractures;i++) {
            for(unsigned int j=i+1;j<numFractures;j++){ // we select two fractures without looking for intersections twice

                cout << i<< endl;
                cout <<j<< endl;
                // we must check if it's possible that two fractures can intersect
                const MatrixXd Fract1_vertices = FractureList.FractVertices[i];
                const MatrixXd Fract2_vertices = FractureList.FractVertices[j];

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
                cout << rho1 << endl;
                cout << rho2 << endl;
                const double meanPointsDistance = sqrt(computeSquaredDistancePoints(Fract1_meanPoint,Fract2_meanPoint));
                cout << meanPointsDistance << endl;
                if(rho1+rho2 < meanPointsDistance) {
                    continue;
                }
                // else, compute traces





            }
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

}
