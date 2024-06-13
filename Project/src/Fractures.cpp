#include "Fractures.hpp"
#include "Sorting.hpp"
#include "Polygons.hpp"
#include "UCDUtilities.hpp"
#include<iostream>
#include "Eigen/Eigen"
#include <fstream>
#include <iomanip>
#include "fstream"
using namespace std;
using namespace Eigen;
using namespace SortingLibrary;
using namespace GeometryLibrary;
namespace DFNLibrary {
//*********************************************************
    bool importFractureList(const string &filepath, Fractures &FractureList) {
        ifstream inputFile(filepath);
        if(inputFile.fail()) {
            cerr << "Something went wrong." << endl;
            return false;
        }
        string line;
        getline(inputFile,line); // skip header
        getline(inputFile,line);
        unsigned int numFractures;
        istringstream convertN(line);
        convertN >> numFractures;
        convertN.clear(); // clear istringstream object
        FractureList.FractVertices.resize(numFractures);
        FractureList.listTraces.resize(numFractures);
        FractureList.FractPlanes.resize(numFractures);
        FractureList.FractMeanPoint.resize(numFractures);
        FractureList.FractMesh.resize(numFractures);

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
            double mean_x = 0.0;
            double mean_y = 0.0;
            double mean_z = 0.0;
            for(unsigned int j=0; j<3; j++) {
                getline(inputFile,line);
                replace(line.begin(),line.end(),';',' ');
                convertN.str(line);
                for(unsigned int k=0;k<numVertices;k++) {
                    convertN >> FractureList.FractVertices[i](j,k);
                    if(j==0){
                        mean_x += FractureList.FractVertices[i](j,k);
                    } else if(j==1){
                        mean_y += FractureList.FractVertices[i](j,k);
                    } else {
                        mean_z += FractureList.FractVertices[i](j,k);
                    }
                }
                convertN.clear(); // clear istringstream object
            }
            // compute mean point of i-th fracture
            Vector3d mean_point = {mean_x, mean_y, mean_z};
            mean_point = mean_point/numVertices;
            FractureList.FractMeanPoint[i] = mean_point;

            // compute plane of i-th fracture
            const Vector3d P1_P0 = FractureList.FractVertices[i].col(1) - FractureList.FractVertices[i].col(0);
            const Vector3d P2_P0 = FractureList.FractVertices[i].col(2) - FractureList.FractVertices[i].col(0);
            Vector3d n = (P1_P0).cross(P2_P0);
            n = n/n.norm();
            const double d = n.dot(FractureList.FractVertices[i].col(0));
            Plane FractPlane(n,d);
            FractureList.FractPlanes[i] = (FractPlane);
        }
        inputFile.close();
        return true;
    }
//*********************************************************
    void computeTraces(Fractures &FractureList, Traces &TracesList, const double& tol) {

        const unsigned int numFractures = FractureList.FractVertices.size();
        unsigned int count = 0; // variable that will tell us the id of the computed trace
        vector<list<unsigned int>> tempTipsTracesList;
        vector<list<unsigned int>> tempNonTipsTracesList;
        list<array<unsigned int,2>> tempTraceIDFractures;
        list<MatrixXd> tempTraceCoordinates;

        tempTipsTracesList.resize(numFractures);
        tempNonTipsTracesList.resize(numFractures);

        for(unsigned int i=0;i<numFractures;i++) {
            for(unsigned int j=i+1;j<numFractures;j++){ // intersection is commutative

                const MatrixXd Fract1_vertices = FractureList.FractVertices[i];
                const MatrixXd Fract2_vertices = FractureList.FractVertices[j];

                const unsigned int Fract1_numVertices = Fract1_vertices.cols();
                const unsigned int Fract2_numVertices = Fract2_vertices.cols();

                if(!checkIntersectionPossibility(Fract1_vertices,Fract2_vertices,FractureList.FractMeanPoint[i],FractureList.FractMeanPoint[j],tol)) {
                    continue;
                }
                // else, compute traces

                // Plane 1:
                const Vector3d n1 = FractureList.FractPlanes[i].Normal;
                const double d1 = FractureList.FractPlanes[i].d;

                // Plane 2:
                const Vector3d n2 = FractureList.FractPlanes[j].Normal;
                const double d2 = FractureList.FractPlanes[j].d;

                // find line of intersection between the two planes
                const Vector3d t = n1.cross(n2);
                if(fabs(t(0)) < tol && fabs(t(1))<tol && fabs(t(2)) < tol) {
                    // cout << "Planes are parallel."<< endl;
                    continue;
                }
                Matrix3d M1;
                M1 << n1(0), n1(1), n1(2),
                     n2(0), n2(1), n2(2),
                     t(0), t(1), t(2);
                Vector3d b;
                b << d1,d2,0.0;
                Vector3d P = M1.fullPivLu().solve(b);

                // we have now succesfully calculated the line of intersection between the planes in which the two fractures lie

                // We now find points of intersection between given line and edges of fracture i
                vector<Vector3d> A_B_C_D = {};
                vector<double> betas;
                A_B_C_D.reserve(4);
                betas.reserve(4);

                MatrixXd M2 = MatrixXd::Zero(3,2);
                for(unsigned int k=0;k<Fract1_numVertices;k++) {
                    const Vector3d E1 = Fract1_vertices.col(k);
                    const Vector3d E2 = Fract1_vertices.col((k+1)%Fract1_numVertices);
                    M2.col(0) = E1-E2;
                    M2.col(1) = -t;
                    const Vector3d c = P - E2;

                    // we check whether the segment and the line are coplanar and non parallel
                    const Vector3d temp_vec = t.cross((E1-E2)/((E1-E2).norm()));
                    if(fabs(temp_vec(0))<tol && fabs(temp_vec(1))<tol && fabs(temp_vec(2))<tol) {
                        continue; // segment and line are parallel
                    }

                    if(fabs((c/(c.norm())).dot(temp_vec)) > tol) {
                        continue; // distance between lines is positive --> no intersection
                    }
                    // else, we must compute the intersection

                    Vector2d alpha_beta = M2.fullPivLu().solve(c);


                    if(alpha_beta(0) > 1.0+tol || alpha_beta(0) < -tol) {
                        continue; // intersection does not belong to the segment
                    }

                    betas.push_back(alpha_beta(1));
                    A_B_C_D.push_back(alpha_beta(0)*E1+(1-alpha_beta(0))*E2);
                    if(A_B_C_D.size()==2) { // optimization: given our hypothesis of convexity, there can't be more than two points of intersection
                        break;
                    }
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

                    // we check whether the segment and the line are coplanar and non parallel
                    const Vector3d temp_vec = t.cross((E1-E2)/((E1-E2).norm()));
                    if(fabs(temp_vec(0))<tol && fabs(temp_vec(1))<tol && fabs(temp_vec(2))<tol) {
                        continue; // segment and line are parallel
                    }

                    if(fabs((c/(c.norm())).dot(temp_vec)) > tol) {
                        continue; // distance between lines is positive --> no intersection
                    }

                    // else, we must compute the intersection

                    Vector2d alpha_beta = M2.fullPivLu().solve(c);

                    if(alpha_beta(0) > 1.0+tol|| alpha_beta(0) < -tol) {
                        continue; // intersection does not belong to the segment
                    }
                    betas.push_back(alpha_beta(1));
                    A_B_C_D.push_back(alpha_beta(0)*E1+(1-alpha_beta(0))*E2);
                    if(A_B_C_D.size()==4) { // same as before
                        break;
                    }
                }
                if(A_B_C_D.size() != 4) {
                    //cout << "No intersection" << endl; // no intersection
                    continue;
                }

                // it's possible that there is a trace

                double temp;
                Vector3d temp_vec;

                // "sort" the points of intersection on the line for the first fracture
                if(betas[0] > betas[1]) {
                    temp = betas[0];
                    betas[0] = betas[1];
                    betas[1] = temp;
                    temp_vec = A_B_C_D[0];
                    A_B_C_D[0] = A_B_C_D[1];
                    A_B_C_D[1] = temp_vec;
                }

                // "sort" the points of intersection on the line for the first fracture
                if(betas[2] > betas[3]) {
                    temp = betas[2];
                    betas[2] = betas[3];
                    betas[3] = temp;
                    temp_vec = A_B_C_D[2];
                    A_B_C_D[2] = A_B_C_D[3];
                    A_B_C_D[3] = temp_vec;
                }


                const int n = findCase(betas,tol);

                if(n==0) {
                    //cout << "Trace is A_B. False for both." << endl;
                    array<bool,2> tips = {false,false};
                    executeCase(count,i,j,0,1,A_B_C_D,tips, tempTipsTracesList, tempNonTipsTracesList, tempTraceIDFractures, tempTraceCoordinates);
                    continue;
                } else if(n == 1 || n == 5) {
                    //cout << "No trace" << endl;
                    continue;
                } else if(n==-1 || n==2) {
                    //cout << "Trace is A_B. False F1. True F2." << endl;
                    array<bool,2> tips = {false,true};
                    executeCase(count,i,j,0,1,A_B_C_D,tips,tempTipsTracesList, tempNonTipsTracesList, tempTraceIDFractures, tempTraceCoordinates);
                    continue;
                } else if(n==3) {
                    //cout << Trace is CB. True both." << endl;
                    array<bool,2> tips = {true,true};
                    executeCase(count,i,j,2,1,A_B_C_D,tips,tempTipsTracesList, tempNonTipsTracesList, tempTraceIDFractures, tempTraceCoordinates);
                    continue;
                } else if(n==4) {
                    //cout << Trace is AD. True both." << endl;
                    array<bool,2> tips = {true,true};
                    executeCase(count,i,j,0,3,A_B_C_D,tips,tempTipsTracesList, tempNonTipsTracesList, tempTraceIDFractures, tempTraceCoordinates);
                    continue;
                } else if(n==6 || n==-2 || n==-3) {
                    //cout << "Trace is CD. F1 true. F2 false." << endl;
                    array<bool,2> tips = {true,false};
                    executeCase(count,i,j,2,3,A_B_C_D,tips,tempTipsTracesList, tempNonTipsTracesList, tempTraceIDFractures, tempTraceCoordinates);
                    continue;
                } else if(n==-4) {
                    //cout << "Trace is AD. false F1. true F2" << endl;
                    array<bool,2> tips = {false,true};
                    executeCase(count,i,j,0,3,A_B_C_D,tips,tempTipsTracesList, tempNonTipsTracesList, tempTraceIDFractures, tempTraceCoordinates);
                    continue;
                }

            }
        }
        for(unsigned int i=0;i<numFractures;i++) {
            const unsigned int numTipsTraces = tempTipsTracesList[i].size();
            const unsigned int numNonTipsTraces = tempNonTipsTracesList[i].size();
            FractureList.listTraces[i][true].reserve(numTipsTraces);
            FractureList.listTraces[i][false].reserve(numNonTipsTraces);
            for(auto& traceid : tempTipsTracesList[i]) {
                FractureList.listTraces[i][false].push_back(traceid);
            }
            for(auto& traceid : tempNonTipsTracesList[i]) {
                FractureList.listTraces[i][true].push_back(traceid);
            }
        }

        const unsigned int numTraces = tempTraceIDFractures.size();
        TracesList.TraceCoordinates.reserve(numTraces);
        TracesList.TraceIDFractures.reserve(numTraces);
        for(auto& ids:tempTraceIDFractures) {
            TracesList.TraceIDFractures.push_back(ids);
        }
        for(auto& coord:tempTraceCoordinates) {
            TracesList.TraceCoordinates.push_back(coord);
        }

        computeTracesSquaredLength(TracesList);
        exportTraces("Traces.txt", TracesList);
        exportFractures("Fractures.txt", FractureList, TracesList);

    }
//*********************************************************
    inline void executeCase(unsigned int& count, const unsigned int& i, const unsigned int& j, const unsigned int& pos1, const unsigned int& pos2,
                            const vector<Vector3d>& A_B_C_D,const array<bool,2>& tips, vector<list<unsigned int>> &tempTipsTracesList,
                            vector<list<unsigned int>> &tempNonTipsTracesList, list<array<unsigned int,2>> &tempTraceIDFractures,
                            list<MatrixXd> &tempTraceCoordinates) {
        MatrixXd trace_coord = MatrixXd::Zero(3,2);
        trace_coord.col(0) = A_B_C_D[pos1];
        trace_coord.col(1) = A_B_C_D[pos2];

        if(tips[0]==true) {
            tempNonTipsTracesList[i].push_back(count);
        } else {
            tempTipsTracesList[i].push_back(count);
        }
        if(tips[1]==true) {
            tempNonTipsTracesList[j].push_back(count);
        } else {
            tempTipsTracesList[j].push_back(count);
        }
        array<unsigned int,2> trace_id = {i,j};
        tempTraceIDFractures.push_back(trace_id);
        tempTraceCoordinates.push_back(trace_coord);
        count++;
    }
//*********************************************************
    bool checkIntersectionPossibility(const MatrixXd &Fract1_vertices, const MatrixXd &Fract2_vertices,
                                      const Vector3d& meanPoint1, const Vector3d& meanPoint2, const double& tol) {
        const unsigned int Fract1_numVertices = Fract1_vertices.cols();
        const unsigned int Fract2_numVertices = Fract2_vertices.cols();

        // We find the maximum distance between Fract[i]_meanPoint and the vertices of the fracture i

        double squared_rho1 = computeSquaredDistancePoints(meanPoint1, Fract1_vertices.col(0));
        double squared_rho2 = computeSquaredDistancePoints(meanPoint2,Fract2_vertices.col(0));
        double temp1  = 0.0;
        double temp2 = 0.0;
        for(unsigned int k=1;k<Fract1_numVertices;k++) {
            temp1 = computeSquaredDistancePoints(meanPoint1,Fract1_vertices.col(k));
            if(temp1 > squared_rho1) {
                squared_rho1 = temp1;
            }
        }
        for(unsigned int k=1;k<Fract2_numVertices;k++) {
            temp2 = computeSquaredDistancePoints(meanPoint2,Fract2_vertices.col(k));
            if(temp2 > squared_rho2) {
                squared_rho2 = temp2;
            }
        }
        const double rho1 = sqrt(squared_rho1);
        const double rho2 = sqrt(squared_rho2);
        const double meanPointsDistance = sqrt(computeSquaredDistancePoints(meanPoint1,meanPoint2));
        if(rho1+rho2 < meanPointsDistance-tol) {
            return false;
        }
        return true;
    }
//*********************************************************
    int findCase(const vector<double> &betas, const double &tol) {
        // we first check if any of these points coincide

        const bool A_equals_C = (fabs(betas[0]-betas[2])/max(max(fabs(betas[0]),fabs(betas[2])),{1}) < tol);
        bool A_equals_D = (fabs(betas[0]-betas[3])/max(max(fabs(betas[0]),fabs(betas[3])),{1}) < tol);
        if(A_equals_D) {
            return 1;
        }
        bool B_equals_C = (fabs(betas[1]-betas[2])/max(max(fabs(betas[1]),fabs(betas[2])),{1}) < tol);
        if(B_equals_C) {
            return 1;
        }
        const bool B_equals_D = (fabs(betas[1]-betas[3])/max(max(fabs(betas[1]),fabs(betas[3])),{1}) < tol);
        if(A_equals_C && B_equals_D) {
            //cout << "Trace is A_B. Tips for both." << endl;
            return 0;
        } else {
            if(A_equals_C) { // && !B_equals_D
                if(betas[1] < betas[3]) {
                    //cout << "Case 6.1 " << endl;
                    //cout << "Trace is AB.passante F1.no passante F2." << endl;
                    return -1;
                } else {
                    //cout << "Case 6.2" << endl;
                    //cout << "Trace is CD. Non passante for F1. passante for F2." << endl;
                    return -2;
                }
            } else if(B_equals_D) { // && !A_equals_C
                if(betas[0] < betas[2]) {
                    //cout << "Case 7.1" << endl;
                    //cout << "Trace is CD. Non passante  F1. passante F2" << endl;
                    return -3;
                } else {
                    //cout << "Case 7.2" << endl;
                    //cout << "Trace is AB. passante F1. non passante F2" << endl;
                    return -4;
                }
            } else { // no corresponding points
                if(betas[1]< betas[2]) {
                    //cout << "Case 1. No trace." << endl;
                    return 1;
                }
                if(betas[2] < betas[0] && betas[1] < betas[3]) {
                    //cout << "Case 2. Trace is AB. Passante F1. Non passante F2." << endl;
                    return 2;
                }
                if(betas[0] < betas[2] && betas[2] < betas[1] && betas[1] < betas[3]) {
                    //cout << "Case 3. Trace is CB. Non passante." << endl;
                    return 3;
                }
                if(betas[2] < betas[0] && betas[0] < betas[3] && betas[3] < betas[1]) {
                    //cout << "Case 4. Trace is AD. Non passante." << endl;
                    return 4;
                }
                if(betas[3] < betas[0]) {
                    //cout << "Case 5. No trace." << endl;
                    return 5;
                }
                if(betas[0] < betas[2] && betas[3] < betas[1]) {
                    //cout << "Case 6. Trace is CD. F1 non passante. F2 passante." << endl;
                    return 6;
                }

            }

        }
        cerr << "Something went wrong" << endl;
        return 100;
    }
//*********************************************************
    void exportTraces(const string& outputFileName, const Traces &TracesList) {
        ofstream file;
        file.open(outputFileName);
        if(file.fail()) {
            cerr << "Something went wrong opening output file." << endl;
        }
        file << "# Number of Traces" << endl;
        const unsigned int num_traces = TracesList.TraceIDFractures.size();
        file << num_traces << endl;
        file << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
        for(unsigned int k=0;k<num_traces;k++) {
            file << k << "; ";
            file << scientific << setprecision(16) << TracesList.TraceIDFractures[k][0] << "; ";
            file << TracesList.TraceIDFractures[k][1] << "; ";
            file << TracesList.TraceCoordinates[k](0,0) << "; ";
            file << TracesList.TraceCoordinates[k](1,0) << "; ";
            file << TracesList.TraceCoordinates[k](2,0) << "; ";
            file << TracesList.TraceCoordinates[k](0,1) << "; ";
            file << TracesList.TraceCoordinates[k](1,1) << "; ";
            file << TracesList.TraceCoordinates[k](2,1) << endl;
        }
    }
//*********************************************************
    void exportFractures(const string &outputFileName, Fractures &FractureList, const Traces &TracesList) {
        ofstream file;
        file.open(outputFileName);
        if(file.fail()) {
            cerr << "Something went wrong opening output file." << endl;
        }
        const unsigned int num_fractures = FractureList.FractVertices.size();
        for(unsigned int k=0; k<num_fractures; k++) {
            file << "#FractureId; NumTraces" << endl;
            const unsigned int num_non_tips_traces = FractureList.listTraces[k][false].size();
            const unsigned int num_tips_traces = FractureList.listTraces[k][true].size();
            file << k << "; " << num_non_tips_traces + num_tips_traces << endl;
            if(num_non_tips_traces + num_tips_traces == 0) {
                continue;
            }
            file << "# TraceId; Tips; Length" << endl;

            if(num_non_tips_traces > 0) {
                MergeSortTraces(FractureList.listTraces[k][false], TracesList.TraceLength); // sort traces by their length
                for(unsigned int l=0;l<FractureList.listTraces[k][false].size();l++) {
                    file << FractureList.listTraces[k][false][l] << "; " << false << "; " << scientific << setprecision(16) << TracesList.TraceLength[FractureList.listTraces[k][false][l]] << ";" << endl;
                }
            }

            if(num_tips_traces > 0) {
                MergeSortTraces(FractureList.listTraces[k][true], TracesList.TraceLength); // sort traces by their length
                for(unsigned int l=0;l<FractureList.listTraces[k][true].size();l++) {
                    file << FractureList.listTraces[k][true][l] << "; " << true << "; " << scientific << setprecision(16) << TracesList.TraceLength[FractureList.listTraces[k][true][l]] << ";" << endl;
                }
            }
        }

    }
//*********************************************************
    void computeTracesSquaredLength(Traces &TracesList) {
        const unsigned int num_traces = TracesList.TraceCoordinates.size();
        TracesList.TraceLength.reserve(num_traces);
        for(unsigned int k=0;k<num_traces;k++) {
            const Vector3d P1 = TracesList.TraceCoordinates[k].col(0);
            const Vector3d P2 = TracesList.TraceCoordinates[k].col(1); // points of first trace
            TracesList.TraceLength.push_back(sqrt(computeSquaredDistancePoints(P1,P2))); // squared length of trace
        }
    }
//*********************************************************
    void computePolygonalMesh(Fractures &FractureList, const Traces &TracesList, const double &tol) {
        const unsigned int numFractures = FractureList.FractVertices.size();

         for(unsigned int i=0; i<numFractures; i++) {
            // initialize polygonal mesh
            cout << "Frattura " << i << endl;
            MatrixXd vertices = FractureList.FractVertices[i];
            unsigned int num_vertices = vertices.cols();
            unsigned int num_edges = num_vertices;
            unsigned int num_polygons = 1; // first iteration, the only polygon is the fracture itself

            FractureList.FractMesh[i].NumberCell0D = num_vertices;
            FractureList.FractMesh[i].Cell0DCoordinates.resize(num_vertices);
            FractureList.FractMesh[i].Cell0DId.resize(num_vertices);

            FractureList.FractMesh[i].NumberCell1D = num_edges;
            FractureList.FractMesh[i].Cell1DId.resize(num_edges);
            FractureList.FractMesh[i].Cell1DVertices.resize(num_edges);

            FractureList.FractMesh[i].NumberCell2D = 1;
            FractureList.FractMesh[i].Cell2DId.resize(1);
            FractureList.FractMesh[i].Cell2DVertices.resize(1);
            FractureList.FractMesh[i].Cell2DEdges.resize(1);

            for(unsigned int j=0; j<num_vertices;j++) {
                FractureList.FractMesh[i].Cell0DId[j] = j;
                FractureList.FractMesh[i].Cell0DCoordinates[j] = vertices.col(j);

                FractureList.FractMesh[i].Cell1DId[j] = j;
                FractureList.FractMesh[i].Cell1DVertices[j] = {j, (j+1)%num_vertices};
            }

            FractureList.FractMesh[i].Cell2DId[0] = 0;
            FractureList.FractMesh[i].Cell2DVertices[0] = FractureList.FractMesh[i].Cell0DId;
            FractureList.FractMesh[i].Cell2DEdges[0] = FractureList.FractMesh[i].Cell1DId;

            // first polygonal mesh initialized
            // we will now update it

            // calculate the polygonal mesh for i-th fracture
            vector<unsigned int> passing_traces = FractureList.listTraces[i][false];
            vector<unsigned int> non_passing_traces = FractureList.listTraces[i][true];

            const unsigned int num_passing_traces = passing_traces.size();
            const unsigned int num_non_passing_traces = non_passing_traces.size();

            if(num_passing_traces > 0) {
                // compute
                FractureList.FractMesh[i].Cell0DId.reserve(2*num_passing_traces);
                FractureList.FractMesh[i].Cell0DCoordinates.reserve(2*num_passing_traces);

                FractureList.FractMesh[i].Cell1DId.reserve(5*num_passing_traces);
                FractureList.FractMesh[i].Cell1DVertices.reserve(5*num_passing_traces);

                FractureList.FractMesh[i].Cell2DId.reserve(num_passing_traces);
                FractureList.FractMesh[i].Cell2DVertices.reserve(num_passing_traces);
                FractureList.FractMesh[i].Cell2DEdges.reserve(num_passing_traces);

                for(unsigned int j=0; j<num_passing_traces; j++) {
                    // iterate on passing traces
                    //cout << j << endl;
                    num_polygons = FractureList.FractMesh[i].NumberCell2D;
                    const unsigned int trace_id = passing_traces[j];
                    const MatrixXd trace_end_points = TracesList.TraceCoordinates[trace_id];
                    const Vector3d Q0 = trace_end_points.col(0);
                    const Vector3d Q1 = trace_end_points.col(1);
                    for(unsigned int l=0;l<num_polygons;l++) { // iterate on every 2D cell

                        // CHECK WHETHER THE TRACE CUTS THE POLYGON, AND IN WHICH POINTS
                        num_edges = FractureList.FractMesh[i].Cell2DEdges[l].size();

                        vector<int> tempVec;
                        tempVec.reserve(2);
                        // at the end of the procedure, tempVec(0) will equal -1 if the first intersection does not coincide with any of the
                        // vertices of the 2D cell. otherwise, it will be equal to the id of the edge in which that intersection is
                        // the starting point. same for tempVec(1).


                        vector<Vector3d> solVec = {};
                        vector<unsigned int> edges_ids_sol = {};
                        solVec.reserve(2); // solVec will contain the coordinates of the 2 points of intersection between the trace
                        // and the edges of the 2D cell, if they exist
                        edges_ids_sol.reserve(2); // edges_ids_sol will contain the ids of the edges in which
                        // the possible two intersections lie
                        bool skip = false;
                        vector<unsigned int> iter;
                        iter.reserve(2);
                        for(unsigned int k=0;k<num_edges;k++) { // iterate on every edge of that 2D cell

                            if(skip) {
                                skip = false;
                                continue;
                            }

                            const Vector3d P0 = FractureList.FractMesh[i].Cell0DCoordinates[FractureList.FractMesh[i].Cell2DVertices[l][k]];
                            const Vector3d P1 = FractureList.FractMesh[i].Cell0DCoordinates[FractureList.FractMesh[i].Cell2DVertices[l][(k+1)%num_edges]];

                            MatrixXd M = MatrixXd::Zero(3,2);
                            M.col(0) = Q0-Q1;
                            M.col(1) = P1-P0;
                            const Vector3d b = P1-Q1;

                            const Vector3d temp = (P1-P0).cross(Q1-Q0);  // check whether segments are parallel or not

                            /// GESTIRE TOLLERANZA
                            if(fabs(temp(0)) < tol && fabs(temp(1)) < tol && fabs(temp(2)) < tol) {
                                continue; // there is no intersection, segments are parallel
                            }

                            FullPivLU<Eigen::MatrixXd> lu_decomp(M);
                            const Vector2d alpha_beta = lu_decomp.solve(b);

                            const double alpha = alpha_beta(0);
                            const double beta = alpha_beta(1);
                            const Vector3d P = alpha*(Q0-Q1)+Q1; // Point of intersection between the two lines
                            /// GESTIRE TOLLERANZE

                            if((alpha < -tol) || (alpha > 1+tol) || (beta < -tol) || (beta >1+tol)) {
                                continue;
                            } else {
                                if(-tol <= beta && beta <= tol) {
                                    // beta is basically equal to zero
                                    // if this happens, the intersection coincides with P1
                                    skip = true;
                                    solVec.push_back(P1);
                                    edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][(k+1)%num_edges]); // because P1 is the starting point of the next edge
                                    tempVec.push_back(FractureList.FractMesh[i].Cell2DEdges[l][(k+1)%num_edges]); // because P1 is the starting point of the next edge
                                    iter.push_back((k+1)%num_edges);
                                } else if(tol < beta && beta < 1-tol) {
                                    // intersection coincides with P
                                    solVec.push_back(P);
                                    edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]);
                                    tempVec.push_back(-1);
                                    iter.push_back(k);
                                } else if(1-tol <= beta && beta <= 1+ tol) {
                                    // beta is basically equal to 1
                                    // if this happens, the intersection coincides with P0
                                    solVec.push_back(P0);
                                    edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]); // because P0 is the starting point of this edge
                                    tempVec.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]); // because P0 is the starting point of this edge
                                    iter.push_back(k);
                                }
                            }

                            if(solVec.size()==2) {
                                break;
                            }

                        } // chiudo il ciclo sui lati

                        // we have now successfully calculated the intersection between the passing trace and the 2D
                        // cell with Id l

                        // we must now update the polygonal mesh

                        if(solVec.size()!=2) {
                            continue; // no update
                        }
                        if((solVec[0]-solVec[1]).norm() < tol) {
                            continue;
                        }

                        cut(FractureList.FractMesh[i], solVec, edges_ids_sol, tempVec, iter,l,tol,false);

                    } // chiudo ciclo sui poligoni

                } // chiudo ciclo sulle tracce

            } // chiudo if sul numero delle tracce passanti

            FractureList.FractMesh[i].Cell0DId.shrink_to_fit();
            FractureList.FractMesh[i].Cell0DCoordinates.shrink_to_fit();

            FractureList.FractMesh[i].Cell1DId.shrink_to_fit();
            FractureList.FractMesh[i].Cell1DVertices.shrink_to_fit();

            FractureList.FractMesh[i].Cell2DId.shrink_to_fit();
            FractureList.FractMesh[i].Cell2DVertices.shrink_to_fit();
            FractureList.FractMesh[i].Cell2DEdges.shrink_to_fit();



            if(num_non_passing_traces > 0) {
                // compute
                FractureList.FractMesh[i].Cell0DId.reserve(2*num_non_passing_traces);
                FractureList.FractMesh[i].Cell0DCoordinates.reserve(2*num_non_passing_traces);

                FractureList.FractMesh[i].Cell1DId.reserve(5*num_non_passing_traces);
                FractureList.FractMesh[i].Cell1DVertices.reserve(5*num_non_passing_traces);

                FractureList.FractMesh[i].Cell2DId.reserve(num_non_passing_traces);
                FractureList.FractMesh[i].Cell2DVertices.reserve(num_non_passing_traces);
                FractureList.FractMesh[i].Cell2DEdges.reserve(num_non_passing_traces);

                for(unsigned int j=0; j<num_non_passing_traces; j++) {
                    // iterate on passing traces
                    //cout << j << endl;
                    num_polygons = FractureList.FractMesh[i].NumberCell2D;
                    const unsigned int trace_id = non_passing_traces[j];
                    const MatrixXd trace_end_points = TracesList.TraceCoordinates[trace_id];
                    const Vector3d Q0 = trace_end_points.col(0);
                    const Vector3d Q1 = trace_end_points.col(1);
                    for(unsigned int l=0;l<num_polygons;l++) { // iterate on every 2D cell
                        // check whether Q0 and Q1 belong to the edges of the polygon
                        // check whether Q0 and Q1 belong to the edges of the polygon
                        int id1 = -1;
                        int id2 = -1;
                        const unsigned int num_edges = FractureList.FractMesh[i].Cell2DEdges[l].size();
                        for(unsigned int k=0;k<num_edges;k++) {
                            const Vector3d P0 = FractureList.FractMesh[i].Cell0DCoordinates[FractureList.FractMesh[i].Cell2DVertices[l][k]];
                            const Vector3d P1 = FractureList.FractMesh[i].Cell0DCoordinates[FractureList.FractMesh[i].Cell2DVertices[l][(k+1)%num_edges]];
                            const double edge_length = sqrt(computeSquaredDistancePoints(P0,P1));
                            if(id1==-1) {
                                const double lengthQ0P0 = sqrt(computeSquaredDistancePoints(Q0,P0));
                                const double lengthQ0P1 = sqrt(computeSquaredDistancePoints(Q0,P1));
                                if(fabs(edge_length-lengthQ0P0-lengthQ0P1)<=tol) {
                                    id1 = k;
                                }
                            }
                            if(id2==-1) {
                                const double lengthQ1P0 = sqrt(computeSquaredDistancePoints(Q1,P0));
                                const double lengthQ1P1 = sqrt(computeSquaredDistancePoints(Q1,P1));
                                if(fabs(edge_length-lengthQ1P0-lengthQ1P1)<=tol) {
                                    id2 = k;
                                }
                            }
                            if(id1!=-1 && id2!=-1) {
                                break;
                            }
                        }

                        // check whether trace lies inside the polygon or not
                        bool inside_cellQ0 = false;
                        bool inside_cellQ1 = false;

                        if(id1==-1) {
                            inside_cellQ0 = pointInsidePolygon(tol,Q0,FractureList.FractMesh[i].Cell2DVertices[l],FractureList.FractMesh[i].Cell0DCoordinates,FractureList.FractPlanes[i].Normal);
                        }
                        if(id2==-1) {
                            inside_cellQ1 = pointInsidePolygon(tol,Q1,FractureList.FractMesh[i].Cell2DVertices[l],FractureList.FractMesh[i].Cell0DCoordinates,FractureList.FractPlanes[i].Normal);
                        }
                        if(id1==-1 && id2==-1) {
                            if(inside_cellQ0 || inside_cellQ1) {

                                vector<int> tempVec;
                                tempVec.reserve(2);
                                // at the end of the procedure, temp(0) will equal -1 if the first intersection does not coincide with any of the
                                // vertices of the 2D cell. otherwise, it will be equal to the id of the edge in which that intersection is
                                // the starting point. same for temp(1).

                                vector<Vector3d> solVec = {};
                                vector<unsigned int> edges_ids_sol = {};
                                solVec.reserve(2); // solVec will contain the coordinates of the 2 points of intersection between the trace
                                // and the edges of the 2D cell, if they exist
                                edges_ids_sol.reserve(2); // edges_ids_sol will contain the ids of the edges in which
                                // the possible two intersections lie
                                bool skip = false;
                                vector<unsigned int> iter;
                                iter.reserve(2);
                                for(unsigned int k=0;k<num_edges;k++) { // iterate on every edge of that 2D cell

                                    if(skip) {
                                        skip = false;
                                        continue;
                                    }

                                    const Vector3d P0 = FractureList.FractMesh[i].Cell0DCoordinates[FractureList.FractMesh[i].Cell2DVertices[l][k]];
                                    const Vector3d P1 = FractureList.FractMesh[i].Cell0DCoordinates[FractureList.FractMesh[i].Cell2DVertices[l][(k+1)%num_edges]];

                                    MatrixXd M = MatrixXd::Zero(3,2);
                                    M.col(0) = Q0-Q1;
                                    M.col(1) = P1-P0;
                                    const Vector3d b = P1-Q1;

                                    const Vector3d temp = (P1-P0).cross(Q1-Q0);  // check whether segments are parallel or not

                                    /// GESTIRE TOLLERANZA
                                    if(fabs(temp(0)) < tol && fabs(temp(1)) < tol && fabs(temp(2)) < tol) {
                                        continue; // there is no intersection, segments are parallel
                                    }

                                    FullPivLU<Eigen::MatrixXd> lu_decomp(M);
                                    const Vector2d alpha_beta = lu_decomp.solve(b);

                                    const double alpha = alpha_beta(0);
                                    const double beta = alpha_beta(1);
                                    const Vector3d P = alpha*(Q0-Q1)+Q1; // Point of intersection between the two lines
                                    /// GESTIRE TOLLERANZE

                                    if( (beta < -tol) || (beta >1+tol)) {
                                        continue;
                                    } else {
                                        if(-tol <= beta && beta <= tol) {
                                            // beta is basically equal to zero
                                            // if this happens, the intersection coincides with P1
                                            skip = true;
                                            solVec.push_back(P1);
                                            edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][(k+1)%num_edges]); // because P1 is the starting point of the next edge
                                            tempVec.push_back(FractureList.FractMesh[i].Cell2DEdges[l][(k+1)%num_edges]); // because P1 is the starting point of the next edge
                                            iter.push_back((k+1)%num_edges);
                                        } else if(tol < beta && beta < 1-tol) {
                                            // intersection coincides with P
                                            solVec.push_back(P);
                                            edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]);
                                            tempVec.push_back(-1);
                                            iter.push_back(k);
                                        } else if(1-tol <= beta && beta <= 1+ tol) {
                                            // beta is basically equal to 1
                                            // if this happens, the intersection coincides with P0
                                            solVec.push_back(P0);
                                            edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]); // because P0 is the starting point of this edge
                                            tempVec.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]); // because P0 is the starting point of this edge
                                            iter.push_back(k);
                                        }
                                    }

                                    if(solVec.size()==2) {
                                        break;
                                    }

                                } // chiudo il ciclo sui lati

                                // we have now successfully calculated the intersection between the passing trace and the 2D
                                // cell with Id l

                                // we must now update the polygonal mesh

                                if(solVec.size()!=2) {
                                    continue; // no update
                                }
                                if((solVec[0]-solVec[1]).norm() < tol) {
                                    continue;
                                }

                                cut(FractureList.FractMesh[i],solVec,edges_ids_sol,tempVec,iter,l,tol,true);
                                if(inside_cellQ0 && inside_cellQ1) {
                                    break;
                                }
                            } else if(!inside_cellQ0 && !inside_cellQ1) {
                                vector<int> tempVec;
                                tempVec.reserve(2);
                                // at the end of the procedure, temp(0) will equal -1 if the first intersection does not coincide with any of the
                                // vertices of the 2D cell. otherwise, it will be equal to the id of the edge in which that intersection is
                                // the starting point. same for temp(1).

                                vector<Vector3d> solVec = {};
                                vector<unsigned int> edges_ids_sol = {};
                                solVec.reserve(2); // solVec will contain the coordinates of the 2 points of intersection between the trace
                                // and the edges of the 2D cell, if they exist
                                edges_ids_sol.reserve(2); // edges_ids_sol will contain the ids of the edges in which
                                // the possible two intersections lie
                                bool skip = false;
                                vector<unsigned int> iter;
                                iter.reserve(2);
                                for(unsigned int k=0;k<num_edges;k++) { // iterate on every edge of that 2D cell

                                    if(skip) {
                                        skip = false;
                                        continue;
                                    }

                                    const Vector3d P0 = FractureList.FractMesh[i].Cell0DCoordinates[FractureList.FractMesh[i].Cell2DVertices[l][k]];
                                    const Vector3d P1 = FractureList.FractMesh[i].Cell0DCoordinates[FractureList.FractMesh[i].Cell2DVertices[l][(k+1)%num_edges]];

                                    MatrixXd M = MatrixXd::Zero(3,2);
                                    M.col(0) = Q0-Q1;
                                    M.col(1) = P1-P0;
                                    const Vector3d b = P1-Q1;

                                    const Vector3d temp = (P1-P0).cross(Q1-Q0);  // check whether segments are parallel or not

                                    /// GESTIRE TOLLERANZA
                                    if(fabs(temp(0)) < tol && fabs(temp(1)) < tol && fabs(temp(2)) < tol) {
                                        continue; // there is no intersection, segments are parallel
                                    }

                                    FullPivLU<Eigen::MatrixXd> lu_decomp(M);
                                    const Vector2d alpha_beta = lu_decomp.solve(b);

                                    const double alpha = alpha_beta(0);
                                    const double beta = alpha_beta(1);
                                    const Vector3d P = alpha*(Q0-Q1)+Q1; // Point of intersection between the two lines
                                    /// GESTIRE TOLLERANZE

                                    if((alpha < -tol) || (alpha > 1+tol) || (beta < -tol) || (beta >1+tol)) {
                                        continue;
                                    } else {
                                        if(-tol <= beta && beta <= tol) {
                                            // beta is basically equal to zero
                                            // if this happens, the intersection coincides with P1
                                            skip = true;
                                            solVec.push_back(P1);
                                            edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][(k+1)%num_edges]); // because P1 is the starting point of the next edge
                                            tempVec.push_back(FractureList.FractMesh[i].Cell2DEdges[l][(k+1)%num_edges]); // because P1 is the starting point of the next edge
                                            iter.push_back((k+1)%num_edges);
                                        } else if(tol < beta && beta < 1-tol) {
                                            // intersection coincides with P
                                            solVec.push_back(P);
                                            edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]);
                                            tempVec.push_back(-1);
                                            iter.push_back(k);
                                        } else if(1-tol <= beta && beta <= 1+ tol) {
                                            // beta is basically equal to 1
                                            // if this happens, the intersection coincides with P0
                                            solVec.push_back(P0);
                                            edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]); // because P0 is the starting point of this edge
                                            tempVec.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]); // because P0 is the starting point of this edge
                                            iter.push_back(k);
                                        }
                                    }

                                    if(solVec.size()==2) {
                                        break;
                                    }

                                } // chiudo il ciclo sui lati

                                // we have now successfully calculated the intersection between the passing trace and the 2D
                                // cell with Id l

                                // we must now update the polygonal mesh

                                if(solVec.size()!=2) {
                                    continue; // no update
                                }
                                if((solVec[0]-solVec[1]).norm() < tol) {
                                    continue;
                                }

                                cut(FractureList.FractMesh[i],solVec,edges_ids_sol,tempVec,iter,l,tol,true);
                            }

                        } else if((id1!=-1 && id2==-1) || (id1==-1 && id2!=-1)) {
                            // one point is outside or inside the polygon, the other one belongs to one of the edges
                            vector<int> tempVec;
                            tempVec.reserve(2);
                            // at the end of the procedure, temp(0) will equal -1 if the first intersection does not coincide with any of the
                            // vertices of the 2D cell. otherwise, it will be equal to the id of the edge in which that intersection is
                            // the starting point. same for temp(1).

                            vector<Vector3d> solVec = {};
                            vector<unsigned int> edges_ids_sol = {};
                            solVec.reserve(2); // solVec will contain the coordinates of the 2 points of intersection between the trace
                            // and the edges of the 2D cell, if they exist
                            edges_ids_sol.reserve(2); // edges_ids_sol will contain the ids of the edges in which
                            // the possible two intersections lie
                            bool skip = false;
                            vector<unsigned int> iter;
                            iter.reserve(2);
                            for(unsigned int k=0;k<num_edges;k++) { // iterate on every edge of that 2D cell

                                if(skip) {
                                    skip = false;
                                    continue;
                                }

                                const Vector3d P0 = FractureList.FractMesh[i].Cell0DCoordinates[FractureList.FractMesh[i].Cell2DVertices[l][k]];
                                const Vector3d P1 = FractureList.FractMesh[i].Cell0DCoordinates[FractureList.FractMesh[i].Cell2DVertices[l][(k+1)%num_edges]];

                                MatrixXd M = MatrixXd::Zero(3,2);
                                M.col(0) = Q0-Q1;
                                M.col(1) = P1-P0;
                                const Vector3d b = P1-Q1;

                                const Vector3d temp = (P1-P0).cross(Q1-Q0);  // check whether segments are parallel or not

                                /// GESTIRE TOLLERANZA
                                if(fabs(temp(0)) < tol && fabs(temp(1)) < tol && fabs(temp(2)) < tol) {
                                    continue; // there is no intersection, segments are parallel
                                }

                                FullPivLU<Eigen::MatrixXd> lu_decomp(M);
                                const Vector2d alpha_beta = lu_decomp.solve(b);

                                const double alpha = alpha_beta(0);
                                const double beta = alpha_beta(1);
                                const Vector3d P = alpha*(Q0-Q1)+Q1; // Point of intersection between the two lines
                                /// GESTIRE TOLLERANZE

                                if((id1!=-1 && !inside_cellQ1) || (!inside_cellQ0 && id2!=-1)) {
                                    if((alpha < -tol) || (alpha > 1+tol) || (beta < -tol) || (beta >1+tol)) {
                                        continue;
                                    } else {
                                        if(-tol <= beta && beta <= tol) {
                                            // beta is basically equal to zero
                                            // if this happens, the intersection coincides with P1
                                            skip = true;
                                            solVec.push_back(P1);
                                            edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][(k+1)%num_edges]); // because P1 is the starting point of the next edge
                                            tempVec.push_back(FractureList.FractMesh[i].Cell2DEdges[l][(k+1)%num_edges]); // because P1 is the starting point of the next edge
                                            iter.push_back((k+1)%num_edges);
                                        } else if(tol < beta && beta < 1-tol) {
                                            // intersection coincides with P
                                            solVec.push_back(P);
                                            edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]);
                                            tempVec.push_back(-1);
                                            iter.push_back(k);
                                        } else if(1-tol <= beta && beta <= 1+ tol) {
                                            // beta is basically equal to 1
                                            // if this happens, the intersection coincides with P0
                                            solVec.push_back(P0);
                                            edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]); // because P0 is the starting point of this edge
                                            tempVec.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]); // because P0 is the starting point of this edge
                                            iter.push_back(k);
                                        }
                                    }
                                } else if((id1!=-1 && inside_cellQ1) || (inside_cellQ0 && id2!=-1)) {
                                    // alpha can be anything in this case
                                    if((beta < -tol) || (beta >1+tol)) {
                                        continue;
                                    } else {
                                        if(-tol <= beta && beta <= tol) {
                                            // beta is basically equal to zero
                                            // if this happens, the intersection coincides with P1
                                            skip = true;
                                            solVec.push_back(P1);
                                            edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][(k+1)%num_edges]); // because P1 is the starting point of the next edge
                                            tempVec.push_back(FractureList.FractMesh[i].Cell2DEdges[l][(k+1)%num_edges]); // because P1 is the starting point of the next edge
                                            iter.push_back((k+1)%num_edges);
                                        } else if(tol < beta && beta < 1-tol) {
                                            // intersection coincides with P
                                            solVec.push_back(P);
                                            edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]);
                                            tempVec.push_back(-1);
                                            iter.push_back(k);
                                        } else if(1-tol <= beta && beta <= 1+ tol) {
                                            // beta is basically equal to 1
                                            // if this happens, the intersection coincides with P0
                                            solVec.push_back(P0);
                                            edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]); // because P0 is the starting point of this edge
                                            tempVec.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]); // because P0 is the starting point of this edge
                                            iter.push_back(k);
                                        }
                                    }
                                }

                                if(solVec.size()==2) {
                                    break;
                                }

                            } // chiudo il ciclo sui lati

                            // we have now successfully calculated the intersection between the passing trace and the 2D
                            // cell with Id l

                            // we must now update the polygonal mesh

                            if(solVec.size()!=2) {
                                continue; // no update
                            }
                            if((solVec[0]-solVec[1]).norm() < tol) {
                                continue;
                            }

                            cut(FractureList.FractMesh[i],solVec,edges_ids_sol,tempVec,iter,l,tol,true);

                            if((id1!=-1 && inside_cellQ1) || (inside_cellQ0 && id2!=-1)) {
                                break;
                            }
                        } else if(id1 != -1 && id2 != -1) {
                            if(id1 == id2) {
                                break;
                            } else {
                                vector<int> tempVec;
                                tempVec.reserve(2);
                                // at the end of the procedure, temp(0) will equal -1 if the first intersection does not coincide with any of the
                                // vertices of the 2D cell. otherwise, it will be equal to the id of the edge in which that intersection is
                                // the starting point. same for temp(1).

                                vector<Vector3d> solVec = {};
                                vector<unsigned int> edges_ids_sol = {};
                                solVec.reserve(2); // solVec will contain the coordinates of the 2 points of intersection between the trace
                                // and the edges of the 2D cell, if they exist
                                edges_ids_sol.reserve(2); // edges_ids_sol will contain the ids of the edges in which
                                // the possible two intersections lie
                                bool skip = false;
                                vector<unsigned int> iter;
                                iter.reserve(2);
                                for(unsigned int k=0;k<num_edges;k++) { // iterate on every edge of that 2D cell

                                    if(skip) {
                                        skip = false;
                                        continue;
                                    }

                                    const Vector3d P0 = FractureList.FractMesh[i].Cell0DCoordinates[FractureList.FractMesh[i].Cell2DVertices[l][k]];
                                    const Vector3d P1 = FractureList.FractMesh[i].Cell0DCoordinates[FractureList.FractMesh[i].Cell2DVertices[l][(k+1)%num_edges]];

                                    MatrixXd M = MatrixXd::Zero(3,2);
                                    M.col(0) = Q0-Q1;
                                    M.col(1) = P1-P0;
                                    const Vector3d b = P1-Q1;

                                    const Vector3d temp = (P1-P0).cross(Q1-Q0);  // check whether segments are parallel or not

                                    /// GESTIRE TOLLERANZA
                                    if(fabs(temp(0)) < tol && fabs(temp(1)) < tol && fabs(temp(2)) < tol) {
                                        continue; // there is no intersection, segments are parallel
                                    }

                                    FullPivLU<Eigen::MatrixXd> lu_decomp(M);
                                    const Vector2d alpha_beta = lu_decomp.solve(b);

                                    const double alpha = alpha_beta(0);
                                    const double beta = alpha_beta(1);
                                    const Vector3d P = alpha*(Q0-Q1)+Q1; // Point of intersection between the two lines
                                    /// GESTIRE TOLLERANZE

                                    if((alpha < -tol) || (alpha > 1+tol) || (beta < -tol) || (beta >1+tol)) {
                                        continue;
                                    } else {
                                        if(-tol <= beta && beta <= tol) {
                                            // beta is basically equal to zero
                                            // if this happens, the intersection coincides with P1
                                            skip = true;
                                            solVec.push_back(P1);
                                            edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][(k+1)%num_edges]); // because P1 is the starting point of the next edge
                                            tempVec.push_back(FractureList.FractMesh[i].Cell2DEdges[l][(k+1)%num_edges]); // because P1 is the starting point of the next edge
                                            iter.push_back((k+1)%num_edges);
                                        } else if(tol < beta && beta < 1-tol) {
                                            // intersection coincides with P
                                            solVec.push_back(P);
                                            edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]);
                                            tempVec.push_back(-1);
                                            iter.push_back(k);
                                        } else if(1-tol <= beta && beta <= 1+ tol) {
                                            // beta is basically equal to 1
                                            // if this happens, the intersection coincides with P0
                                            solVec.push_back(P0);
                                            edges_ids_sol.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]); // because P0 is the starting point of this edge
                                            tempVec.push_back(FractureList.FractMesh[i].Cell2DEdges[l][k]); // because P0 is the starting point of this edge
                                            iter.push_back(k);
                                        }
                                    }

                                    if(solVec.size()==2) {
                                        break;
                                    }

                                } // chiudo il ciclo sui lati

                                // we have now successfully calculated the intersection between the passing trace and the 2D
                                // cell with Id l

                                // we must now update the polygonal mesh

                                if(solVec.size()!=2) {
                                    continue; // no update
                                }
                                if((solVec[0]-solVec[1]).norm() < tol) {
                                    continue;
                                }

                                cut(FractureList.FractMesh[i],solVec,edges_ids_sol,tempVec,iter,l,tol,true);
                            }

                        }


                    }


                }

            }
            //FractureList.FractMesh[i] = mesh;
        }
    }
//*********************************************************
    void cut(PolygonalMesh &mesh, vector<Vector3d> &solVec, const vector<unsigned int> &edges_ids_sol, vector<int> &tempVec, vector<unsigned int> &iter, const unsigned int &l, const double &tol, const bool& tips) {

        if(tempVec[0]==-1 && tempVec[1]==-1) {
            const unsigned int id_edge_first_intersection = edges_ids_sol[0];
            const Vector3d first_intersection = solVec[0];

            const unsigned int id_edge_second_intersection = edges_ids_sol[1];
            const Vector3d second_intersection = solVec[1];

            vector<Vector3d> newCell0DCoordinates = mesh.Cell0DCoordinates;
            newCell0DCoordinates.reserve(2);
            newCell0DCoordinates.push_back(first_intersection);
            newCell0DCoordinates.push_back(second_intersection);


            vector<unsigned int> IdNewVertices1 = {};
            vector<unsigned int> IdNewVertices2 = {};
            const unsigned int current_number_of_vertices = mesh.Cell2DVertices[l].size();
            IdNewVertices1.reserve(current_number_of_vertices+1);
            IdNewVertices2.reserve(current_number_of_vertices+1); // allocate memory for maximum number of vertices

            vector<unsigned int> IdOldEdges = mesh.Cell2DEdges[l]; // current edges of 2D Cell with id l
            vector<unsigned int> IdNewEdges1 = {};
            vector<unsigned int> IdNewEdges2 = {};
            IdNewEdges1.reserve(IdOldEdges.size()+1);
            IdNewEdges2.reserve(IdOldEdges.size()+1);

            unsigned int k1 = 0;
            unsigned int m1 = 0;

            for(unsigned int k=0;k<IdOldEdges.size();k++) {
                if(IdOldEdges[k]==id_edge_first_intersection) {
                    IdNewVertices1.push_back(mesh.Cell2DVertices[l][k]);
                    k1 = k;
                    break;
                } else {
                    IdNewEdges1.push_back(IdOldEdges[k]);
                    IdNewVertices1.push_back(mesh.Cell2DVertices[l][k]);
                }
            }

            IdNewEdges1.push_back(id_edge_first_intersection);
            IdNewEdges1.push_back(mesh.NumberCell1D+1);
            IdNewVertices1.push_back(mesh.NumberCell0D);
            IdNewVertices1.push_back(mesh.NumberCell0D+1);

            IdNewEdges2.push_back(mesh.NumberCell1D);
            IdNewVertices2.push_back(mesh.NumberCell0D);

            for(unsigned int m=k1+1;m<IdOldEdges.size();m++) {
                if(IdOldEdges[m]==id_edge_second_intersection) {
                    IdNewVertices2.push_back(mesh.Cell2DVertices[l][m]);
                    m1 = m;
                    break;
                } else {
                    IdNewVertices2.push_back(mesh.Cell2DVertices[l][m]);
                    IdNewEdges2.push_back(IdOldEdges[m]);
                }
            }
            IdNewVertices2.push_back(mesh.NumberCell0D+1);
            IdNewEdges2.push_back(id_edge_second_intersection);
            IdNewEdges2.push_back(mesh.NumberCell1D+1);
            IdNewEdges1.push_back(mesh.NumberCell1D+2);
            for(unsigned int n=m1+1;n<IdOldEdges.size();n++) {
                IdNewVertices1.push_back(mesh.Cell2DVertices[l][n]);
                IdNewEdges1.push_back(IdOldEdges[n]);;
            }

            IdNewVertices1.shrink_to_fit();
            IdNewVertices2.shrink_to_fit();
            IdNewEdges1.shrink_to_fit();
            IdNewEdges2.shrink_to_fit();

            /// BEFORE UPDATING THE MESH, CHECK WHETHER THE AREA IS LARGE ENOUGH

            const double area1 = computePolygonsArea({IdNewVertices1},newCell0DCoordinates)[0];
            const double area2 = computePolygonsArea({IdNewVertices2},newCell0DCoordinates)[0];

            if(fabs(area1) <= tol || fabs(area2) <= tol) {
                return; // NO UPDATES, NEW 2D CELL WOULD HAVE AREA APPROXIMATELY EQUAL TO ZERO
            }
            // else, do the updates

            // remove the edges with id id_edge_first_intersection and id_edge_second_intersection
            bool first_int_found = false;
            bool second_int_found = false;
            unsigned int start;
            if(!tips) {
                start = l+1;
            } else {
                start = 0;
            }
            for(unsigned int l1 = start; l1<mesh.NumberCell2D;l1++) {
                // we check every 2D cell
                if(l1==l){
                    continue;
                }
                const vector<unsigned int> edges = mesh.Cell2DEdges[l1];
                for(unsigned int k=0;k<edges.size();k++) {
                    if(!first_int_found) {
                        if(edges[k] == id_edge_first_intersection) {
                            // modify the edges of that cell
                            first_int_found = true;
                            const unsigned int idP0 = mesh.Cell2DVertices[l1][k]; // initial point of that edge
                            const unsigned int idP1 = mesh.Cell2DVertices[l1][(k+1)%edges.size()]; // end point of that edge
                            if(mesh.Cell2DVertices[l][k1]==idP0) {
                                mesh.Cell2DEdges[l1][k] = id_edge_first_intersection;
                                mesh.Cell2DEdges[l1].insert(mesh.Cell2DEdges[l1].begin()+k+1,mesh.NumberCell1D);
                            } else if(mesh.Cell2DVertices[l][k1]==idP1) {
                                mesh.Cell2DEdges[l1][k] = mesh.NumberCell1D;
                                mesh.Cell2DEdges[l1].insert(mesh.Cell2DEdges[l1].begin()+k+1, id_edge_first_intersection);
                            }

                            mesh.Cell2DVertices[l1].insert(mesh.Cell2DVertices[l1].begin()+k+1, mesh.NumberCell0D); // added new vertex to that cell
                            break;
                        }
                    }

                    if(!second_int_found) {
                        if(edges[k] == id_edge_second_intersection) {
                            // modify the edges of that cell
                            second_int_found = true;
                            const unsigned int idP0 = mesh.Cell2DVertices[l1][k]; // initial point of that edge
                            const unsigned int idP1 = mesh.Cell2DVertices[l1][(k+1)%edges.size()]; // end point of that edge
                            if(mesh.Cell2DVertices[l][m1]==idP0) {
                                mesh.Cell2DEdges[l1][k] = id_edge_second_intersection;
                                mesh.Cell2DEdges[l1].insert(mesh.Cell2DEdges[l1].begin()+k+1,mesh.NumberCell1D+2);
                            } else if(mesh.Cell2DVertices[l][m1]==idP1) {
                                mesh.Cell2DEdges[l1][k] = mesh.NumberCell1D+2;
                                mesh.Cell2DEdges[l1].insert(mesh.Cell2DEdges[l1].begin()+k+1, id_edge_second_intersection);
                            }
                            mesh.Cell2DVertices[l1].insert(mesh.Cell2DVertices[l1].begin()+k+1, mesh.NumberCell0D+1); // added new vertex to that cell
                            break;
                        }
                    }

                }
                if(first_int_found && second_int_found) {
                    break;
                }
            }

            //update Cell0Ds
            mesh.Cell0DId.push_back(mesh.NumberCell0D);
            mesh.Cell0DId.push_back(mesh.NumberCell0D+1);
            mesh.Cell0DCoordinates.push_back(first_intersection);
            mesh.Cell0DCoordinates.push_back(second_intersection);

            // update Cell1Ds
            for(unsigned int k=0;k<3;k++) {
                mesh.Cell1DId.push_back(mesh.NumberCell1D+k);
            }
            mesh.Cell1DVertices[id_edge_first_intersection] = {mesh.Cell0DId[mesh.Cell2DVertices[l][iter[0]]],mesh.NumberCell0D};
            mesh.Cell1DVertices.push_back({mesh.NumberCell0D, mesh.Cell0DId[mesh.Cell2DVertices[l][(iter[0]+1)%mesh.Cell2DVertices[l].size()]]});
            mesh.Cell1DVertices.push_back({mesh.NumberCell0D,mesh.NumberCell0D+1});
            mesh.Cell1DVertices[id_edge_second_intersection] = {mesh.Cell0DId[mesh.Cell2DVertices[l][iter[1]]],mesh.NumberCell0D+1};
            mesh.Cell1DVertices.push_back({mesh.NumberCell0D+1, mesh.Cell0DId[mesh.Cell2DVertices[l][(iter[1]+1)%mesh.Cell2DVertices[l].size()]]});

            // update Cell2Ds
            mesh.Cell2DId.push_back(mesh.NumberCell2D);

            mesh.Cell2DVertices[l] = IdNewVertices1;
            mesh.Cell2DVertices.push_back(IdNewVertices2);

            mesh.Cell2DEdges[l] = IdNewEdges1;
            mesh.Cell2DEdges.push_back(IdNewEdges2);

            mesh.NumberCell0D += 2;
            mesh.NumberCell1D += 3;
            mesh.NumberCell2D += 1;

        } else if(tempVec[0]!=-1 && tempVec[1]!=-1){


            const unsigned int id_edge_first_intersection = edges_ids_sol[0];

            const unsigned int id_edge_second_intersection = edges_ids_sol[1];


            // update edges and vertices

            vector<unsigned int> IdNewVertices1 = {};
            vector<unsigned int> IdNewVertices2 = {};
            const unsigned int current_number_of_vertices = mesh.Cell2DVertices[l].size();
            IdNewVertices1.reserve(current_number_of_vertices-1); // allocate memory for maximum number of vertices for the new 2D cell
            IdNewVertices2.reserve(current_number_of_vertices-1);

            vector<unsigned int> IdOldEdges = mesh.Cell2DEdges[l]; // current edges of 2D Cell with id l
            vector<unsigned int> IdNewEdges1 = {};
            vector<unsigned int> IdNewEdges2 = {};
            IdNewEdges1.reserve(IdOldEdges.size()-1); // allocate memory for maximum number of vertices for the new 2D cell
            IdNewEdges2.reserve(IdOldEdges.size()-1);

            unsigned int k1 = 0;
            unsigned int m1 = 0;

            for(unsigned int k=0;k<IdOldEdges.size();k++) {
                if(IdOldEdges[k]==id_edge_first_intersection || IdOldEdges[k] == id_edge_second_intersection) {
                    IdNewVertices1.push_back(mesh.Cell2DVertices[l][k]);
                    IdNewVertices2.push_back(mesh.Cell2DVertices[l][k]);
                    k1 = k;
                    break;
                } else {
                    IdNewEdges1.push_back(IdOldEdges[k]);
                    IdNewVertices1.push_back(mesh.Cell2DVertices[l][k]);
                }
            }

            IdNewEdges1.push_back(mesh.NumberCell1D);

            IdNewEdges2.push_back(IdOldEdges[k1]);
            for(unsigned int m=k1+1;m<IdOldEdges.size();m++) {
                if(IdOldEdges[m]==id_edge_second_intersection) {
                    IdNewVertices1.push_back(mesh.Cell2DVertices[l][m]);
                    IdNewVertices2.push_back(mesh.Cell2DVertices[l][m]);
                    m1 = m;
                    break;
                } else {
                    IdNewEdges2.push_back(IdOldEdges[m]);
                    IdNewVertices2.push_back(mesh.Cell2DVertices[l][m]);
                }
            }
            IdNewEdges2.push_back(mesh.NumberCell1D);
            IdNewEdges1.push_back(IdOldEdges[m1]);
            for(unsigned int n=m1+1;n<IdOldEdges.size();n++) {
                IdNewVertices1.push_back(mesh.Cell2DVertices[l][n]);
                IdNewEdges1.push_back(IdOldEdges[n]);;
            }

            IdNewVertices1.shrink_to_fit();
            IdNewVertices2.shrink_to_fit();
            IdNewEdges1.shrink_to_fit();
            IdNewEdges2.shrink_to_fit();

            /// BEFORE UPDATING THE MESH, CHECK WHETHER THE AREA IS LARGE ENOUGH

            const double area1 = computePolygonsArea({IdNewVertices1},mesh.Cell0DCoordinates)[0];
            const double area2 = computePolygonsArea({IdNewVertices2},mesh.Cell0DCoordinates)[0];

            if(fabs(area1) <= tol || fabs(area2) <= tol) {
                return; // NO UPDATES, NEW 2D CELL WOULD HAVE AREA APPROXIMATELY EQUAL TO ZERO
            }
            // else, do the updates

            // no updates on Cell0Ds

            // update Cell1Ds
            mesh.Cell1DId.push_back(mesh.NumberCell1D);
            mesh.Cell1DVertices.push_back({mesh.Cell0DId[mesh.Cell2DVertices[l][iter[0]]],mesh.Cell0DId[mesh.Cell2DVertices[l][iter[1]]]});

            // update Cell2Ds
            mesh.Cell2DId.push_back(mesh.NumberCell2D);

            mesh.Cell2DVertices[l] = IdNewVertices1;
            mesh.Cell2DVertices.push_back(IdNewVertices2);

            mesh.Cell2DEdges[l] = IdNewEdges1;
            mesh.Cell2DEdges.push_back(IdNewEdges2);

            mesh.NumberCell1D += 1;
            mesh.NumberCell2D += 1;

        } else if((tempVec[0]!=-1 && tempVec[1] == -1)) {

            const unsigned int id_edge_first_intersection = edges_ids_sol[0];

            const unsigned int id_edge_second_intersection = edges_ids_sol[1];
            const Vector3d second_intersection = solVec[1];

            vector<Vector3d> newCell0DCoordinates = mesh.Cell0DCoordinates;
            newCell0DCoordinates.reserve(1);
            newCell0DCoordinates.push_back(second_intersection);

            vector<unsigned int> IdNewVertices1 = {};
            vector<unsigned int> IdNewVertices2 = {};
            const unsigned int current_number_of_vertices = mesh.Cell2DVertices[l].size();
            IdNewVertices1.reserve(current_number_of_vertices);
            IdNewVertices2.reserve(current_number_of_vertices); // allocate memory for maximum number of vertices

            vector<unsigned int> IdOldEdges = mesh.Cell2DEdges[l]; // current edges of 2D Cell with id l
            vector<unsigned int> IdNewEdges1 = {};
            vector<unsigned int> IdNewEdges2 = {};
            IdNewEdges1.reserve(IdOldEdges.size());
            IdNewEdges2.reserve(IdOldEdges.size());

            unsigned int k1 = 0;
            unsigned int m1 = 0;

            for(unsigned int k=0;k<IdOldEdges.size();k++) {
                if(IdOldEdges[k]==id_edge_first_intersection) {
                    IdNewVertices1.push_back(mesh.Cell2DVertices[l][k]);
                    k1 = k;
                    break;
                } else {
                    IdNewEdges1.push_back(IdOldEdges[k]);
                    IdNewVertices1.push_back(mesh.Cell2DVertices[l][k]);
                }
            }

            IdNewEdges1.push_back(mesh.NumberCell1D+1);
            IdNewVertices1.push_back(mesh.NumberCell0D);

            IdNewEdges2.push_back(IdOldEdges[k1]);
            IdNewVertices2.push_back(mesh.Cell2DVertices[l][k1]);

            for(unsigned int m=k1+1;m<IdOldEdges.size();m++) {
                if(IdOldEdges[m]==id_edge_second_intersection) {
                    IdNewVertices2.push_back(mesh.Cell2DVertices[l][m]);
                    m1 = m;
                    break;
                } else {
                    IdNewVertices2.push_back(mesh.Cell2DVertices[l][m]);
                    IdNewEdges2.push_back(IdOldEdges[m]);
                }
            }
            IdNewVertices2.push_back(mesh.NumberCell0D);
            IdNewEdges2.push_back(id_edge_second_intersection);
            IdNewEdges2.push_back(mesh.NumberCell1D+1);
            IdNewEdges1.push_back(mesh.NumberCell1D);
            for(unsigned int n=m1+1;n<IdOldEdges.size();n++) {
                IdNewVertices1.push_back(mesh.Cell2DVertices[l][n]);
                IdNewEdges1.push_back(IdOldEdges[n]);;
            }

            IdNewVertices1.shrink_to_fit();
            IdNewVertices2.shrink_to_fit();
            IdNewEdges1.shrink_to_fit();
            IdNewEdges2.shrink_to_fit();

            /// BEFORE UPDATING THE MESH, CHECK WHETHER THE AREA IS LARGE ENOUGH

            const double area1 = computePolygonsArea({IdNewVertices1},newCell0DCoordinates)[0];
            const double area2 = computePolygonsArea({IdNewVertices2},newCell0DCoordinates)[0];

            if(fabs(area1) <= tol || fabs(area2) <= tol) {
                return; // NO UPDATES, NEW 2D CELL WOULD HAVE AREA APPROXIMATELY EQUAL TO ZERO
            }
            // else, do the updates

            // remove the edges with id id_edge_first_intersection and id_edge_second_intersection
            bool second_int_found = false;
            unsigned int start;
            if(!tips) {
                start = l+1;
            } else {
                start = 0;
            }
            for(unsigned int l1 = start; l1<mesh.NumberCell2D;l1++) {
                if(l1==l) {
                    continue;
                }
                const vector<unsigned int> edges = mesh.Cell2DEdges[l1];
                for(unsigned int k=0;k<edges.size();k++) {
                    if(!second_int_found) {
                        if(edges[k] == id_edge_second_intersection) {
                            // modify the edges of that cell
                            second_int_found = true;
                            const unsigned int idP0 = mesh.Cell2DVertices[l1][k]; // initial point of that edge
                            const unsigned int idP1 = mesh.Cell2DVertices[l1][(k+1)%edges.size()]; // end point of that edge
                            if(mesh.Cell2DVertices[l][m1]==idP0) {
                                mesh.Cell2DEdges[l1][k] = id_edge_second_intersection;
                                mesh.Cell2DEdges[l1].insert(mesh.Cell2DEdges[l1].begin()+k+1,mesh.NumberCell1D);
                            } else if(mesh.Cell2DVertices[l][m1]==idP1) {
                                mesh.Cell2DEdges[l1][k] = mesh.NumberCell1D;
                                mesh.Cell2DEdges[l1].insert(mesh.Cell2DEdges[l1].begin()+k+1, id_edge_second_intersection);
                            }
                            mesh.Cell2DVertices[l1].insert(mesh.Cell2DVertices[l1].begin()+k+1, mesh.NumberCell0D); // added new vertex to that cell
                            break;
                        }
                    }

                }
                if(second_int_found) {
                    break;
                }
            }

            //update Cell0Ds
            mesh.Cell0DId.push_back(mesh.NumberCell0D);
            mesh.Cell0DCoordinates.push_back(second_intersection);

            // update Cell1Ds
            for(unsigned int k=0;k<2;k++) {
                mesh.Cell1DId.push_back(mesh.NumberCell1D+k);
            }
            mesh.Cell1DVertices[id_edge_second_intersection] = {mesh.Cell0DId[mesh.Cell2DVertices[l][iter[1]]],mesh.NumberCell0D};
            mesh.Cell1DVertices.push_back({mesh.NumberCell0D, mesh.Cell0DId[mesh.Cell2DVertices[l][(iter[1]+1)%mesh.Cell2DVertices[l].size()]]});
            mesh.Cell1DVertices.push_back({mesh.NumberCell0D,mesh.Cell0DId[mesh.Cell2DVertices[l][iter[0]]]});

            // update Cell2Ds
            mesh.Cell2DId.push_back(mesh.NumberCell2D);

            mesh.Cell2DVertices[l] = IdNewVertices1;
            mesh.Cell2DVertices.push_back(IdNewVertices2);

            mesh.Cell2DEdges[l] = IdNewEdges1;
            mesh.Cell2DEdges.push_back(IdNewEdges2);

            mesh.NumberCell0D += 1;
            mesh.NumberCell1D += 2;
            mesh.NumberCell2D += 1;

        } else if(tempVec[0]==-1 && tempVec[1]!=-1) {

            const unsigned int id_edge_first_intersection = edges_ids_sol[0];

            const unsigned int id_edge_second_intersection = edges_ids_sol[1];
            const Vector3d first_intersection = solVec[0];

            vector<Vector3d> newCell0DCoordinates = mesh.Cell0DCoordinates;
            newCell0DCoordinates.reserve(1);
            newCell0DCoordinates.push_back(first_intersection);

            vector<unsigned int> IdNewVertices1 = {};
            vector<unsigned int> IdNewVertices2 = {};
            const unsigned int current_number_of_vertices = mesh.Cell2DVertices[l].size();
            IdNewVertices1.reserve(current_number_of_vertices);
            IdNewVertices2.reserve(current_number_of_vertices); // allocate memory for maximum number of vertices

            vector<unsigned int> IdOldEdges = mesh.Cell2DEdges[l]; // current edges of 2D Cell with id l
            vector<unsigned int> IdNewEdges1 = {};
            vector<unsigned int> IdNewEdges2 = {};
            IdNewEdges1.reserve(IdOldEdges.size());
            IdNewEdges2.reserve(IdOldEdges.size());

            unsigned int k1 = 0;
            unsigned int m1 = 0;

            for(unsigned int k=0;k<IdOldEdges.size();k++) {
                if(IdOldEdges[k]==id_edge_first_intersection) {
                    IdNewVertices1.push_back(mesh.Cell2DVertices[l][k]);
                    k1 = k;
                    break;
                } else {
                    IdNewEdges1.push_back(IdOldEdges[k]);
                    IdNewVertices1.push_back(mesh.Cell2DVertices[l][k]);
                }
            }

            IdNewEdges1.push_back(id_edge_first_intersection);
            IdNewEdges1.push_back(mesh.NumberCell1D+1);
            IdNewVertices1.push_back(mesh.NumberCell0D);

            IdNewVertices2.push_back(mesh.NumberCell0D);
            IdNewEdges2.push_back(mesh.NumberCell1D);

            for(unsigned int m=k1+1;m<IdOldEdges.size();m++) {
                if(IdOldEdges[m]==id_edge_second_intersection) {
                    IdNewVertices2.push_back(mesh.Cell2DVertices[l][m]);
                    m1 = m;
                    break;
                } else {
                    IdNewVertices2.push_back(mesh.Cell2DVertices[l][m]);
                    IdNewEdges2.push_back(IdOldEdges[m]);
                }
            }
            IdNewEdges2.push_back(mesh.NumberCell1D+1);
            IdNewEdges1.push_back(IdOldEdges[m1]);
            IdNewVertices1.push_back(mesh.Cell2DVertices[l][m1]);
            for(unsigned int n=m1+1;n<IdOldEdges.size();n++) {
                IdNewVertices1.push_back(mesh.Cell2DVertices[l][n]);
                IdNewEdges1.push_back(IdOldEdges[n]);;
            }

            IdNewVertices1.shrink_to_fit();
            IdNewVertices2.shrink_to_fit();
            IdNewEdges1.shrink_to_fit();
            IdNewEdges2.shrink_to_fit();

            /// BEFORE UPDATING THE MESH, CHECK WHETHER THE AREA IS LARGE ENOUGH

            const double area1 = computePolygonsArea({IdNewVertices1},newCell0DCoordinates)[0];
            const double area2 = computePolygonsArea({IdNewVertices2},newCell0DCoordinates)[0];

            if(fabs(area1) <= tol || fabs(area2) <= tol) {
                return; // NO UPDATES, NEW 2D CELL WOULD HAVE AREA APPROXIMATELY EQUAL TO ZERO
            }
            // else, do the updates

            // remove the edges with id id_edge_first_intersection and id_edge_second_intersection
            bool first_int_found = false;
            unsigned int start;
            if(!tips) {
                start = l+1;
            } else {
                start = 0;
            }
            for(unsigned int l1 = start; l1<mesh.NumberCell2D;l1++) {
                // we check every 2D cell
                if(l1==l) {
                    continue;
                }
                const vector<unsigned int> edges = mesh.Cell2DEdges[l1];
                for(unsigned int k=0;k<edges.size();k++) {
                    if(!first_int_found) {
                        if(edges[k] == id_edge_first_intersection) {
                            // modify the edges of that cell
                            first_int_found = true;
                            const unsigned int idP0 = mesh.Cell2DVertices[l1][k]; // initial point of that edge
                            const unsigned int idP1 = mesh.Cell2DVertices[l1][(k+1)%edges.size()]; // end point of that edge
                            if(mesh.Cell2DVertices[l][k1]==idP0) {
                                mesh.Cell2DEdges[l1][k] = id_edge_first_intersection;
                                mesh.Cell2DEdges[l1].insert(mesh.Cell2DEdges[l1].begin()+k+1,mesh.NumberCell1D);
                            } else if(mesh.Cell2DVertices[l][k1]==idP1) {
                                mesh.Cell2DEdges[l1][k] = mesh.NumberCell1D;
                                mesh.Cell2DEdges[l1].insert(mesh.Cell2DEdges[l1].begin()+k+1, id_edge_first_intersection);
                            }
                            mesh.Cell2DVertices[l1].insert(mesh.Cell2DVertices[l1].begin()+k+1, mesh.NumberCell0D); // added new vertex to that cell
                            break;
                        }
                    }

                }
                if(first_int_found) {
                    break;
                }
            }

            //update Cell0Ds
            mesh.Cell0DId.push_back(mesh.NumberCell0D);
            mesh.Cell0DCoordinates.push_back(first_intersection);

            // update Cell1Ds
            for(unsigned int k=0;k<2;k++) {
                mesh.Cell1DId.push_back(mesh.NumberCell1D+k);
            }
            mesh.Cell1DVertices[id_edge_first_intersection] = {mesh.Cell0DId[mesh.Cell2DVertices[l][iter[0]]],mesh.NumberCell0D};
            mesh.Cell1DVertices.push_back({mesh.NumberCell0D, mesh.Cell0DId[mesh.Cell2DVertices[l][(iter[0]+1)%mesh.Cell2DVertices[l].size()]]});
            mesh.Cell1DVertices.push_back({mesh.NumberCell0D,mesh.Cell0DId[mesh.Cell2DVertices[l][iter[1]]]});

            // update Cell2Ds
            mesh.Cell2DId.push_back(mesh.NumberCell2D);

            mesh.Cell2DVertices[l] = IdNewVertices1;
            mesh.Cell2DVertices.push_back(IdNewVertices2);

            mesh.Cell2DEdges[l] = IdNewEdges1;
            mesh.Cell2DEdges.push_back(IdNewEdges2);

            mesh.NumberCell0D += 1;
            mesh.NumberCell1D += 2;
            mesh.NumberCell2D += 1;


        }
    }
//*********************************************************
    void exportParaview(const string &outputFileName, const Fractures& FractureList) {
        Gedim::UCDUtilities exporter;

        vector<vector<unsigned int>> triangles;
        VectorXi materials;

        // find total number of points
        unsigned int num_points = 0;
        for(unsigned int i=0;i<FractureList.FractVertices.size();i++) {
            num_points += FractureList.FractVertices[i].cols();
        }
        MatrixXd VerticesCoordinates = MatrixXd::Zero(3,num_points);
        unsigned int temp = 0;
        for(unsigned int k=0;k<FractureList.FractVertices.size();k++) {
            for(unsigned int l=0;l<FractureList.FractVertices[k].cols();l++) {
                VerticesCoordinates.col(l+temp) = FractureList.FractVertices[k].col(l);
            }
            temp += FractureList.FractVertices[k].cols();
        }
        temp = 0;
        vector<vector<unsigned int>> listVertices;
        listVertices.resize(FractureList.FractVertices.size());
        for(unsigned int i=0;i<listVertices.size();i++) {
            const unsigned int n = FractureList.FractVertices[i].cols();
            listVertices[i].resize(n);
            for(unsigned int k=0;k<n;k++){
                listVertices[i][k] = k+temp;
            }
            temp+=n;
        }

        const unsigned int numPolygons = listVertices.size();
        vector<vector<vector<unsigned int>>> triangleList = TriangulatePolygons(listVertices);

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

        exporter.ExportPolygons(outputFileName,
                                VerticesCoordinates,
                                triangles,
                                {},
                                {},
                                materials);
    }
//*********************************************************
    void exportFractureMesh(const string& outputFileName, const Fractures& FractureList, const unsigned int& fractureID) {
        PolygonalMesh mesh = FractureList.FractMesh[fractureID];
        Gedim::UCDUtilities exporter;
        vector<vector<unsigned int>> triangles;
        VectorXi materials;
        FractMeshGedimInterface(mesh.Cell2DVertices, triangles, materials);
        const unsigned int numPoints = mesh.Cell0DCoordinates.size();
        MatrixXd CoordinatesMatrix = MatrixXd(3,numPoints);
        for(unsigned int i=0; i<numPoints;i++) {
            CoordinatesMatrix.col(i) = mesh.Cell0DCoordinates[i];
        }

        exporter.ExportPolygons(outputFileName,
                                CoordinatesMatrix,
                                triangles,
                                {},
                                {},
                                materials);

    }
//*********************************************************
    bool printFractureMesh(const string &outputFileName1, const string &outputFileName2, const string &outputFileName3, const Fractures &FractureList) {
        ofstream file1;
        file1.open(outputFileName1);
        if(file1.fail()) {
            cerr << "Something went wrong opening output file." << endl;
            return false;
        }
        ofstream file2;
        file2.open(outputFileName2);
        if(file2.fail()) {
            cerr << "Something went wrong opening output file." << endl;
            return false;
        }
        ofstream file3;
        file3.open(outputFileName3);
        if(file3.fail()) {
            cerr << "Something went wrong opening output file." << endl;
            return false;
        }
        file1 << scientific << setprecision(16);
        const unsigned int num_fractures = FractureList.FractMeanPoint.size();
        for(unsigned int i=0; i<num_fractures;i++) {
            file1 << "Fracture id;" << i << ";Number Cell0Ds;" << FractureList.FractMesh[i].NumberCell0D << ";" << endl;
            file1 << "Id;X;Y;Z;" << endl;
            for(unsigned int j=0;j<FractureList.FractMesh[i].NumberCell0D;j++) {
                file1 << j << ";" << FractureList.FractMesh[i].Cell0DCoordinates[j](0) <<";"
                      << FractureList.FractMesh[i].Cell0DCoordinates[j](1) <<";"
                      << FractureList.FractMesh[i].Cell0DCoordinates[j](2) <<";" << endl;
            }
            file1 << endl;

            file2 << "Fracture id;" << i << ";Number Cell1Ds;" << FractureList.FractMesh[i].NumberCell1D << ";" << endl;
            file2 << "Id;Origin;End;" << endl;
            for(unsigned int j=0;j<FractureList.FractMesh[i].NumberCell1D;j++) {
                file2 << j <<";" << FractureList.FractMesh[i].Cell1DVertices[j](0) << ";" << FractureList.FractMesh[i].Cell1DVertices[j](1) << ";" << endl;
            }
            file2 << endl;

            file3 << "Fracture id;" << i << ";Number Cell2Ds;" << FractureList.FractMesh[i].NumberCell2D << ";" << endl;
            file3 << "Id;NumVertices;Vertices;NumEdges;Edges;" << endl;
            for(unsigned int j=0; j<FractureList.FractMesh[i].NumberCell2D;j++) {
                file3 << j <<";" << FractureList.FractMesh[i].Cell2DVertices[j].size() << ";";
                for(unsigned int k=0;k<FractureList.FractMesh[i].Cell2DVertices[j].size();k++) {
                    file3 << FractureList.FractMesh[i].Cell2DVertices[j][k] << ";";
                }
                file3 << FractureList.FractMesh[i].Cell2DEdges[j].size() << ";";
                for(unsigned int k=0;k<FractureList.FractMesh[i].Cell2DEdges[j].size();k++) {
                    file3 << FractureList.FractMesh[i].Cell2DEdges[j][k] << ";";
                }
                file3 << endl;
            }
            file3 << endl;
        }
        return true;
    }
//*********************************************************
}

