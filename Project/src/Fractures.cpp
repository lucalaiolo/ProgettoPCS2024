#include "Fractures.hpp"
#include "Sorting.hpp"
#include "UCDUtilities.hpp"
#include<iostream>
#include "Eigen/Eigen"
#include <fstream>
#include <iomanip>
using namespace std;
using namespace Eigen;
using namespace SortingLibrary;
namespace DFNLibrary {

    void computeTraces(Fractures &FractureList, Traces &TracesList) {
        const unsigned int numFractures = FractureList.FractVertices.size();
        unsigned int count = 0; // variable that will tell us the id of the computed traces
        for(unsigned int i=0;i<numFractures;i++) {
            for(unsigned int j=i+1;j<numFractures;j++){ // intersection is commutative
                // cout << "Fracture " << i<< endl;
                // cout << "Fracture " << j<< endl << endl;
                // we must check if it's possible that two fractures can intersect

                // const vector<unsigned int> Fract1_list_vertices_id = listVertices[i];
                // const vector<unsigned int> Fract2_list_vertices_id = listVertices[j];
                // const unsigned int num_vertices_frac1 = Fract1_list_vertices_id.size();
                // const unsigned int num_vertices_frac2 = Fract2_list_vertices_id.size();

                // const MatrixXd Fract1_vertices = MatrixXd::Zero(3,num_vertices_frac1);
                // const MatrixXd Fract2_vertices = MatrixXd::Zero(3,num_vertices_frac2);
                // for(unsigned int k=0;k<num_vertices_frac1;k++) {
                //      Fract1_vertices.col(k)=listVertices[Fract1_list_vertices_id[k]];
                // }
                // for(unsigned int k=0;k<num_vertices_frac2;k++) {
                //      Fract2_vertices.col(k)=listVertices[Fract2_list_vertices_id[k]];
                // }

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


                const int n = findCase(betas);

                if(n==0) {
                    //cout << "Trace is A_B. False for both." << endl;
                    array<bool,2> tips = {false,false};
                    executeCase(count,i,j,0,1,A_B_C_D,tips,FractureList,TracesList);
                    continue;
                } else if(n == 1 || n == 5) {
                    //cout << "No trace" << endl;
                    continue;
                } else if(n==-1 || n==2) {
                    //cout << "Trace is A_B. False F1. True F2." << endl;
                    array<bool,2> tips = {false,true};
                    executeCase(count,i,j,0,1,A_B_C_D,tips,FractureList,TracesList);
                    continue;
                } else if(n==3) {
                    //cout << "Case 3. Trace is CB. True both." << endl;
                    array<bool,2> tips = {true,true};
                    executeCase(count,i,j,2,1,A_B_C_D,tips,FractureList,TracesList);
                    continue;
                } else if(n==4) {
                    //cout << "Case 4. Trace is AD. True both." << endl;
                    array<bool,2> tips = {true,true};
                    executeCase(count,i,j,0,3,A_B_C_D,tips,FractureList,TracesList);
                    continue;
                } else if(n==6 || n==-2 || n==-3) {
                    //cout << "Case 6. Trace is CD. F1 true. F2 false." << endl;
                    array<bool,2> tips = {true,false};
                    executeCase(count,i,j,2,3,A_B_C_D,tips,FractureList,TracesList);
                    continue;
                } else if(n==-4) {
                    //cout << "Case 7.2" << endl;
                    //cout << "Trace is AD. false F1. true F2" << endl;
                    array<bool,2> tips = {false,true};
                    executeCase(count,i,j,0,3,A_B_C_D,tips,FractureList,TracesList);
                    continue;
                }

            }
        }
        computeTracesSquaredLength(TracesList);
        exportTraces("Traces.txt", TracesList);
        exportFractures("Fractures.txt", FractureList, TracesList);

    }

    bool importFractureList(const string &filepath, Fractures &FractureList, Traces& TracesList) {
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
        computeTraces(FractureList,TracesList);
        return true;
    }

    inline double computeSquaredDistancePoints(const Vector3d Point1, const Vector3d Point2) {
        double dist = 0.0;
        for(unsigned int i=0;i<3;i++) {
            dist += (Point1(i)-Point2(i))*(Point1(i)-Point2(i));
        }
        return dist;
    }

    inline void executeCase(unsigned int &count, const unsigned int &i, const unsigned int &j,
                            const unsigned int &pos1, const unsigned int &pos2, const vector<Vector3d> &A_B_C_D,
                            const array<bool, 2> &tips, Fractures &FractureList, Traces &TracesList) {
        MatrixXd trace_coord = MatrixXd::Zero(3,2);
        trace_coord.col(0) = A_B_C_D[pos1];
        trace_coord.col(1) = A_B_C_D[pos2];
        //
        //
        // MIGLIORARE QUI PUSH BACK
        //
        //
        auto ret1 = FractureList.listTraces[i].insert({tips[0], {count}});
        if(!ret1.second)
            (*(ret1.first)).second.push_back(count);
        auto ret2 = FractureList.listTraces[j].insert({tips[1], {count}});
        if(!ret2.second)
            (*(ret2.first)).second.push_back(count);
        array<unsigned int,2> trace_id = {i,j};
        TracesList.TraceIDFractures.push_back(trace_id);
        TracesList.TraceCoordinates.push_back(trace_coord);
        count++;
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

    int findCase(const vector<double> &betas) {
        // we first check if any of these points coincide
        const double tol = 2.2e-16;
        const bool A_equals_C = (fabs(betas[0]-betas[2])/max(max(fabs(betas[0]),fabs(betas[2])),{1}) < tol);
        //bool A_equals_D = (fabs(betas[0]-betas[3])/max(max(fabs(betas[0]),fabs(betas[3])),{1}) < tol); // not possible under our hypotheses
        // bool B_equals_C = (fabs(betas[1]-betas[2])/max(max(fabs(betas[1]),fabs(betas[2])),{1}) < tol); // not possible under our hypotheses
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
            vector<unsigned int> non_tips_traces = FractureList.listTraces[k][false];
            if(non_tips_traces.size() > 0) {
                MergeSortTraces(non_tips_traces, TracesList.TraceLength); // sort traces by their length
                for(unsigned int l=0;l<non_tips_traces.size();l++) {
                    file << non_tips_traces[l] << "; " << false << "; " << scientific << setprecision(16) << sqrt(TracesList.TraceLength.at(non_tips_traces[l])) << ";" << endl;
                }
            }
            vector<unsigned int> tips_traces = FractureList.listTraces[k][true];
            if(tips_traces.size() > 0) {
                MergeSortTraces(tips_traces, TracesList.TraceLength); // sort traces by their length
                for(unsigned int l=0;l<tips_traces.size();l++) {
                    file << tips_traces[l] << "; " << true << "; " << scientific << setprecision(16) << sqrt(TracesList.TraceLength[tips_traces[l]]) << ";" << endl;
                }
            }
        }

    }

    void computeTracesSquaredLength(Traces &TracesList) {
        const unsigned int num_traces = TracesList.TraceCoordinates.size();
        TracesList.TraceLength.reserve(num_traces);
        for(unsigned int k=0;k<num_traces;k++) {
            const Vector3d P1 = TracesList.TraceCoordinates[k].col(0);
            const Vector3d P2 = TracesList.TraceCoordinates[k].col(1); // points of first trace
            TracesList.TraceLength.push_back(computeSquaredDistancePoints(P1,P2)); // squared length of trace
        }
    }


    void exportParaview(const string &outputFileName, const Fractures& FractureList) {
        Gedim::UCDUtilities exporter;
        vector<vector<unsigned int>> quadrilaterals;
        quadrilaterals.resize(FractureList.FractVertices.size());
        for(unsigned int k=0;k<quadrilaterals.size();k++) {
            quadrilaterals[k] = {4*k,4*k+1,4*k+2,4*k+3};
        }
        VectorXi materials = VectorXi::Zero(FractureList.FractVertices.size());
        for(unsigned int k=0;k<materials.size();k++) {
            materials(k) = k;
        }
        MatrixXd VerticesCoordinates = MatrixXd::Zero(3,FractureList.FractVertices.size()*4);
        for(unsigned int k=0;k<FractureList.FractVertices.size();k++) {
            for(unsigned int l=0;l<FractureList.FractVertices[k].cols();l++) {
                VerticesCoordinates.col(4*k+l) = FractureList.FractVertices[k].col(l);
            }
        }
        exporter.ExportPolygons(outputFileName,
                                VerticesCoordinates,
                                quadrilaterals,
                                {},
                                {},
                                materials);
    }
}

