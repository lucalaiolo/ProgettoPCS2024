#pragma once
#include "Fractures.hpp"
#include<iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;
using namespace DFNLibrary;

namespace SortingLibrary {

void MergeTraces(vector<unsigned int>& traces_ids, vector<double> traces_length,
           const unsigned int& sx,
           const unsigned int& cx,
           const unsigned int& dx){

    unsigned int i = sx;
    unsigned int j = cx + 1;

    vector<unsigned int> b;
    b.reserve(dx - sx + 1);

    while( i <= cx && j <= dx)
    {
        const double length1 = traces_length[traces_ids[i]]; // length of first trace
        const double length2 = traces_length[traces_ids[j]]; // length of second trace
        if (length2 <= length1)
            b.push_back(traces_ids[i++]);
        else
            b.push_back(traces_ids[j++]);
    }

    if (i <= cx)
        b.insert(b.end(), traces_ids.begin() + i, traces_ids.begin() + cx + 1);
    if ( j <= dx)
        b.insert(b.end(), traces_ids.begin() + j, traces_ids.begin() + dx + 1);

    copy(b.begin(), b.end(), traces_ids.begin() + sx);

}

void MergeSortTraces(vector<unsigned int>& traces_ids, vector<double> traces_length,
               const unsigned int& sx,
               const unsigned int& dx){

    if (sx < dx)
    {
        unsigned int cx = (sx + dx) / 2;
        MergeSortTraces(traces_ids,traces_length, sx, cx);
        MergeSortTraces(traces_ids,traces_length, cx + 1, dx);

        MergeTraces(traces_ids,traces_length, sx, cx, dx);
    }

}


void MergeSortTraces(vector<unsigned int>& traces_ids, const vector<double> traces_length){
    MergeSortTraces(traces_ids,traces_length, 0, traces_ids.size()-1);
}


}
