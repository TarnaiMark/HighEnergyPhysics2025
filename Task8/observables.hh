#pragma once
#pragma warning(push)
#pragma warning(disable : 4127)

// used headers and/or libraries
// custom
#include "auxiliary.hh"

using namespace std;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// YOUR CODE GOES HERE
// Polyakov loop
complex<double> PolyakovLoop(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims, vector<int> const &spatialSite);

// average Polyakov loop
complex<double> AveragePolyakovLoop(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims);

#pragma warning(pop)
