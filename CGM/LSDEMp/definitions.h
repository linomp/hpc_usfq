/*
 * definitions.h
 *
 *  Created on: May 30, 2014
 *      Author: Reid
 */

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

// needed for Eigen to work with Windows (who uses windows anymore?)
//#ifdef _WIN32
//	#define EIGEN_DONT_ALIGN_STATICALLY
//#endif

// C libraries
#include <float.h>		// DBL_MAX, INT_MAX
#include <stdio.h>		// printf, fgets
#include <stdlib.h>		// atoi
#include <omp.h>			// OpenMP
#include <mpi.h>			// MPI

// C++ libraries
#include <cmath>
using std::sqrt;
using std::sin;
using std::cos;
using std::isnan; // already included????

#include <cstdlib>
using std::size_t;

#include <string>		// strings
using std::string;

#include <iostream>		// terminal output
using std::cout;
using std::endl;

#include <fstream>		// file io
using std::stringstream;
using std::ifstream;
using std::getline;
using std::istringstream;

#include <algorithm>	// math stuff
using std::min;
using std::max;

#include <vector>		// standard vector
using std::vector;

#include<tr1/array>
using std::tr1::array;	// standard array

#include <Eigen/Dense> 
using Eigen::MatrixXd; // shellWalls
using Eigen::Vector3d; // position, velocity, etc
using Eigen::Vector2i; // for some stuff with 2 integers as elements
using Eigen::Vector3i; // 3-integer vector
using Eigen::Vector4d; // quaternions
using Eigen::Matrix3d; // rotation matrices
typedef Eigen::Matrix<double,1,6> Vector6d; // voigt stress 
typedef Eigen::Matrix<int,1,6> Vector6i; // 6-integer vector









#endif /* DEFINITIONS_H_ */
