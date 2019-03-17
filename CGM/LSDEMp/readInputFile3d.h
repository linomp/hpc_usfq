/*
 * readInputFile3d.h
 *
 *  Created on: July 15, 2014
 *      Author: Reid
 */

#ifndef READINPUTFILE3D_H_
#define READINPUTFILE3D_H_

#include "definitions.h"
#include "Levelset3d.h"
#include "Grain3d.h"
#include "Morphology3d.h"

extern vector<size_t> readIntegerFile(string filename, size_t num) {
	ifstream file(filename.c_str());
	string   line;
	vector<size_t> integers(num);
	for (size_t n = 0; n < num; n++) {
		getline(file, line);
		integers[n]= atoi(line.c_str());
	}
	return integers;
}

extern vector<double> readDoubleFile(string filename, size_t num) {
	ifstream file(filename.c_str());
	string   line;
	vector<double> doubles(num);
	for (size_t n = 0; n < num; n++) {
		getline(file, line);
		doubles[n]= atof(line.c_str());
	}
	return doubles;
}


extern vector<Vector3d> readPositionFile(string filename, size_t ngrains) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	vector<Vector3d> positions;
	positions.resize(ngrains);
	for (size_t grainidx = 0; grainidx < ngrains; grainidx++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		positions[grainidx](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		positions[grainidx](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		positions[grainidx](2) = atof(partial.c_str());
		iss.clear();
	}
	return positions;
} // end generateGrainsFromFile

extern vector<Vector4d> readQuaternionFile(string filename, size_t ngrains) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	vector<Vector4d> quaternions;
	quaternions.resize(ngrains);
	for (size_t grainidx = 0; grainidx < ngrains; grainidx++) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		quaternions[grainidx](0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		quaternions[grainidx](1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		quaternions[grainidx](2) = atof(partial.c_str());
		getline(iss, partial, ' ');
		quaternions[grainidx](3) = atof(partial.c_str());
		iss.clear();
	}
	return quaternions;
}

// get properties.  order: density, kn, ks, mu
extern vector<double> readPropertiesFile(string filename) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	vector<double> properties(4);
	properties.resize(4);
	for (size_t p = 0; p < 4; p++) {
		getline(file, line);
		properties[p]= atof(line.c_str());
	}
	return properties;
}


// creates a vector of grain objects from a file
extern vector<Grain3d> generateGrainsFromFile(string filename) {

	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	// read first line: number of particles
	getline(file, line);
	size_t numberOfGrains = atoi(line.c_str());

	vector<Grain3d> grainList(numberOfGrains);

	// temp stuff
	Vector3d position;
	Vector3d velocity;
	Vector3d momentOfInertia;
	Vector4d quat;
	Vector3d omega;
	Vector3d point;
	Vector3d cmLset;
	
	for (size_t grainidx = 0; grainidx < numberOfGrains; grainidx++) {
		// volume, also assumes a density of 1 so volume = mass, this can be chagned by changing the density
		getline(file, line);
        double mass = atof(line.c_str()); 

		// position
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		position(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		position(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		position(2) = atof(partial.c_str());
		iss.clear();

		// velocity
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		velocity(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		velocity(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		velocity(2) = atof(partial.c_str());
		iss.clear();

		// moment of inertia
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		momentOfInertia(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		momentOfInertia(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		momentOfInertia(2) = atof(partial.c_str());
		iss.clear();

		// theta (quaternion)
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		quat(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		quat(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		quat(2) = atof(partial.c_str());
		getline(iss, partial, ' ');
		quat(3) = atof(partial.c_str());
		iss.clear();

		// omega
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		omega(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		omega(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		omega(2) = atof(partial.c_str());
		iss.clear();

		// cmLset
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		cmLset(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		cmLset(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		cmLset(2) = atof(partial.c_str());
		iss.clear();
		
		// npoints (INTEGER)
		getline(file, line);
		int npoints = atoi(line.c_str());

		// points [x1, y1, z1, x2, y2, z2, etc]
		getline(file, line);
		vector<Vector3d> pointList(npoints);
		iss.str(line);
		for (int ptidx = 0; ptidx < npoints; ptidx++) {
			getline(iss, partial, '\t');
			point(0) = atof(partial.c_str());
			getline(iss, partial, '\t');
			point(1) = atof(partial.c_str());
			getline(iss, partial, '\t');
			point(2) = atof(partial.c_str());
			pointList[ptidx] = point;
		}
		iss.clear();

		// bbox radius
		getline(file, line);
		double bboxRadius = atof(line.c_str());

		// level set dimensions (INTEGERS)
		getline(file, line);
		iss.str(line);
		getline(iss, partial, '\t');
		int xdim = atoi(partial.c_str());
		getline(iss, partial, '\t');
		int ydim = atoi(partial.c_str());
		getline(iss, partial, '\t');
		int zdim = atoi(partial.c_str());
		iss.clear();

		// level set 
		getline(file, line);
		vector<float> lsetvec(xdim*ydim*zdim);
		iss.str(line);
		for (int i = 0; i < xdim*ydim*zdim; i++) {
			getline(iss, partial, '\t');
			lsetvec[i] = atof(partial.c_str());
		}
		iss.clear();

		// kn
		getline(file, line);
		double kn = atof(line.c_str());

		// ks
		getline(file, line);
		double ks = atof(line.c_str());

		// mu
		getline(file, line);
		double mu = atof(line.c_str());

		// create level set object
		Levelset3d * lset = new Levelset3d(lsetvec, xdim, ydim, zdim);
		
		// update grain object in the vector that was created at the beginning of this function
		grainList[grainidx] = Grain3d(mass, position, velocity, momentOfInertia, quat, omega, cmLset, pointList, bboxRadius, lset, kn, ks, mu, grainidx);
        
        grainList[grainidx].changeDensity(7.5346e-11); // hardcoded modification of grain density.
        
	}
	return grainList;
} // end generateGrainsFromFile


extern Grain3d generateGrainFromFiles(string morphfile, string posfile, string rotfile, string propfile, size_t grainidx) {

	ifstream file(morphfile.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	// read first line: number of particles
	getline(file, line);


	// temp stuff
	Vector3d momentOfInertia;
	Vector3d point;
	Vector3d cmLset;
	
	// mass
	getline(file, line);
    double mass = atof(line.c_str());

	// moment of inertia
	getline(file, line);
	iss.str(line);
	getline(iss, partial, ' ');
	momentOfInertia(0) = atof(partial.c_str());
	getline(iss, partial, ' ');
	momentOfInertia(1) = atof(partial.c_str());
	getline(iss, partial, ' ');
	momentOfInertia(2) = atof(partial.c_str());
	iss.clear();
	
	// cmLset
	getline(file, line);
	iss.str(line);
	getline(iss, partial, ' ');
	cmLset(0) = atof(partial.c_str());
	getline(iss, partial, ' ');
	cmLset(1) = atof(partial.c_str());
	getline(iss, partial, ' ');
	cmLset(2) = atof(partial.c_str());
	iss.clear();
	
	// npoints (INTEGER)
	getline(file, line);
	int npoints = atoi(line.c_str());

	// points
	getline(file, line);
	vector<Vector3d> pointList(npoints);
	iss.str(line);
	for (int ptidx = 0; ptidx < npoints; ptidx++) {
		getline(iss, partial, ' ');
		point(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		point(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		point(2) = atof(partial.c_str());
		pointList[ptidx] = point;
	}
	iss.clear();

	// bbox radius
	getline(file, line);
	double bboxRadius = atof(line.c_str());

	// level set dimensions (INTEGERS)
	getline(file, line);
	iss.str(line);
	getline(iss, partial, ' ');
	int xdim = atoi(partial.c_str());
	getline(iss, partial, ' ');
	int ydim = atoi(partial.c_str());
	getline(iss, partial, ' ');
	int zdim = atoi(partial.c_str());
	iss.clear();

	// level set
	getline(file, line);
	vector<float> lsetvec(xdim*ydim*zdim);
	iss.str(line);
	for (int i = 0; i < xdim*ydim*zdim; i++) {
		getline(iss, partial, ' ');
		lsetvec[i] = atof(partial.c_str());
	}
	iss.clear();
	
	vector<Vector3d> position = readPositionFile(posfile,1);
	vector<Vector4d> rotation = readQuaternionFile(rotfile,1);
	vector<double> props = readPropertiesFile(propfile);

	// create level set object
	Levelset3d * lset = new Levelset3d(lsetvec, xdim, ydim, zdim);
	
	// update grain object in the vector that was created at the beginning of this function
	Grain3d grain(mass, position[0], Vector3d(0,0,0), momentOfInertia, rotation[0], Vector3d(0,0,0), cmLset, pointList, bboxRadius, lset, props[1], props[2], props[3], grainidx);
	grain.changeDensity(props[0]);
		
	return grain;
} // end generateGrainsFromFile

// generates a grain without any position and rotation
extern Grain3d generateGrainFromFiles(string morphfile, string propfile, size_t grainidx) {

	ifstream file(morphfile.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	// read first line: number of particles
//	getline(file, line);

	// temp stuff
	Vector3d momentOfInertia;
	Vector3d point;
	Vector3d cmLset;
	
	// mass
	getline(file, line);
	double mass = atof(line.c_str());
    

	// moment of inertia
	getline(file, line);
	iss.str(line);
	getline(iss, partial, ' ');
	momentOfInertia(0) = atof(partial.c_str());
	getline(iss, partial, ' ');
	momentOfInertia(1) = atof(partial.c_str());
	getline(iss, partial, ' ');
	momentOfInertia(2) = atof(partial.c_str());
	iss.clear();
	
	// cmLset
	getline(file, line);
	iss.str(line);
	getline(iss, partial, ' ');
	cmLset(0) = atof(partial.c_str());
	getline(iss, partial, ' ');
	cmLset(1) = atof(partial.c_str());
	getline(iss, partial, ' ');
	cmLset(2) = atof(partial.c_str());
	iss.clear();
	
	// npoints (INTEGER)
	getline(file, line);
	int npoints = atoi(line.c_str());

	// points
	getline(file, line);
	vector<Vector3d> pointList(npoints);
	iss.str(line);
	for (int ptidx = 0; ptidx < npoints; ptidx++) {
		getline(iss, partial, ' ');
		point(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		point(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		point(2) = atof(partial.c_str());
		pointList[ptidx] = point;
	}
	iss.clear();

	// bbox radius
	getline(file, line);
	double bboxRadius = atof(line.c_str());

	// level set dimensions (INTEGERS)
	getline(file, line);
	iss.str(line);
	getline(iss, partial, ' ');
	int xdim = atoi(partial.c_str());
	getline(iss, partial, ' ');
	int ydim = atoi(partial.c_str());
	getline(iss, partial, ' ');
	int zdim = atoi(partial.c_str());
	iss.clear();

	// level set
	getline(file, line);
	vector<float> lsetvec(xdim*ydim*zdim);
	iss.str(line);
	for (int i = 0; i < xdim*ydim*zdim; i++) {
		getline(iss, partial, ' ');
		lsetvec[i] = atof(partial.c_str());
	}
	iss.clear();
	
	vector<double> props = readPropertiesFile(propfile);
	// create level set object
	Levelset3d * lset = new Levelset3d(lsetvec, xdim, ydim, zdim);
	// update grain object in the vector that was created at the beginning of this function
	Grain3d grain(mass, Vector3d(0,0,0), Vector3d(0,0,0), momentOfInertia, Vector4d(0,0,0,1), Vector3d(0,0,0), cmLset, pointList, bboxRadius, lset, props[1], props[2], props[3], grainidx);
	grain.changeDensity(props[0]);
	return grain;
} // end generateGrainsFromFile

// generates a pointer to a vector of morphologies
extern vector<Morphology3d> * generateMorphologiesFromFile(string filename) {

	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	// read first line: number of particles
	getline(file, line);
	size_t numberOfGrains = atoi(line.c_str());

	vector<Morphology3d> * grainList = new vector<Morphology3d>(numberOfGrains);

	// temp stuff
	Vector3d momentOfInertia;
	Vector3d point;
	Vector3d cmLset;
	
	for (size_t grainidx = 0; grainidx < numberOfGrains; grainidx++) {
		// mass
		getline(file, line);
		double mass = atof(line.c_str());
		
		// moment of inertia
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		momentOfInertia(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		momentOfInertia(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		momentOfInertia(2) = atof(partial.c_str());
		iss.clear();
		
		// cmLset
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		cmLset(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		cmLset(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		cmLset(2) = atof(partial.c_str());
		iss.clear();
		
		// npoints (INTEGER)
		getline(file, line);
		int npoints = atoi(line.c_str());
		
		// points
		getline(file, line);
		vector<Vector3d> pointList(npoints);
		iss.str(line);
		for (int ptidx = 0; ptidx < npoints; ptidx++) {
			getline(iss, partial, ' ');
			point(0) = atof(partial.c_str());
			getline(iss, partial, ' ');
			point(1) = atof(partial.c_str());
			getline(iss, partial, ' ');
			point(2) = atof(partial.c_str());
			pointList[ptidx] = point;
		}
		iss.clear();

		// bbox radius
		getline(file, line);
		double bboxRadius = atof(line.c_str());

		// level set dimensions (INTEGERS)
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		int xdim = atoi(partial.c_str());
		getline(iss, partial, ' ');
		int ydim = atoi(partial.c_str());
		getline(iss, partial, ' ');
		int zdim = atoi(partial.c_str());
		iss.clear();

		// level set
		getline(file, line);
		vector<float> lsetvec(xdim*ydim*zdim);
		iss.str(line);
		for (int i = 0; i < xdim*ydim*zdim; i++) {
			getline(iss, partial, ' ');
			lsetvec[i] = atof(partial.c_str());
		}
		iss.clear();
		
		// create level set object
		Levelset3d lset(lsetvec, xdim, ydim, zdim);
		
		// update grain object in the vector that was created at the beginning of this function
		grainList->at(grainidx) = Morphology3d(mass, momentOfInertia, cmLset, pointList, bboxRadius, lset, grainidx);
	}
	return grainList;
} // end generateGrainsFromFile

void readShearHistFile(string filename, vector<Grain3d> & grains) {
	ifstream file(filename.c_str());
	string   line;
	string 	partial;
	istringstream iss;
	size_t grainidx;
	size_t nodeidx;
	size_t contactingidx;
	Vector3d shear;
	Vector3d normal;
	while (true) {
		getline(file, line);
		iss.str(line);
		getline(iss, partial, ' ');
		grainidx = atoi(partial.c_str());
		getline(iss, partial, ' ');
		nodeidx = atoi(partial.c_str());
		getline(iss, partial, ' ');
		contactingidx = atoi(partial.c_str());
		getline(iss, partial, ' ');
		shear(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		shear(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		shear(2) = atof(partial.c_str());
		getline(iss, partial, ' ');
		normal(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		normal(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		normal(2) = atof(partial.c_str());
		grains[grainidx].changeShearHist(nodeidx, contactingidx, shear, normal);
		iss.clear();
		if (file.eof()) {
			break;
		}
	}
}

// generates a vector of grains with only base morphologies
//extern vector<Morphology3d> generateMorphologiesFromFile2(string filename) {
//
//	ifstream file(filename.c_str());
//	string   line;
//	string 	partial;
//	istringstream iss;
//	// read first line: number of particles
//	getline(file, line);
//	size_t numberOfGrains = atoi(line.c_str());
//
//	vector<Morphology3d> grainList(numberOfGrains);
//
//	// temp stuff
//	Vector3d momentOfInertia;
//	Vector3d point;
//	Vector3d cmLset;
//	
//	for (size_t grainidx = 0; grainidx < numberOfGrains; grainidx++) {
//		// mass
//		getline(file, line);
//		double mass = atof(line.c_str());
//		
//		// moment of inertia
//		getline(file, line);
//		iss.str(line);
//		getline(iss, partial, ' ');
//		momentOfInertia(0) = atof(partial.c_str());
//		getline(iss, partial, ' ');
//		momentOfInertia(1) = atof(partial.c_str());
//		getline(iss, partial, ' ');
//		momentOfInertia(2) = atof(partial.c_str());
//		iss.clear();
//		
//		// cmLset
//		getline(file, line);
//		iss.str(line);
//		getline(iss, partial, ' ');
//		cmLset(0) = atof(partial.c_str());
//		getline(iss, partial, ' ');
//		cmLset(1) = atof(partial.c_str());
//		getline(iss, partial, ' ');
//		cmLset(2) = atof(partial.c_str());
//		iss.clear();
//		
//		// npoints (INTEGER)
//		getline(file, line);
//		int npoints = atoi(line.c_str());
//		
//		// points
//		getline(file, line);
//		vector<Vector3d> pointList(npoints);
//		iss.str(line);
//		for (int ptidx = 0; ptidx < npoints; ptidx++) {
//			getline(iss, partial, ' ');
//			point(0) = atof(partial.c_str());
//			getline(iss, partial, ' ');
//			point(1) = atof(partial.c_str());
//			getline(iss, partial, ' ');
//			point(2) = atof(partial.c_str());
//			pointList[ptidx] = point;
//		}
//		iss.clear();
//
//		// bbox radius
//		getline(file, line);
//		double bboxRadius = atof(line.c_str());
//
//		// level set dimensions (INTEGERS)
//		getline(file, line);
//		iss.str(line);
//		getline(iss, partial, ' ');
//		int xdim = atoi(partial.c_str());
//		getline(iss, partial, ' ');
//		int ydim = atoi(partial.c_str());
//		getline(iss, partial, ' ');
//		int zdim = atoi(partial.c_str());
//		iss.clear();
//
//		// level set
//		getline(file, line);
//		vector<float> lsetvec(xdim*ydim*zdim);
//		iss.str(line);
//		for (int i = 0; i < xdim*ydim*zdim; i++) {
//			getline(iss, partial, ' ');
//			lsetvec[i] = atof(partial.c_str());
//		}
//		iss.clear();
//		
//		// create level set object
//		Levelset3d lset(lsetvec, xdim, ydim, zdim);
//		
//		// update grain object in the vector that was created at the beginning of this function
//		grainList[grainidx] = Morphology3d(mass, momentOfInertia, cmLset, pointList, bboxRadius, lset, grainidx);
//	}
//	return grainList;
//} // end generateGrainsFromFile



//extern vector<double> readWallPosFile(string filename) {
//	vector<double> wallPos;
//	
//	return wallPos;
//}

//extern array<double,2> readWallPosFile(string filename) {
//	ifstream file(filename.c_str());
//	string   line;
//	string 	partial;
//	istringstream iss;
//	array<double,2> wallPos;
//	
//	
//}














#endif /* READINPUTFILE3D_H_ */
