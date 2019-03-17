/*
 * generationMethods.h
 *
 *  Created on: Sep 16, 2015
 *      Author: reid
 */

#ifndef GENERATIONMETHODS_H_
#define GENERATIONMETHODS_H_

#include "definitions.h"
#include "Levelset3d.h"
#include "Grain3d.h"
#include "Morphology3d.h"
 
//vector<Morphology3d> morphs;
extern vector<Grain3d> genGrainsFromFiles(string morphFile, string posFile, string rotFile, string propFile, size_t ngrains) {
	vector<Grain3d> grains(ngrains);
	
	vector<Morphology3d> * morphs = generateMorphologiesFromFile(morphFile);
	
	size_t nmorphs = morphs->size();
//	cout << "number of morphs: " << nmorphs << endl;
	vector<Vector3d> positions = readPositionFile(posFile,ngrains);
	vector<Vector4d> quaternions = readQuaternionFile(rotFile,ngrains);
	vector<double> props = readPropertiesFile(propFile);
	
	for (size_t i = 0; i < ngrains; i++) {
		size_t n = i%nmorphs;
		grains[i] = Grain3d(morphs->at(n).getMass(), positions[i], Vector3d(0,0,0), morphs->at(n).getMomentInertia(), quaternions[i], Vector3d(0,0,0), 
								  morphs->at(n).getCmLset(), morphs->at(n).getRefPointList(), morphs->at(n).getRadius(), &(morphs->at(n).getLset()), props[1], props[2], props[3], i);
		grains[i].changeDensity(props[0]);
	}
	return grains;
}

//extern vector<Grain3d> genGrainsFromFiles(string morphFile, vector<Vector3d> positions, vector<Vector4d> quaternions, string propFile, size_t ngrains) {
//	vector<Grain3d> grains(ngrains);
//	
//	vector<Morphology3d> * morphs = generateMorphologiesFromFile(morphFile);
//	
//	size_t nmorphs = morphs->size();
////	cout << "number of morphs: " << nmorphs << endl;
////	vector<Vector3d> positions = readPositionFile(posFile,ngrains);
//	vector<Vector4d> quaternions = readQuaternionFile(rotFile,ngrains);
//	vector<double> props = readPropertiesFile(propFile);
//	
//	for (size_t i = 0; i < ngrains; i++) {
//		size_t n = i%nmorphs;
//		grains[i] = Grain3d(morphs->at(n).getMass(), positions[i], Vector3d(0,0,0), morphs->at(n).getMomentInertia(), quaternions[i], Vector3d(0,0,0), 
//								  morphs->at(n).getCmLset(), morphs->at(n).getRefPointList(), morphs->at(n).getRadius(), &(morphs->at(n).getLset()), props[1], props[2], props[3], i);
//		grains[i].changeDensity(props[0]);
//	}
//	return grains;
//}

//extern vector<Grain3d> genGrainsFromFiles2(string morphFile, string posFile, string rotFile, string propFile, size_t ngrains) {
//	vector<Grain3d> grains(ngrains);
//	
//	vector<Morphology3d> morphs = generateMorphologiesFromFile2(morphFile);
//	
//	size_t nmorphs = morphs.size();
////	cout << "number of morphs: " << nmorphs << endl;
//	vector<Vector3d> positions = readPositionFile(posFile,ngrains);
//	vector<Vector4d> quaternions = readQuaternionFile(rotFile,ngrains);
//	vector<double> props = readPropertiesFile(propFile);
//	
//	for (size_t i = 0; i < ngrains; i++) {
//		size_t n = i%nmorphs;
//		grains[i] = Grain3d(morphs[n].getMass(), positions[i], Vector3d(0,0,0), morphs[n].getMomentInertia(), quaternions[i], Vector3d(0,0,0), 
//								  morphs[n].getCmLset(), morphs[n].getRefPointList(), morphs[n].getRadius(), &(morphs[n].getLset()), props[1], props[2], props[3], i);
//		grains[i].changeDensity(props[0]);
//	}
//	return grains;
//}


#endif /* GENERATIONMETHODS_H_ */
