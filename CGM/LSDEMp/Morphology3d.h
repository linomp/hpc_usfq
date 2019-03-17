/*
 * Morphology3d.h
 *
 *  Created on: Aug 26, 2015
 *      Author: reid
 */

#ifndef MORPHOLOGY3D_H_
#define MORPHOLOGY3D_H_

#include "definitions.h"
#include "Levelset3d.h"



class Morphology3d {
public:
	// constructors
	Morphology3d() {
		_radius = 0; _mass = 0; _id = 0;
	}

	Morphology3d(const double & mass, const Vector3d & momentInertia, const Vector3d & cmLset, 
			        const vector<Vector3d> & refPointList, const double & radius, const Levelset3d & lset, 
			        const size_t & id):
			  _mass(mass), _momentInertia(momentInertia), _cmLset(cmLset), _refPointList(refPointList), 
			  _radius(radius), _lset(lset), _id(id) {
		
		// clean up grain because i'm an idiot and my characterization method sucks (sometimes)
		for (size_t p = 0; p < _refPointList.size(); p++) {
			if (_refPointList[p].norm() > 1e3) {
				double radius = 0;
				_refPointList[p] << 0,0,0;
				for (size_t q = 0; q < _refPointList.size(); q++) {
					if (_refPointList[q].norm() > radius) {
						radius = _refPointList[q].norm();
					}
				}
				_radius = radius;
//				cout << _radius << endl;
			}
		}
	}
	
	// get methods
	const double & getMass() const {
		return _mass;
	}
	const Vector3d & getMomentInertia() const {
		return _momentInertia;
	}
	const Vector3d & getCmLset() const {
		return _cmLset;
	}
	const double & getRadius() const {
		return _radius;
	}
	Levelset3d & getLset() {
		return _lset;
	}
	const vector<Vector3d> & getRefPointList() const {
		return _refPointList;
	}
	const size_t & getId() const {
		return _id;
	}
	 
	// non const
//	vector<Vector3d> & getNodeShearsNonConst() {
//		return _nodeShears;
//	}
//	vector<size_t> & getNodeContactNonConst() {
//		return _nodeContact;
//	}
	
private:
	double 		_mass;
	Vector3d 	_momentInertia;// moment of inertia in principal frame (purely diagonal terms)
	Vector3d		_cmLset; 		// center of mass wrt the level set reference configuration
	vector<Vector3d> _refPointList; 	// list of points in reference config (center of mass is at (0,0,0) and I is diagonal)
	double 		_radius; // radius of bounding sphere
	Levelset3d 	_lset;	// level set of grain
	size_t _id;
	
};


#endif /* MORPHOLOGY3D_H_ */
