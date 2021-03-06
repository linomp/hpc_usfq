/*
 * Levelset3d.h
 *
 *  Created on: July 14, 2014
 *      Author: Reid
 */

#ifndef LEVELSET3D_H_
#define LEVELSET3D_H_

#include "definitions.h"

class Levelset3d {

public:
	// default constructor
	Levelset3d() {
		_xdim = 0; _ydim = 0; _zdim = 0;
	}
	// actual constructor
	Levelset3d(const vector<float> & levelset, const size_t & xdim, const size_t & ydim, const size_t & zdim):
		_levelset(levelset), _xdim(xdim), _ydim(ydim), _zdim(zdim) {
		if (_xdim*_ydim*_zdim != _levelset.size()) {
			cout << "ERROR: Level set size not consistent with dimensions" << endl;
		}
	}
	// Checks if there is penetration.  If there is, returns true and finds the penetration amount and the normalized gradient
	// and stores them in input fields value and gradient.  If not, returns false.
	bool findPenetration(const Vector3d & point, double & value, Vector3d & gradient) const {
		double x = point(0);
		double y = point(1);
		double z = point(2);
		// check if the point exists in the level set, if not, return false
		if (x+1 > (double)_xdim || y+1 > (double)_ydim || z+1 > (double)_zdim || x < 0 || y < 0 || z < 0 ){
			return false;
		}
		size_t xr = (size_t)round(x);
		size_t yr = (size_t)round(y);
		size_t zr = (size_t)round(z);
		// check if the point is close to the surface, if not, return false
		if (getGridValue(xr,yr,zr) > 1) {
			return false;
		} 
		// if the point is close to the surface, find the value by performing trilinear interpolation
		size_t x0 	= (size_t)floor(x);
		size_t y0 	= (size_t)floor(y);
		size_t z0 	= (size_t)floor(z);
		size_t x1 	= (size_t)ceil(x);
		size_t y1 	= (size_t)ceil(y);
		size_t z1 	= (size_t)ceil(z);
		double p000 = getGridValue(x0, y0, z0); 
		double p001 = getGridValue(x0,y0,z1);
		double p010 = getGridValue(x0, y1, z0);
		double p011 = getGridValue(x0,y1,z1);
		double p101 = getGridValue(x1,y0,z1);
		
		double xm 	= getGridValue(x1, y0, z0) - p000;
		double ym 	= p010 - p000;
		double zm	= p001 - p000;
		double xym	= -xm - p010 + getGridValue(x1,y1,z0);
		double xzm	= -xm - p001 + p101;
		double yzm	= -ym - p001 + p011;
		double xyzm = -xym + p001 - p101 - p011 + getGridValue(x1,y1,z1);
		double dx 	= x - double(x0);
		double dy 	= y - double(y0);
		double dz	= z - double(z0);
		value = p000 + xm*dx + ym*dy + zm*dz + xym*dx*dy + xzm*dx*dz + yzm*dy*dz + xyzm*dx*dy*dz;
		// if the point lies outside the surface, return false
		if (value > 0) {
			return false;
		}
//		cout << point.transpose() << endl;
		// finally, if the point lies in the surface, find the gradient, normalize, and return true
		gradient << xm + xym*dy + xzm*dz + xyzm*dy*dz, 
						ym + xym*dx + yzm*dz + xyzm*dx*dz, 
						zm + xzm*dx + yzm*dy + xyzm*dx*dy;
		gradient /= gradient.norm();
		return true;
	}
	
//	// public get methods (not actually used, for debugging only)
	double getXdim() const {
		return _xdim;
	}
//	double getYdim() const {
//		return _ydim;
//	}
//	vector<double> getLevelset() const {
//		return _levelset;
//	}
//	double getGridValue(size_t & x, size_t & y, size_t & z) const {
//		return _levelset[z*_ydim*_xdim + y*_xdim + x];
//	}
	
private:
	vector<float> _levelset;
	size_t _xdim;
	size_t _ydim;
	size_t _zdim;
	
	double getGridValue(size_t & x, size_t & y, size_t & z) const {
		return (double)_levelset[z*_ydim*_xdim + y*_xdim + x];
	}

};

#endif /* LEVELSET3D_H_ */
