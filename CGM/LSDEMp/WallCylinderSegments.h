/*
 * WallCylinderSegments.h
 *
 *  Created on: Apr 15, 2015
 *      Author: Reid
 */

#ifndef WALLCYLINDERSEGMENTS_H_
#define WALLCYLINDERSEGMENTS_H_


#include "definitions.h"
#include "Grain3d.h"


class WallCylinderSegments {
	
public:
	WallCylinderSegments() { _height = 0; _kn = 0; _mu = 0; _pressure = 0; _surfaceArea = 0; _nsegs = 0;}
	WallCylinderSegments(const double & height, const double & radius, const double & kn, const double & mu, const double & pressure, const size_t & nsegs):
	_height(height), _kn(kn), _mu(mu), _pressure(pressure), _nsegs(nsegs) {
		_surfaceArea = 2*M_PI*radius*_height;
		_segRadii.resize(_nsegs);
		for (size_t i = 0; i < _nsegs; i++) {
			_segRadii[i] = radius;
		}
	}
	
	 
	bool bCircleContact(const Grain3d & grain) const {
		Vector3d pos; double gRadius, gRadial;
		pos = grain.getPosition();
		gRadius = grain.getRadius();
		gRadial = sqrt(pos(0)*pos(0) + pos(1)*pos(1));
		if (gRadial+gRadius > _radius) {
			return true;
		}
		return false;
	}
	
	bool findWallForce(Grain3d & grain, Vector3d & force, Vector3d & moment, size_t & ncontacts, Vector6d & stress) const {
		// zero the inputs
		force << 0,0,0; // radial force pushing the wall outwards
		moment << 0,0,0;
		ncontacts = 0;			// number of wall contacts for this particular grain
		stress << 0,0,0,0,0,0;
		// temp vars
		bool checkflag = false;
		double gRadial;	// radial component of grain's x-y coordinates
//		double gTheta;		// theta component of grain's x-y coordinates
//		double gHeight;	// height of grain's z-coordinate
		Vector3d pos;		// position of grain points
		double gRadius;	// radius of the grain
		double penetration; // penetration amount into wall
		Vector3d ptcm;		// location of point wrt the grain center of mass
		Vector3d df;	// incremental force on the grain
		pos = grain.getPosition();
		gRadius = grain.getRadius();
		gRadial = sqrt(pos(0)*pos(0) + pos(1)*pos(1));
		
		
		// naive check
		if (gRadial+gRadius > _radius) {
			// iterate through points
			for (size_t p = 0; p < grain.getPointList().size(); p++) {
				pos = grain.getPointList()[p];
				ptcm = pos - grain.getPosition();
				penetration = sqrt(pos(0)*pos(0) + pos(1)*pos(1)) - _radius;
				if (penetration > 0.) {
					ncontacts++;
					checkflag = true;
					df << -_kn*penetration*pos(0)/gRadial, -_kn*penetration*pos(1)/gRadial, 0;
					force += df;
					moment += -ptcm.cross(df);
					
					stress(0) -= df(0)*ptcm(0);
					stress(1) -= df(1)*ptcm(1);
					stress(2) -= df(2)*ptcm(2);
					stress(3) -= 0.5*(df(1)*ptcm(2) + df(2)*ptcm(1));
					stress(4) -= 0.5*(df(2)*ptcm(0) + df(0)*ptcm(2));
					stress(5) -= 0.5*(df(1)*ptcm(0) + df(0)*ptcm(1));
				}

			}
		}
		return checkflag;
	}
	
	
	void takeTimestep(const double & wallForce, const size_t & ncontacts) {

		static const double alpha = 0.7;
		if (ncontacts < 5)
			_radius += alpha*(wallForce-_pressure*_surfaceArea)/_kn/5.;
		else
			_radius += alpha*(wallForce-_pressure*_surfaceArea)/_kn/double(ncontacts);
		_surfaceArea = 2*M_PI*_radius*_height;
	}
	
	// edit methods
	void moveHeight(const double & amount) {
		_height += amount;
		_surfaceArea = 2*M_PI*_radius*_height;
	}
	void changeHeight(const double & newHeight) {
		_height = newHeight;
		_surfaceArea = 2*M_PI*_radius*_height;
	}
	void changePressure(const double & pressure) {
		_pressure = pressure;
	}
	void moveRadius(const double & amount) {
		_radius += amount;
		_surfaceArea = 2*M_PI*_radius*_height;
	}
	void changeRadius(const double & radius) {
		_radius = radius;
		_surfaceArea = 2*M_PI*_radius*_height;
	}
	void changeKn(const double & newkn) {
		_kn = newkn;
	}
	void changeMu(const double & mu) {
		_mu = mu;
	}
	
	 
	
	// get methods
	const double & getHeight() const {
		return _height;
	}
	const double & getRadius() const {
		return _radius;
	}
	
	
private:
	
	double _height;	// height of wall
	double _kn;			// normal stiffness of wall
	double _mu;			// friction of wall
	double _pressure;	// radial pressure pushing the wall inwards
	double _surfaceArea; // surface area of cylinder (without top and bottom caps)
	size_t _nsegs;		// number of segments
	vector<double> _segRadii; // radius of each segment
	
};


#endif /* WALLCYLINDERSEGMENTS_H_ */
