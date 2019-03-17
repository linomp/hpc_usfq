/*
 * WallCylinder.h
 *
 *  Created on: Oct 24, 2014
 *      Author: Reid
 */

#ifndef WALLCYLINDER_H_
#define WALLCYLINDER_H_

#include "definitions.h"
#include "Grain3d.h"


class WallCylinder {
	
public:
	WallCylinder() { _height = 0; _radius = 0; _kn = 0; _mu = 0; _pressure = 0; _surfaceArea = 0; _id = 0;}
	WallCylinder(const double & height, const double & radius, const double & kn, const double & mu, const double & pressure, const size_t & id):
	_height(height), _radius(radius), _kn(kn), _mu(mu), _pressure(pressure), _id(id) {
		_surfaceArea = 2*M_PI*_radius*_height;
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
		
//		Vector3d normal(pos(0)/gRadial, pos(1)/gRadial, 0);
//		Vector3d v;
//		Vector3d ds;
//		Vector3d Fs;
//		double Fsmag;
		
		// naive check
		if (gRadial+gRadius > _radius) {
			// iterate through points
			for (size_t p = 0; p < grain.getPointList().size(); p++) {
				pos = grain.getPointList()[p];
				ptcm = pos - grain.getPosition();
				gRadial = sqrt(pos(0)*pos(0) + pos(1)*pos(1));
				penetration = gRadial - _radius;
				if (penetration > 0.) {
					ncontacts++;
					checkflag = true;
					df << -_kn*penetration*pos(0)/gRadial, -_kn*penetration*pos(1)/gRadial, 0;
					force += df;
					moment += -ptcm.cross(df);
					// friction stuff
//					v = grain.getVelocity() + grain.getOmegaGlobal().cross(ptcm); // eq (6)
//					ds = (v - v.dot(normal)*normal)*dt; // eq (7) and (9)
//					
//					grain.getNodeShearsNonConst()[p] -= ds*_kn; // technically ks
//					grain.getNodeContactNonConst()[p] = _id;
//					Fsmag = min(df.norm()*_mu, grain.getNodeShearsNonConst()[p].norm() );
//					
//					if (Fsmag > 0) {
//						Fs = Fsmag*grain.getNodeShearsNonConst()[p]/grain.getNodeShearsNonConst()[p].norm();
//						df += Fs;
//						fs += Fs;
//						moment  += ptcm.cross(Fs);
//					}
					
					stress(0) -= df(0)*ptcm(0);
					stress(1) -= df(1)*ptcm(1);
					stress(2) -= df(2)*ptcm(2);
					stress(3) -= 0.5*(df(1)*ptcm(2) + df(2)*ptcm(1));
					stress(4) -= 0.5*(df(2)*ptcm(0) + df(0)*ptcm(2));
					stress(5) -= 0.5*(df(1)*ptcm(0) + df(0)*ptcm(1));
				}
//				else if (grain.getNodeContactNonConst()[p] == _id ){
//					grain.getNodeContactNonConst()[p] = 0;
//					grain.getNodeShearsNonConst()[p] << 0,0,0;
//				}
			}
		}
		return checkflag;
	}
	
	
	void takeTimestep(const double & wallForce, const size_t & ncontacts) {
		// update _velocity and _radius using the wallForce
//		_velocity = 1/(1+gDamping*dt/2)*( (1-gDamping*dt/2)*_velocity + dt*(wallForce - _pressure*_surfaceArea)/_mass   );
//		_radius += dt*_velocity;
//		// compute the new _surfacearea
//		_surfaceArea = 2*M_PI*_radius*_height;
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
	double _radius;	// radius of wall
	double _kn;			// normal stiffness of wall
	double _mu;			// friction of wall
	double _pressure;	// radial pressure pushing the wall inwards
	double _surfaceArea; // surface area of cylinder (without top and bottom caps)
	size_t _id;
	
};



#endif /* WALLCYLINDER_H_ */
