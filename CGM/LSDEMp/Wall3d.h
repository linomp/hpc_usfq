/*
 * Wall3d.h
 *
 *  Created on: August 25, 2014
 *      Author: Reid
 */

#ifndef WALL3D_H_
#define WALL3D_H_

#include "definitions.h"
#include "Grain3d.h"

class Wall3d {
	
public:
	Wall3d() {
		_d = 0; _kn = 0; _mu = 0; _id = INT_MAX-1;
	}
	Wall3d(const Vector3d & normal, const Vector3d & position, const double & kn, const double & mu): 
		_normal(normal), _position(position), _kn(kn), _mu(mu)  {
		_normal = _normal/_normal.norm();
		_d = -(normal.dot(position) );
		_velocity << 0,0,0;
		_id = INT_MAX-1;
		_ks = 0.8*_kn;
	}
	
	bool bCircleContact(const Grain3d & grain) const {
		if	 ( _normal.dot(grain.getPosition()) + _d > grain.getRadius() ) {
			return false;
		}
		return true;
	}
	
	// frictionless walls - for now
	bool findWallForceMoment(const Grain3d & grain, Vector3d & force, Vector3d & moment, size_t & ncontacts, Vector6d & stress) const {
		// zero the input vars
		force << 0., 0., 0.;
		moment << 0., 0., 0.;
		stress << 0,0,0,0,0,0;
		ncontacts = 0;
//		if	 ( _normal.dot(grain.getPosition()) + _d > grain.getRadius() ) {
//			return false;
//		}
		bool 		checkflag = false;	// changes to true if penetration exists
		double 	penetration; 			// penetration amount
		Vector3d df;						// force increment on grain
		Vector3d v;							// relative velocity
		Vector3d sdot;						// projection of relative velocity into tangential direction
		Vector3d Fs;
		Vector3d ptcm;
		for (size_t ptidx = 0; ptidx < grain.getPointList().size(); ptidx++) {
			penetration = _normal.dot(grain.getPointList()[ptidx]) + _d;
			if ( penetration < 0 ) {
				ncontacts++;
				ptcm = grain.getPointList()[ptidx] - grain.getPosition();
				ptcm = ptcm + penetration*ptcm/ptcm.norm();
				checkflag = true;
				df = -penetration*_normal*_kn;
				force += df;
				moment += ptcm.cross(df);
				
				stress(0) -= df(0)*ptcm(0);
				stress(1) -= df(1)*ptcm(1);
				stress(2) -= df(2)*ptcm(2);
				stress(3) -= 0.5*(df(1)*ptcm(2) + df(2)*ptcm(1));
				stress(4) -= 0.5*(df(2)*ptcm(0) + df(0)*ptcm(2));
				stress(5) -= 0.5*(df(1)*ptcm(0) + df(0)*ptcm(1));
				
				// friction stuff
//				v = grain.getVelocity() + grain.getOmegaGlobal().cross(ptcm); // eq (6)
//				ds = (v - v.dot(normal)*normal)*dt; // eq (7) and (9)
//				
//				grain.getNodeShearsNonConst()[p] -= ds*_kn; // technically ks
//				grain.getNodeContactNonConst()[p] = _id;
//				Fsmag = min(df.norm()*_mu, grain.getNodeShearsNonConst()[p].norm() );
//				
//				if (Fsmag > 0) {
//					Fs = Fsmag*grain.getNodeShearsNonConst()[p]/grain.getNodeShearsNonConst()[p].norm();
//					df += Fs;
//					fs += Fs;
//					moment  += ptcm.cross(Fs);
//				}

			}
//			else if (grain.getNodeContactNonConst()[p] == _id ){
//				grain.getNodeContactNonConst()[p] = 0;
//				grain.getNodeShearsNonConst()[p] << 0,0,0;
//			}
		}
		return checkflag;
	}
	
	// frictionless walls - for now
	bool findWallForceMomentFriction(Grain3d & grain, Vector3d & force, Vector3d & moment, size_t & ncontacts, Vector6d & stress, const double & dt) {
		// zero the input vars
		force << 0., 0., 0.;
		moment << 0., 0., 0.;
		stress << 0,0,0,0,0,0;
		ncontacts = 0;
//		if	 ( _normal.dot(grain.getPosition()) + _d > grain.getRadius() ) {
//			return false;
//		}
		bool 		checkflag = false;	// changes to true if penetration exists
		double 	penetration; 			// penetration amount
		Vector3d df;						// force increment on grain
		Vector3d v;							// relative velocity
		Vector3d sdot;						// projection of relative velocity into tangential direction
		Vector3d Fs;
		Vector3d ptcm;
		Vector3d ds;
		double Fsmag;
		Vector3d fs;
		for (size_t ptidx = 0; ptidx < grain.getPointList().size(); ptidx++) {
			penetration = _normal.dot(grain.getPointList()[ptidx]) + _d;
			if ( penetration < 0 ) {
				ncontacts++;
				ptcm = grain.getPointList()[ptidx] - grain.getPosition();
				ptcm = ptcm + penetration*ptcm/ptcm.norm();
				checkflag = true;
				df = -penetration*_normal*_kn;
				force += df;
				moment += ptcm.cross(df);
				
				stress(0) -= df(0)*ptcm(0);
				stress(1) -= df(1)*ptcm(1);
				stress(2) -= df(2)*ptcm(2);
				stress(3) -= 0.5*(df(1)*ptcm(2) + df(2)*ptcm(1));
				stress(4) -= 0.5*(df(2)*ptcm(0) + df(0)*ptcm(2));
				stress(5) -= 0.5*(df(1)*ptcm(0) + df(0)*ptcm(1));
				
				// friction stuff
				v = grain.getVelocity() + grain.getOmegaGlobal().cross(ptcm) - _velocity; // eq (6)
				ds = (v - v.dot(_normal)*_normal)*dt; // eq (7) and (9)
				
				grain.getNodeShearsNonConst()[ptidx] -= ds*_ks;
				grain.getNodeContactNonConst()[ptidx] = _id;
				Fsmag = min(df.norm()*_mu, grain.getNodeShearsNonConst()[ptidx].norm() );
				
				if (Fsmag > 0) {
					Fs = Fsmag*grain.getNodeShearsNonConst()[ptidx]/grain.getNodeShearsNonConst()[ptidx].norm();
					grain.getNodeShearsNonConst()[ptidx] = Fs;
					force += Fs;
//					moment += ptcm.cross(Fs);
				}
			}
			else if (grain.getNodeContactNonConst()[ptidx] == _id ){
				grain.getNodeContactNonConst()[ptidx] = INT_MAX;
				grain.getNodeShearsNonConst()[ptidx] << 0,0,0;
			}
		}
		return checkflag;
	}
	
	
	bool isContact(const Grain3d & grain) const {
		double 	penetration; 			// penetration amount
		for (size_t ptidx = 0; ptidx < grain.getPointList().size(); ptidx++) {
			penetration = _normal.dot(grain.getPointList()[ptidx]) + _d;
			if ( penetration < 0 ) {
				return true;
			}
		}
		return false;
	}
	
	
	// methods to move/rotate the wall
	// moves wall by amount
	void moveWall(const Vector3d & amount) {
		_position += amount;
		_d = -(_normal.dot(_position) );
	}
	// rotates wall about _position
	void rotateWall(const Matrix3d & R) {
		_normal = R*_normal;
		_d = -(_normal.dot(_position) );
	}
	
	// newton-raphson timestep in the normal direction
	void takeTimestep(const double & netForce, const size_t & ncontacts) {
		const double _alpha = 0.7;
		if (ncontacts < 1)
			moveWall( _normal*_alpha*netForce/_kn/1. );
		else
			moveWall( _normal*_alpha*netForce/_kn/double(ncontacts) );
	}
	
	// inertia-based timestep returns the displacement of the wall
	Vector3d takeTimestepInertia(const Vector3d & force, const double & mass, const double & gDamping, const double & dt) {
//		Vector3d tanForce = force - force*_normal.dot(force);
		_velocity = 1./(1.+gDamping*dt/2.)*( (1.-gDamping*dt/2.)*_velocity + dt*force/mass   );
		_position += _velocity*dt;
		return _velocity*dt;
	}
	
	// inertia-based timestep rotationally
//	void takeTimestepRotation(const Vector3d & moment, const double & moin, const double & gDamping, const double & dt) {
//		 
//	}
	 
	void changeKn(const double & newkn) {
		_kn = newkn;
	}
	void changeKs(const double & newks) {
		_ks = newks;
	}
	
	void changeVelocity(const Vector3d & velocity) {
		_velocity = velocity;
	}
	
	// get methods
	const Vector3d & getPosition() const {
		return _position;
	}
	const Vector3d & getNormal() const {
		return _normal;
	}
	const Vector3d & getVelocity() const {
		return _velocity;
	}
	
	
private:
	Vector3d _normal;		// outward normal of wall (normalized to magnitude 1 in constructor)
	Vector3d _position;	// 'center' of wall
	double 	_kn;			// wall stiffness
	double	_ks;			// shear stiffness
	double	_mu;
	Vector3d _velocity;	// velocity in the normal direction
	double 	_d;			// such that the plane can be written in form ax + by + cz + d = 0 where (a,b,c) = _normal
	size_t 	_id;

};


# endif
