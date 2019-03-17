/*
 * WallShell.h
 *
 *  Created on: Nov 17, 2015
 *      Author: reid
 */

#ifndef WALLSHELL_H_
#define WALLSHELL_H_


#include "definitions.h"
#include "WorldStates.h"
#include "Grain3d.h"


class WallShell {
	
public:
	WallShell() {
		_kn = 0; _mu = 0; _exists = false; _nvert = 0; _nangular = 0;
	}
	WallShell(const double & height, const double & initRadius, const size_t & nvert, const size_t & nangular, const double & pressure, const double & kn, const double & mu): 
		_height(height), _nvert(nvert), _nangular(nangular), _pressure(pressure), _kn(kn), _mu(mu)  {
		_exists = true;
		
		_normals.resize(_nangular*_nvert);
		_areas.resize(_nangular*_nvert);
		_positions.resize(_nangular*_nvert);
		_radii.resize(_nangular*_nvert);
		_dvals.resize(_nangular*_nvert);
		
		_theta = 2.*M_PI/double(_nangular);
		
		// compute area
		double area = 2*initRadius*tan(M_PI/_nangular)*_height/_nvert;
		for (size_t t = 0; t < _nangular; t++) {
			
			double curtheta = double(t)/double(_nangular)*2.*M_PI;
//			cout << curtheta << endl;
			double x = initRadius*cos(curtheta);
			double y = initRadius*sin(curtheta);
			
			for (size_t h = 0; h < _nvert; h++) {
				
				size_t wallidx = t*_nvert + h;
				_normals[wallidx] << -cos(curtheta), -sin(curtheta), 0;
				_dvals[wallidx] = -(_normals[wallidx](0)*x + _normals[wallidx](1)*y);
//				cout << _normals[wallidx].transpose() << endl;
//				_positions[t*_nangular + h] << x, y, _height*h/_nvert;
				_positions[wallidx] << x, y, 0;
				_areas[wallidx] = area;
				_radii[wallidx] = initRadius;
			}
		}
	}
	
	
	vector<Vector2i> bCircleContact(const Grain3d & grain) const {
		vector<Vector2i> possWalls;
		Vector3d pos = grain.getPosition();
		double radius = grain.getRadius();
		for (size_t t = 0; t < _nangular; t++) {
			for (size_t h = 0; h < _nvert; h++) {
				double botHeight = _height*double(h)/double(_nvert);
				double topHeight = _height*double(h+1)/double(_nvert);
				size_t wallidx = t*_nvert+h;
				if (_normals[wallidx].dot(pos) + _dvals[wallidx] < grain.getRadius()) {
					if ( (pos(2) < topHeight && pos(2) > botHeight) ||
						  (pos(2) + radius > botHeight && pos(2) < botHeight) || 
						  (pos(2) - radius < topHeight && pos(2) > topHeight) ) {
//						cout << "grain " << grain.getId() << " possibly contacting with wall " << t << ", " << h << endl;
						possWalls.push_back(Vector2i(t,h));
						
					}
				}
			}
		}
		return possWalls;
	}
	
	// frictionless walls - for now
	bool findWallForceMoment(const Grain3d & grain, const vector<Vector2i> & possWalls, Vector3d & force, Vector3d & moment, Vector6d & stress, WallState3d & wallState) {
		// zero the input vars
		force << 0., 0., 0.;
		moment << 0., 0., 0.;
		stress << 0,0,0,0,0,0;
		bool 		checkflag = false;	// changes to true if penetration exists
		double 	penetration; 			// penetration amount
		Vector3d df;						// force increment on grain
//		Vector3d v;							// relative velocity
//		Vector3d sdot;						// projection of relative velocity into tangential direction
//		Vector3d Fs;
		Vector3d ptcm;
		
		for (size_t i = 0; i < possWalls.size(); i++) {
			size_t t = possWalls[i](0);
			size_t h = possWalls[i](1);
//			cout << "second check of grain " << grain.getId() << " possibly contacting with wall " << t << ", " << h << endl;
			double botHeight = _height*double(h)/double(_nvert);
			double topHeight = _height*double(h+1)/double(_nvert);
			size_t wallidx = t*_nvert+h;
			for (size_t ptidx = 0; ptidx < grain.getPointList().size(); ptidx++) {
				Vector3d point = grain.getPointList()[ptidx];
				penetration = _normals[wallidx].dot(point) + _dvals[wallidx];
				if ( penetration < 0 && point(2) < topHeight && point(2) > botHeight) {
//					cout << "shell wall penetration" << endl;
					ptcm = grain.getPointList()[ptidx] - grain.getPosition();
					ptcm = ptcm + penetration*ptcm/ptcm.norm();
					checkflag = true;
					df = -penetration*_normals[wallidx]*_kn;
					force += df;
					moment -= ptcm.cross(df);
					stress(0) -= df(0)*ptcm(0);
					stress(1) -= df(1)*ptcm(1);
					stress(2) -= df(2)*ptcm(2);
					stress(3) -= 0.5*(df(1)*ptcm(2) + df(2)*ptcm(1));
					stress(4) -= 0.5*(df(2)*ptcm(0) + df(0)*ptcm(2));
					stress(5) -= 0.5*(df(1)*ptcm(0) + df(0)*ptcm(1));
					// update wall's info
					wallState._shellWallContacts[wallidx]++;
					wallState._shellWallForces[wallidx] -= penetration*_kn;
//					_ncontacts[wallidx]++;
//					_forces[wallidx] -= penetration*_kn;
				}
			}
		}
		return checkflag;
	}
	
	// moves wall given by wallidx by amount, updates values
	void moveWall(const size_t & wallidx, const Vector3d & amount) {
		_positions[wallidx] += amount;
		_dvals[wallidx] = -(_normals[wallidx].dot(_positions[wallidx]));
		_radii[wallidx] = _positions[wallidx].norm();
	}
	
	void takeTimestep(const WallState3d & wallState) {
		const double _alpha = 0.7;
		for (size_t t = 0; t < _nangular; t++) {
			// prescribe the top and bottom rings to not move to simulate platen-grain behavior
//			for (size_t h = 1; h < _nvert-1; h++) {
			for (size_t h = 0; h < _nvert; h++) {
				size_t wallidx = t*_nvert+h;
				size_t ncontacts = wallState._shellWallContacts[wallidx];
				double netForce = _pressure*_areas[wallidx] - wallState._shellWallForces[wallidx];
				if (ncontacts < 5) {
					moveWall( wallidx, _alpha*netForce/_kn/5.*_normals[wallidx] );
				}
				else {
					moveWall( wallidx, _alpha*netForce/_kn/double(ncontacts)*_normals[wallidx] );
				}
			}
		}
		// update areas after all walls are moved
		double segHeight = _height/double(_nvert);
		for (size_t t = 0; t < _nangular; t++) {
			size_t tNext = (t+1)%_nangular;
			size_t tPrev = (t+_nangular-1)%_nangular;
			for (size_t h = 0; h < _nvert; h++) {
				size_t idxCurr = t    *_nvert+h;
				size_t idxNext = tNext*_nvert+h;
				size_t idxPrev = tPrev*_nvert+h;
				double aPrev = (_radii[idxPrev] - _radii[idxCurr]*cos(_theta))/sin(_theta);
				double aNext = (_radii[idxNext] - _radii[idxCurr]*cos(_theta))/sin(_theta);
				_areas[idxCurr] = segHeight*(aPrev+aNext);
//				cout << _areas[idxCurr] << endl;
			}
		}
	}
	
	void changeHeight(const double & newHeight) {
		for (size_t i = 0; i < _nangular*_nvert; i++) {
			_areas[i] *= newHeight/_height;
		}
		_height = newHeight;
	}
	
	void changeKn(const double & newkn) {
		_kn = newkn;
	}
	

	
	double findVolume() const {
		double volume = 0;
		for (size_t i = 0; i < _nangular*_nvert; i++) {
			volume += _areas[i]*_radii[i]/2.;
		}
		return volume;
	}
	
	double findTopArea() const {
		double topArea = 0;
		for (size_t i = _nangular*(_nvert-2); i < _nangular*(_nvert-1); i++) {
			topArea += _areas[i]*_radii[i]/2.;
		}
		return topArea/(_height/_nvert);
	}
	
	// get methods
	bool exists() const {
		return _exists;
	}
	size_t getNumWalls() const {
		return _nvert*_nangular;
	}
	
	
private:
	// these values are never changed
	bool _exists;			// does we even have a WallShell?
	double _height;		// total height of wall (bottom is always at z=0 deal w/ it)
	size_t _nvert;			// number of walls in the vertical direction
	size_t _nangular;		// number of walls around
	double _pressure;		// pressure applied to wall shells
	vector<Vector3d> _normals;	// precomputed vector of normals of size _nvert*_nangular
	double	_theta;		// angle between each wall
	double 	_kn;			// wall stiffness
	double	_mu;
	
	// must be updated at each timestep
	vector<double> _areas;
	vector<Vector3d> _positions; // positions of the centers of each wall
	vector<double> _radii;		// distance of each wall from (0,0)
	vector<double> _dvals;		// vector of values of d (such that the plane can be written in form ax + by + cz + d = 0 where (a,b,c) = _normal) of size _nvert*_nangular
	

	
};



#endif /* WALLSHELL_H_ */
