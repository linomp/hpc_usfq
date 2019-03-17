/*
 * WallBalls.h
 *
 *  Created on: Nov 17, 2015
 *      Author: reid
 */

#ifndef WALLBALLS_H_
#define WALLBALLS_H_


#include "definitions.h"
#include "WorldStates.h"
#include "Grain3d.h"


class WallBalls { 
	
public:
	// default constructor
	WallBalls() {
		_kn = 0; _gdamping = 0; _exists = false; _nvert = 0; _nangular = 0; _pressure = 0; _height = 0; _l0=0; _kms = 0; _ballRadius=0; _ballMass=0; _dt=0; _segHeight=0; _ballMoin=0; _kmn=0; _nballs=0;
	}
	// constructor if you actually want a working wallBalls
	WallBalls(const double & height, const double & initRadius, const double & ballRadius, const double & pressure, const double & density, const double & kn, const double & kmn, const double & kms, const double & dt, const double & gdamping): 
		_pressure(pressure), _kn(kn), _kmn(kmn), _kms(kms), _dt(dt), _gdamping(gdamping) {
		_exists = true;
		// how to arithmetic?
		double n = M_PI*(initRadius+ballRadius)/ballRadius;
		_nangular = floor(n);
		double a = 1-cos(2.*M_PI/double(_nangular));
		_ballRadius = initRadius*(a+sqrt(2*a))/(2.-a);
		_ballMass = 4./3.*M_PI*_ballRadius*_ballRadius*_ballRadius*density;
		_ballMoin = 0.4*_ballMass*_ballRadius*_ballRadius;
		_l0 = 2*_ballRadius;
		double htheta = M_PI/double(_nangular);
		// note: for some reason this doesn't quite lead to all springs being the same starting length fixme
		_segHeight = sqrt(4.*_ballRadius*_ballRadius-(initRadius+_ballRadius)*(initRadius+_ballRadius)*((1.-cos(htheta))*(1.-cos(htheta))+sin(htheta)*sin(htheta)));
		_nvert = ceil(height/_segHeight);
		_nballs = _nvert*_nangular;
		// set up vectors and zero out vectors that need to be zeroed
		_velocities.resize(_nballs);
		_omegas.resize(_nballs);
		_omegaGlobals.resize(_nballs);
		_dOmegas.resize(_nballs);
		_quats.resize(_nballs);
		_rotMatrices.resize(_nballs);
		_shears.resize(_nballs);
		_normals.resize(_nballs);
		for (size_t i = 0; i < _velocities.size(); i++) {
			_velocities[i] << 0., 0., 0.;
			_omegas[i] << 0, 0, 0;
			_omegaGlobals[i] << 0,0,0;
			_dOmegas[i] << 0,0,0;
			_quats[i] << 0, 0, 0, 0;
			_rotMatrices[i] = Matrix3d::Zero();
			for (size_t j = 0; j < 6; j++) {
				_shears[i][j] << 0., 0., 0.;
				_normals[i][j] << 0., 0., 0.;
			}
		} 
		_positions.resize(_nballs);
		_ballForces.resize(_nballs);
		_ballMoments.resize(_nballs);
		// fill ball vectors
		for (size_t t = 0; t < _nangular; t++) {
			double curtheta = (double(t)/double(_nangular))*2.*M_PI;
			double x = (initRadius+_ballRadius)*cos(curtheta);
			double y = (initRadius+_ballRadius)*sin(curtheta);
			for (size_t h = 0; h < _nvert; h+=2) {
				size_t wallidx = t*_nvert + h;
				_positions[wallidx] << x, y, _ballRadius+_segHeight*double(h);
			}
			// offset theta by a half-angle for odd-numbered layers
			curtheta = ((double(t)+.5)/double(_nangular))*2.*M_PI;
			x = (initRadius+_ballRadius)*cos(curtheta);
			y = (initRadius+_ballRadius)*sin(curtheta);
			for (size_t h = 1; h < _nvert; h+=2) {
				size_t wallidx = t*_nvert + h;
				_positions[wallidx] << x, y, _ballRadius+_segHeight*double(h);
			}
		}
		_height = (double)(_nvert-1)*_segHeight+2*_ballRadius;
		moveHeight(height-_height);
		
		// fill faces vector
		_faces.resize( (_nvert-1)*2*_nangular );
		int upsidedown = 0; // are the faces upside-down?
		int t = 0; int h = 0;
		for (size_t f = 0; f < _faces.size(); f++) {
			if (upsidedown == 0 && h%2 == 0) {
				_faces[f] << t*_nvert+h, ((t+1)%_nangular)*_nvert+h, t*_nvert+h+1;
			}
			else if (upsidedown == 1 && h%2 == 1){
				_faces[f] << t*_nvert+h, ((t+1)%_nangular)*_nvert+h-1, ((t+1)%_nangular)*_nvert+h;
			}
			else if (upsidedown == 0 && h%2 == 1) {
				_faces[f] << ((t+0)%_nangular)*_nvert+h, ((t+1)%_nangular)*_nvert+h, ((t+1)%_nangular)*_nvert+h+1;
			}
			else { // upsidedown == 1 && h%2 == 0
				_faces[f] << t*_nvert+h,  ((t+0)%_nangular)*_nvert+h-1, ((t+1)%_nangular)*_nvert+h;
			}
			t = (t+1)%_nangular;
			if (t == 0) {
				if	(upsidedown == 0) {
					h++;
				}
				upsidedown = 1-upsidedown;
			}
		}
		// fill neighbors vector
		_neighbors.resize(_nballs);
		for (size_t t = 0; t < _nangular; t++) {
			for (size_t h = 0; h < _nvert; h++) {
				size_t n = t*_nvert+h;
				size_t nr = ((t+1)%_nangular)*_nvert+h;
				size_t nur = (n+1 + (h%2)*_nvert)%(_nvert*_nangular);
				size_t nul = (n+1 - ((h+1)%2)*_nvert + _nvert*_nangular)%(_nvert*_nangular);
				size_t nl = ((t+_nangular-1)%_nangular)*_nvert+h;
				size_t nbl = (n-1 - ((h+1)%2)*_nvert + _nvert*_nangular)%(_nvert*_nangular);
				size_t nbr = (n-1+(h%2)*_nvert)%(_nvert*_nangular);
				if (h==0)
					_neighbors[n] << nr, nur, nul, nl, -1, -1;
				else if (h==_nvert-1)
					_neighbors[n] << nr, nl, nbl, nbr, -1, -1;
				else
					_neighbors[n] << nr, nur, nul, nl, nbl, nbr;
			}
		}
		
	} // end constructor
	
	vector<Vector2i> bCircleContact(const Grain3d & grain) const {
		vector<Vector2i> possWalls;
//		Vector3d pos = grain.getPosition();
//		double radius = grain.getRadius();
		for (size_t t = 0; t < _nangular; t++) {
			for (size_t h = 0; h < _nvert; h++) {
				size_t wallidx = t*_nvert+h;
				if ( (grain.getPosition()-_positions[wallidx]).squaredNorm() <
				     (grain.getRadius()+_ballRadius)*(grain.getRadius()+_ballRadius) ) {
						possWalls.push_back(Vector2i(t,h));
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
		Vector3d ptcm;
//		size_t ncontacts = 0;
		for (size_t i = 0; i < possWalls.size(); i++) {
			size_t t = possWalls[i](0);
			size_t h = possWalls[i](1);
//			cout << "second check of grain " << grain.getId() << " possibly contacting with wall " << t << ", " << h << endl;
			size_t wallidx = t*_nvert+h;
//			Vector3d radial(_radii[wallidx](0), _radii[wallidx](1), 0);
//			radial = radial/radial.norm(); // unit length vector from (0,0,h) to 
			for (size_t ptidx = 0; ptidx < grain.getPointList().size(); ptidx++) {
				Vector3d point = grain.getPointList()[ptidx];
				penetration = (point-_positions[wallidx]).norm() - _ballRadius;
				if ( penetration < 0) {
//					ncontacts++;
//					cout << "penetration exists" << endl;
					// outward of the ball
					Vector3d normal = (point-_positions[wallidx])/(point-_positions[wallidx]).norm();
					ptcm = grain.getPointList()[ptidx] - grain.getPosition();
					ptcm = ptcm + penetration*ptcm/ptcm.norm();
					checkflag = true;
					df = -penetration*normal*_kn;
					force += df;
					moment += ptcm.cross(df);
					stress(0) -= df(0)*ptcm(0);
					stress(1) -= df(1)*ptcm(1);
					stress(2) -= df(2)*ptcm(2);
					stress(3) -= 0.5*(df(1)*ptcm(2) + df(2)*ptcm(1));
					stress(4) -= 0.5*(df(2)*ptcm(0) + df(0)*ptcm(2));
					stress(5) -= 0.5*(df(1)*ptcm(0) + df(0)*ptcm(1));
					// update wall's info (make sure we're not trying to update the same wallState's members at the same time)
					#pragma omp critical
					{
						wallState._ballWallContacts[wallidx]++;
						wallState._ballWallForces[wallidx] += penetration*normal*_kn; 
					}
				}
			}
		}
		return checkflag;
	}
	
	void takeTimestep(const WallState3d & wallState) { //const WallState3d & wallState
		clearBallForcesMoments();
		// go through wallState's wallForces and distribute to each ball
		for (size_t b = 0; b < _nballs; b++) {
			_ballForces[b] += wallState._ballWallForces[b];
		}
		// go through neighbors and distribute forces from normal springs and shear springs
		Vector3d force(0,0,0), distvec(0,0,0); double dist = 0;
		#pragma omp parallel for schedule(static,1) firstprivate(force,distvec,dist)
		for (size_t b = 0; b < _nballs; b++) {
			for (size_t n = 0; n < 6; n++) {
				int bn = _neighbors[b](n);
				if (bn < 0) 
					break;
				// so we don't double-count the spring force when b and bn are switched and bn > b (which will eventually happen if it didn't already)
				// parallelization edit: to not have memory races, don't update bn and don't enforce bn < b (but bn will be updated at some other point since we don't enforce bn < b anymore)
//				if (bn < (int)b) {
					// normal contribution
					distvec = _positions[bn]-_positions[b]; // pointing from b to bn (outward normal of b)
					dist = distvec.norm();
					force = _kmn*(dist-_l0)*distvec/dist;
					_ballForces[b] += force;
//					_ballForces[bn] -= force;
					// shear contribution TODO
					Vector3d normal = distvec/dist;
					Vector3d ds = _dt*( _velocities[b] - _velocities[bn] + _omegaGlobals[b].cross(normal*_ballRadius) - _omegaGlobals[bn].cross(-normal*_ballRadius) ); // eq (6)
					if (_normals[b][n].norm() > 0) {
						Vector3d k = normal.cross(_normals[b][n]);
						double sint = k.norm(); // compute the sin and cos of the rotation
						double cost = sqrt(1-sint*sint);
						if (sint > 0) {
							k = k/sint; // normalize k by its norm
						}
						_shears[b][n] = _shears[b][n]*cost + k.cross(_shears[b][n])*sint + k*k.dot(_shears[b][n])*(1.-cost); // rodriguez rotation
					}
					_normals[b][n] = normal;
					_shears[b][n] -= ds*_kms; // eq (8) and (12)
					_ballForces[b] += _shears[b][n];
//					_ballForces[bn] -= _shears[b][n];
					_ballMoments[b] += _ballRadius*normal.cross(_shears[b][n]); // eq (16) 
//					_ballMoments[bn] += _ballRadius*normal.cross(-_shears[b][n]);; // eq (17) 
//				}
			}
		}
		
		// go through faces and distribute pressures
		for (size_t f = 0; f < _faces.size(); f++) {
			// compute area and normal of face
			Vector3i face = _faces[f]; 
			Vector3d crossp = (_positions[face(2)]-_positions[face(0)]).cross(_positions[face(1)]-_positions[face(0)]);
//			double area = crossp.norm()/2.;
//			Vector3d normal = crossp/(area*2); // pointing out of the wall (into the center of the specimen)
			Vector3d ballForce = crossp*_pressure/6.;
			for (size_t b = 0; b < 3; b++) { 
				_ballForces[face(b)] += ballForce;
			}
		}
		
		// make top/bottom layers unable to move at all (because they wrap around the platen in the experiment)
		for (size_t t = 0; t < _nangular; t++) {
			// top layer
			size_t b = t*_nvert + _nvert-1;
//			_ballForces[b](2) = 0;
			_ballForces[b] << 0,0,0;
			_velocities[b] << 0,0,0;
			b = t*_nvert;
//			_ballForces[b](2) = 0;
			_ballForces[b] << 0,0,0;
			_velocities[b] << 0,0,0;
		}
			
		// time integration
		#pragma omp parallel for schedule(static,1)
		for (size_t b = 0; b < _nballs; b++) {
			// center of mass
			_velocities[b] = 1./(1.+_gdamping*_dt/2.)*( (1.-_gdamping*_dt/2.)*_velocities[b] + _dt*_ballForces[b]/_ballMass );
			_positions[b] += _dt*_velocities[b];
			// rotations
			Vector3d principleMoment = _rotMatrices[b].transpose()*_ballMoments[b];
			Vector3d omegaN = _omegas[b] + _dOmegas[b]/2;
			for (int i = 0; i < 5; i++) {
				_dOmegas[b](0) = (principleMoment(0)  - _gdamping*_ballMoin*omegaN(0) )*_dt/_ballMoin;
				_dOmegas[b](1) = (principleMoment(1)  - _gdamping*_ballMoin*omegaN(1) )*_dt/_ballMoin;
				_dOmegas[b](2) = (principleMoment(2)  - _gdamping*_ballMoin*omegaN(2) )*_dt/_ballMoin;
				omegaN = _omegas[b] + _dOmegas[b]/2.;
			}
			_omegas[b] += _dOmegas[b];
			_omegaGlobals[b] = _rotMatrices[b]*_omegas[b];
			// don't need to compute/update rotations since shear increments are computed only by integrating angular velocity
//			_quats[b] << -_quats[b](0), -_quats[b](1), -_quats[b](2), _quats[b](3) ;
//			double Bx = _dt/4.*_omegas[b](0); 
//			double By = _dt/4.*_omegas[b](1); 
//			double Bz = _dt/4.*_omegas[b](2);
//			double c1 =  _quats[b](0) + Bz*_quats[b](1) - Bx*_quats[b](2) - By*_quats[b](3);
//			double c2 = -Bz*_quats[b](0) + _quats[b](1) - By*_quats[b](2) + Bx*_quats[b](3);
//			double c3 =  Bx*_quats[b](0) + By*_quats[b](1) + _quats[b](2) + Bz*_quats[b](3);
//			double c4 =  By*_quats[b](0) - Bx*_quats[b](1) - Bz*_quats[b](2) + _quats[b](3);
//			double detB = 1 + 2*Bx*Bx + 2*By*By + 2*Bz*Bz + 2*Bx*Bx*By*By + 2*By*By*Bz*Bz + 2*Bx*Bx*Bz*Bz + 
//							  Bx*Bx*Bx*Bx + By*By*By*By + Bz*Bz*Bz*Bz;
//			double bfac = (1 + Bx*Bx + By*By + Bz*Bz)/detB;
//			_quats[b](0) = ( c1 + c2*Bz - c3*Bx - c4*By)*bfac;
//			_quats[b](1) = (-c1*Bz + c2 - c3*By + c4*Bx)*bfac;
//			_quats[b](2) = ( c1*Bx + c2*By + c3 + c4*Bz)*bfac;
//			_quats[b](3) = ( c1*By - c2*Bx - c3*Bz + c4)*bfac;
//			// normalize _quat and use it to update rotation matrix
//			_quats[b] = _quats[b]/_quats[b].norm();
//			// convert quaternion back (going from reference to global, which is what it usually is)
//			_quats[b] << -_quats[b](0), -_quats[b](1), -_quats[b](2), _quats[b](3) ;
//			_rotMatrices[b] = qToR(_quats[b]);
		}
	}
	
	// finds the force at which the membrane pulls the top platen
	Vector3d computeInternalForceOnTopLayer() const {
		Vector3d topLayerForce(0,0,0);
		Vector3d distvec, shear;
		double dist;
		for (size_t i = 0; i < _nangular; i++) {
			size_t b = _nballs - i*_nvert - 1;
			for (size_t n = 0; n < 6; n++) {
				int bn = _neighbors[b](n);
				if (bn < 0) 
					break;
				// so we don't double-count the spring force when b and bn are switched and bn > b (which will eventually happen if it didn't already)
				// parallelization edit: to not have memory races, don't update bn and don't enforce bn < b (but bn will be updated at some other point since we don't enforce bn < b anymore)
		//				if (bn < (int)b) {
					// normal contribution
					distvec = _positions[bn]-_positions[b]; // pointing from b to bn (outward normal of b)
					dist = distvec.norm();
					topLayerForce += _kmn*(dist-_l0)*distvec/dist;
//					cout << _l0 << " " << distvec.norm() << endl;
		//					_ballForces[bn] -= force;
					// shear contribution TODO
					Vector3d normal = distvec/dist;
					Vector3d ds = _dt*( _velocities[b] - _velocities[bn] + _omegaGlobals[b].cross(normal*_ballRadius) - _omegaGlobals[bn].cross(-normal*_ballRadius) ); // eq (6)
					shear = _shears[b][n] - ds*_kms; // eq (8) and (12)
					topLayerForce += shear;
			}
		}
		return topLayerForce;
	}
	
	// moves balls in top layer given by rotation R about pos
	void rotateTopLayer(const Matrix3d & R, const Vector3d & pos, const double & dt) {
		for (size_t t = 0; t < _nangular; t++) {
			// top layer
			size_t b = t*_nvert + _nvert-1;
			// because the position of the top ball doesn't exactly coincide with the position of the top wall
			Vector3d trans = _positions[b] - pos;
//			double z = _positions[b](2);
			_positions[b] = R*trans + pos;
//			cout << _positions[b].transpose() << endl;
//			_positions[b](2) += z;
			// update velocity of ball (notably, this will be integrated when shear increments are computed but will be overwritten to 0 prior to time integration)
//			_velocities[b] = R*trans/dt;
		}
	}
	
	// compresses wall vertically and moves balls correspondingly
	void moveHeight(const double & amount) {
		for (size_t i = 0; i < _nangular*_nvert; i++) {
			_positions[i](2) += amount*_positions[i](2)/_height;
		}
		_height += amount; // handling this as a displacement may affect the shear springs (which are calculated from velocity), need to check
		_segHeight += amount/(double)_nvert;
	}
	
	void adjustHeights() {
		for (size_t h = 0; h < _nvert; h++) {
			double ringHeight = _segHeight*double(h) + _ballRadius;
			for (size_t t = 0; t < _nangular; t++) {
				size_t wallidx = t*_nvert+h;
				_positions[wallidx] << _positions[wallidx](0), _positions[wallidx](1), ringHeight;
			}
		}
	}
	
	void changeKn(const double & newkn) {
		_kn = newkn;
	}
	
	void changePressure(const double & pressure) {
		_pressure = pressure;
	}
	
	void moveBall(const size_t & idx, const Vector3d & pos) {
		_positions[idx] = pos;
	}
	
	// changes ONLY the height (use this)
	void changeHeightOnly(const double & height) {
		_height = height;
	}
	void moveHeightOnly(const double & amt) {
		_height += amt;
	}

	double findTopArea() const {
		double topArea = 0;
		for (size_t t = 0; t < _nangular; t++) {
			size_t idx0 = t*_nvert + _nvert-2;
			size_t idx1 = ((t+1)%_nangular)*_nvert + _nvert-2;
			topArea += _positions[idx0](0)*_positions[idx1](1)-_positions[idx0](1)*_positions[idx1](0);
		}
		return topArea/2.;
	}
	
	double findVolume() const {
		double volume = 0;
		for (size_t h = 0; h < _nvert; h++) {
			for (size_t t = 0; t < _nangular; t++) {
				size_t idx0 = t*_nvert+h;
				size_t idx1 = ((t+1)%_nangular)*_nvert+h;
				volume += _positions[idx0](0)*_positions[idx1](1)-_positions[idx0](1)*_positions[idx1](0); // double the area of a slice
			}
		}
		return volume/double(_nvert)/2.*_height;
	}
	
	void changel0(const double & l0) {
		_l0 = l0;
	}
	
	// get methods
	const bool & exists() const {
		return _exists;
	}
	const double & getBallRadius() const {
		return _ballRadius;
	}
	const double & getKmn() const {
		return _kmn;
	}
	const double & getKms() const {
		return _kms;
	}
	const double & getL0() const {
		return _l0;
	}
	const size_t & getNumWalls() const {
		return _nballs;
	}
	const size_t & getNangular() const {
		return _nangular;
	}
	const size_t & getNvert() {
		return _nvert;
	}
	const Vector3d & getPosition(const size_t & idx) const {
		return _positions[idx];
	}
	const vector<Vector3d> & getPositions() const {
		return _positions;
	}
	const vector<Vector3i> & getFaces() const {
		return _faces;
	}
	const vector<Vector6i> & getNeighbors() const {
		return _neighbors;
	}
	const double & getPressure() const{
                return _pressure;
        }
	
private:
	
	void clearBallForcesMoments() {
		for (size_t i = 0; i < _ballForces.size(); i++) {
			_ballForces[i] << 0., 0., 0.;
			_ballMoments[i] << 0., 0., 0.;
		}
	}
	
	Matrix3d qToR(const Vector4d & quat) {
		Matrix3d rotMatrix;
		rotMatrix(0,0) = -quat(0)*quat(0) + quat(1)*quat(1) - quat(2)*quat(2) + quat(3)*quat(3);
		rotMatrix(0,1) = -2*(quat(0)*quat(1) - quat(2)*quat(3));
		rotMatrix(0,2) =  2*(quat(1)*quat(2) + quat(0)*quat(3));
		rotMatrix(1,0) = -2*(quat(0)*quat(1) + quat(2)*quat(3));
		rotMatrix(1,1) =  quat(0)*quat(0) - quat(1)*quat(1) - quat(2)*quat(2) + quat(3)*quat(3);
		rotMatrix(1,2) = -2*(quat(0)*quat(2) - quat(1)*quat(3));
		rotMatrix(2,0) =  2*(quat(1)*quat(2) - quat(0)*quat(3));
		rotMatrix(2,1) = -2*(quat(0)*quat(2) + quat(1)*quat(3));
		rotMatrix(2,2) = -quat(0)*quat(0) - quat(1)*quat(1) + quat(2)*quat(2) + quat(3)*quat(3);
		return rotMatrix;
	}
	
	// must be updated together
	double _height;		// total height of wall (bottom is always at z=0 deal w/ it)
	double _segHeight;	// height of segment (starts at ball height but decreases during compression);
	
	// must be updated together
	vector<Vector3d> _velocities; // velocity of each ball
	vector<Vector3d> _positions;	// (x,y,z) position of each ball
	vector<Vector3d> _omegas;		// angular velocities of each ball
	vector<Vector3d> _omegaGlobals;
	vector<Vector3d> _dOmegas;
	vector<Vector4d> _quats;		// quaternions of each ball
	vector<Matrix3d> _rotMatrices;
	vector<array<Vector3d, 6> > _shears; //
	vector<array<Vector3d, 6> > _normals; 
	
	// mesh-like stuff
	vector<Vector3d> _ballForces; // list of forces on the balls
	vector<Vector3d> _ballMoments; // list of moments on balls
	vector<Vector3i> _faces;		// list of faces (each face has the index of 3 balls that make up its vertices)
	vector<Vector6i> _neighbors;	// list of neighbors of each ball
	double _l0;							// initial spring length
	
	// these values are never changed
	bool _exists;			// do we even have a WallBalls?
	size_t _nvert;			// number of balls in the vertical direction
	size_t _nangular;		// number of balls around
	size_t _nballs;		// total number of balls
	double _pressure;		// radial pressure (applied to each face)
	double _kn;				// ball normal stiffness
	double _kmn;			// mesh normal stiffness
	double _kms;			// mesh shear stiffness
	double _dt;				// timestep interval
	double _gdamping;		// global damping
	double _ballRadius;	// radius of each ball
	double _ballMass;		// mass of each ball
	double _ballMoin;		// moment of inertia of each ball
	


	
};



#endif /* WALLBALLS_H_ */
