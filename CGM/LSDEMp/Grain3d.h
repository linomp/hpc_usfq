/*
 * Grain3d.h
 *
 *  Created on: July 15, 2014
 *      Author: Reid
 */

#ifndef GRAIN3D_H_
#define GRAIN3D_H_

#include "definitions.h"
#include "Levelset3d.h"


class Grain3d {
public:
	// constructors
	Grain3d() {
		_radius = 0; _rsq = 0; _kn = 0; _ks = 0; _mass = 0; _mu = 0; _id = 0; _density = 1.0; _lset = NULL; 
		_quat << 0,0,0,0; _omega << 0,0,0; _momentInertia << 0,0,0; 
		_position << 0,0,0; _velocity << 0,0,0; 
	}

	Grain3d(const double & mass, const Vector3d & position, const Vector3d & velocity,
			  const Vector3d & momentInertia, const Vector4d & quat, const Vector3d & omega,
			  const Vector3d & cmLset, const vector<Vector3d> & refPointList, const double & radius, 
			  Levelset3d * lset, const double kn, const double ks, const double mu, const size_t & id):
			  _mass(mass), _position(position), _velocity(velocity),  _momentInertia(momentInertia),
			  _quat(quat), _omega(omega), _cmLset(cmLset), _refPointList(refPointList), _radius(radius), 
			  _lset(lset), _kn(kn), _ks(ks), _mu(mu), _id(id) {
		_density = 1.0;
		_rsq = _radius*_radius;
		// clean up grain because the characterization method sucks (sometimes)
		for (size_t p = 0; p < _refPointList.size(); p++) {
			if (_refPointList[p].norm() > 1e3 || isnan(_refPointList[p](0))) {
				_refPointList[p] = _refPointList.back();
				_refPointList.pop_back();
				p--;
//				_refPointList[p] << 0,0,0;
			}
		}
		// recheck radius
		double radius2 = 0;
		for (size_t q = 0; q < _refPointList.size(); q++) {
			if (_refPointList[q].norm() > radius2) {
				radius2 = _refPointList[q].norm();
			}
		}
		_radius = radius2;
		
		_dOmega << 0., 0., 0.;
		updateRotationMatrix(quat);
		_omegaGlobal = _rotMatrix*_omega;
		_pointList.resize( _refPointList.size());
		_nodeShears.resize(_pointList.size()) ;
		_nodeContact.resize(_pointList.size());
		_nodeNormals.resize(_pointList.size());
		// calculate _pointList based on _refPointList and the rotation matrix and the center of mass
		// also zero out the node information stuff
		for (size_t i = 0; i < _refPointList.size(); i++) {
			_nodeNormals[i] << 0,0,0;
			_pointList[i] = _rotMatrix*_refPointList[i] + _position;
			_nodeShears[i] << 0,0,0;
			_nodeContact[i] = 0;
		}
	}
 
	bool bCircleCheck(const Grain3d & other) const {
		Vector3d d = other.getPosition() - _position;
		double rsum = other.getRadius() + _radius;		
		return d.squaredNorm() < rsum*rsum;
	}
	
	// Checks contact between *this and other.  If there is no contact, returns false.
	// Compares points of *this to the level set of other.
	// If there is contact, returns true and updates force, which is the force on *this,
	// thisMoment, which is the moment on *this, and otherMoment, the moment on other.
	bool findInterparticleForceMoment(const Grain3d & other, const double & dt, Vector3d & force, Vector3d & thisMoment, Vector3d & otherMoment) {
		// Zero the outputs force, thisMoment, and otherMoment
		force << 0., 0., 0.;
		thisMoment << 0., 0., 0.;
		otherMoment << 0., 0., 0.;
		
		bool checkflag = false; // return flag (assigned to true if contact exists)
		Vector3d ptThisCM; 		// point wrt the center of mass of *this in real space
		Vector3d ptOtherCM; 		// point wrt the center of mass of other in real space
		Vector3d ptOtherLset; 	// point in the reference config of other's level set
		double	penetration;	// penetration amount
		Vector3d normal; 			// surface normal pointing out of other in the reference config (initially) and then pointing out of *this in real space (after)
		Vector3d Fn; 				// normal force from a single point
		Vector3d v; 				// velocity of *this relative to other at a point of contact
		Vector3d Fs; 				// vector of frictional/shear force
		Vector3d ds;				// tangential increment (projection of relative velocity v in tangengial direction)
		double Fsmag;
		Vector3d k; 				// axis of rotation to rotate old normal to new normal
		double sint, cost;
		size_t ncontacts = 0;	// number of contact points
		// iterate through all of the points of *this and check for contact for each one
		for (size_t ptidx = 0; ptidx < _pointList.size(); ptidx++) {
			ptOtherCM = _pointList[ptidx] - other.getPosition();
			// bounding radius check
			if ( ptOtherCM.squaredNorm()<other.getRsq()){
				ptOtherLset = other.getRotMatrix().transpose()*ptOtherCM;
				ptOtherLset += other.getCmLset();
				// check the point against the level set of other
				if (other.getLset()->findPenetration(ptOtherLset, penetration, normal) ) {
					ptThisCM = _pointList[ptidx] - _position;
					ncontacts++;
					checkflag = true;
					// rotate the normal from the reference config of other's level set to real space
					// and change the direction of the normal to point out of *this to match KWL's 3d GEM paper
					normal = -other.getRotMatrix()*normal;
					// update force: normal force contribution
					// note: penetration is negative which negates the negative sign in eq (2) of KWL
					Fn = penetration*normal*_kn;
					force += Fn;
					// update moments: eccentric loading contribution
					thisMoment  += ptThisCM.cross(Fn); // eq (4)
					otherMoment += ptOtherCM.cross(-Fn); // eq (5)
	//				cout << otherMoment(2) << endl;
					// force/moment calculations based on friction 
					v = _velocity - other.getVelocity() + _omegaGlobal.cross(ptThisCM) - other.getOmegaGlobal().cross(ptOtherCM); // eq (6)
					ds = (v - v.dot(normal)*normal)*dt; // eq (7) and (9)
					
					// rotate the shear force into the new tangent direction
					if (_nodeNormals[ptidx].norm() > 0) {
						k = normal.cross(_nodeNormals[ptidx]); // axis of rotation from n_old to n_cur
						sint = k.norm(); // compute the sin and cos of the rotation
						cost = sqrt(1-sint*sint);
						k = k/(sint+DBL_MIN); // normalize k to unit magnitude, also don't divide by 0
						_nodeShears[ptidx] = _nodeShears[ptidx]*cost + k.cross(_nodeShears[ptidx])*sint + k*k.dot(_nodeShears[ptidx])*(1.-cost);
					}
					_nodeNormals[ptidx] = normal;
					_nodeShears[ptidx] -= ds*_ks; // eq (8) and (12)
					_nodeContact[ptidx] = other.getId();
					
					Fsmag = min(Fn.norm()*_mu, _nodeShears[ptidx].norm() ); // eq (14)
					if (Fsmag > 0) {
						Fs = Fsmag*_nodeShears[ptidx]/_nodeShears[ptidx].norm(); // eq (13)
						_nodeShears[ptidx] = Fs;
						force += Fs;
						thisMoment  += ptThisCM.cross(Fs); // eq (16) 
						otherMoment += ptOtherCM.cross(-Fs); // eq (17) 
					}
				}
				 
			}
			// if there is no contact between the point and other, reset the shear force if other was the last contact
			else if (_nodeContact[ptidx] == other.getId() ){
				_nodeContact[ptidx] = INT_MAX;
				_nodeShears[ptidx] << 0,0,0;
				_nodeNormals[ptidx] << 0,0,0;
			}
//			else if (_conMap.count(ptidx) > 0 ){
//				if (_conMap[ptidx]._otherId == other.getId()) {
//					_conMap.erase(ptidx);
//				}
//			}
		}
		// prevent instability by limiting effective number of contact points based on mass
		double ncmax = 16;
		if ((double)ncontacts > ncmax) {
			force *= ncmax/(double)ncontacts;
			thisMoment *= ncmax/(double)ncontacts;
			otherMoment *= ncmax/(double)ncontacts;
		}
		return checkflag;
	}
	
	// Finds contact data between *this and other.
	CData findContactData(const Grain3d & other) const {
		CData cData;				// output
		Vector3d force;			// total force
		Vector3d ptThisCM; 		// point wrt the center of mass of *this in real space
		Vector3d ptOtherCM; 		// point wrt the center of mass of other in real space
		Vector3d ptOtherLset; 	// point in the reference config of other's level set
		double	penetration;	// penetration amount
		Vector3d normal; 			// surface normal pointing out of other in the reference config (initially) and then pointing out of *this in real space (after)
		double Fsmag;

		// iterate through all of the points/nodes of *this and check for contact for each one
		for (size_t ptidx = 0; ptidx < _pointList.size(); ptidx++) {
			ptThisCM = _pointList[ptidx] - _position;
			ptOtherCM = _pointList[ptidx] - other.getPosition();
			ptOtherLset = other.getRotMatrix().transpose()*ptOtherCM;
			ptOtherLset += other.getCmLset();
			if ( other.getLset()->findPenetration(ptOtherLset, penetration, normal) ) {
				normal = -other.getRotMatrix()*normal;
				force = penetration*normal*_kn;
				Fsmag = min(force.norm()*_mu, _nodeShears[ptidx].norm() ); // eq (14)
				if (Fsmag > 0) {
					force += Fsmag*_nodeShears[ptidx]/_nodeShears[ptidx].norm(); // eq (13)
				}
				cData._cpairs.push_back(Vector2i(_id, other.getId()));
				cData._nodes.push_back(ptidx);
				cData._forces.push_back(force);
				cData._normals.push_back(normal);
				cData._clocs.push_back(_pointList[ptidx]);
			}
		}
		return cData;
	}
	 
	// takes a timestep which affects all related quantities EXCEPT updating the points
	void takeTimestep(const Vector3d & force, const Vector3d & moment, const double & gDamping, const double & dt) {
		// update velocity and position of the center of mass
		_velocity = 1./(1.+gDamping*dt/2.)*( (1.-gDamping*dt/2.)*_velocity + dt*force/_mass   );
		// update angle/rotation of the grain
		// first rotate moment into the principle frame
		Vector3d principleMoment = _rotMatrix.transpose()*moment;
		// predictor-corrector iteration for angular velocities (in principle frame, 5 iterations) 
		Vector3d omegaN;
		omegaN = _omega + _dOmega/2;
//		cout << omegaN.transpose() << endl;
		for (int i = 0; i < 5; i++) {
			_dOmega(0) = (principleMoment(0) + omegaN(1)*omegaN(2)*(_momentInertia(1)-_momentInertia(2)) - gDamping*_momentInertia(0)*omegaN(0) )*dt/_momentInertia(0);
			_dOmega(1) = (principleMoment(1) + omegaN(2)*omegaN(0)*(_momentInertia(2)-_momentInertia(0)) - gDamping*_momentInertia(1)*omegaN(1) )*dt/_momentInertia(1);
			_dOmega(2) = (principleMoment(2) + omegaN(0)*omegaN(1)*(_momentInertia(0)-_momentInertia(1)) - gDamping*_momentInertia(2)*omegaN(2) )*dt/_momentInertia(2);
			omegaN = _omega + _dOmega/2.;
		}
		_omega += _dOmega;
		// rotate angular velocities back to global frame
		_omegaGlobal = _rotMatrix*_omega;
		// everything after this should be correct
		// update quaternions and rotation matrix
		// convert quat to go from global to reference (the opposite of what it usually is, but that's needed for this bit here)
		_quat << -_quat(0), -_quat(1), -_quat(2), _quat(3) ;
		double Bx = dt/4.*_omega(0); 
		double By = dt/4.*_omega(1); 
		double Bz = dt/4.*_omega(2);
		double c1 =  _quat(0) + Bz*_quat(1) - Bx*_quat(2) - By*_quat(3);
		double c2 = -Bz*_quat(0) + _quat(1) - By*_quat(2) + Bx*_quat(3);
		double c3 =  Bx*_quat(0) + By*_quat(1) + _quat(2) + Bz*_quat(3);
		double c4 =  By*_quat(0) - Bx*_quat(1) - Bz*_quat(2) + _quat(3);
		double detB = 1 + 2*Bx*Bx + 2*By*By + 2*Bz*Bz + 2*Bx*Bx*By*By + 2*By*By*Bz*Bz + 2*Bx*Bx*Bz*Bz + 
						  Bx*Bx*Bx*Bx + By*By*By*By + Bz*Bz*Bz*Bz;
		double bfac = (1 + Bx*Bx + By*By + Bz*Bz)/detB;
		_quat(0) = ( c1 + c2*Bz - c3*Bx - c4*By)*bfac;
		_quat(1) = (-c1*Bz + c2 - c3*By + c4*Bx)*bfac;
		_quat(2) = ( c1*Bx + c2*By + c3 + c4*Bz)*bfac;
		_quat(3) = ( c1*By - c2*Bx - c3*Bz + c4)*bfac;
		// normalize _quat and use it to update rotation matrix
		_quat = _quat/_quat.norm();
		// convert quaternion back (going from reference to global, which is what it usually is)
		_quat << -_quat(0), -_quat(1), -_quat(2), _quat(3) ;
		updateRotationMatrix(_quat);
		// update center of mass
		_position += dt*_velocity;
	}
	
	// updates positions of points based on _refPointList and the grain's _position and _rotation
	void updatePoints() {
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] = _rotMatrix*_refPointList[ptid] + _position;
		}
	}

	
	// calculation methods
	double computeKineticEnergy() const {
		// KE = .5*m*v^2 + .5*w'*I*w
		double ke = 0;
		// translational component
		ke += .5*_mass*_velocity.squaredNorm();
		// rotational component
		for (int i = 0; i < 3; i++) {
			ke += .5*_omega(i)*_momentInertia(i)*_omega(i);
		}
		return ke;
	}
	
	/* Methods that change various properties of the grain */
	// changes density
	void changeDensity(const double & density) {
		_mass *= density/_density;
		_momentInertia *= density/_density;
		_density = density;
	}
	
	// Methods which move/rotate the grain BY a displacement/rotation 
	// rotates the grain about point pt by rotation R
	void rotateGrain(const Matrix3d & R, const Vector3d & pt) {
		_rotMatrix = R*_rotMatrix;
		updateQuat(_rotMatrix);
		_position = R*(_position-pt) + pt;
		// update points
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] = R*(_pointList[ptid]-pt) + pt;
		}
	}
	// moves grain by amount
	void moveGrain(const Vector3d & amount) {
		_position = _position + amount;
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] += amount;
		}
	}
	// Methods which move/rotate the grain to a PRESCRIBED location/rotation
	// changes the position to a prescribed location loc
	void changeLocation(const Vector3d & loc) {
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
 			_pointList[ptid] += loc-_position;
		}
		_position = loc;
	}
	// changes the rotation to a prescribed quaternion quat
	void changeRotation(const Vector4d & quat) {
		_quat = quat/quat.norm();
		updateRotationMatrix(_quat);
		for (size_t ptid = 0; ptid < _refPointList.size(); ptid++) {
			_pointList[ptid] = _rotMatrix*_refPointList[ptid] + _position;
		}
	}
	// change stiffness
	void changeKn(const double & kn) {
		_kn = kn;
	}
	void changeMu(const double & mu) {
		_mu = mu;
	}
	void changeKs(const double & ks) {
		_ks = ks;
	}
	void scaleMass(const double & scale) {
		_mass *= scale;
	}
	void scaleMI(const double & scale) {
		_momentInertia *= scale;
	}
	void changeShearHist(const size_t & nodeidx, const size_t & contactingidx, const Vector3d & shear, const Vector3d & normal) {
		_nodeShears[nodeidx] = shear;
		_nodeContact[nodeidx] = contactingidx;
		_nodeNormals[nodeidx] = normal;
	}
	// erase friction history
	void clearFriction() {
		for (size_t i = 0; i < _nodeShears.size(); i++) {
			_nodeShears[i] << 0,0,0;
			_nodeNormals[i] << 0,0,0;
			_nodeContact[i] = 0;
		}
	}
	
	// get methods
	const double & getMass() const {
		return _mass;
	}
	const Vector3d & getPosition() const {
		return _position;
	}
	const Vector3d & getVelocity() const {
		return _velocity;
	}
	const Vector4d & getQuat() const {
		return _quat;
	}
	const Matrix3d & getRotMatrix() const {
		return _rotMatrix;
	}
	const Vector3d & getOmega() const {
		return _omega;
	}
	const Vector3d & getOmegaGlobal() const {
		return _omegaGlobal;
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
	const double & getRsq() const {
		return _rsq;
	}
	const Levelset3d * getLset() const {
		return _lset;
	}
	const double & getKn() const {
		return _kn;
	}
	const double & getKs() const {
		return _ks;
	}
	const double & getMu() const {
		return _mu;
	}
	const double & getDensity() const {
		return _density;
	}
	const vector<Vector3d> & getPointList() const {
		return _pointList;
	}
	const size_t & getId() const {
		return _id;
	}
	const vector<size_t> & getNodeContact() const {
		return _nodeContact;
	}
	const vector<Vector3d> & getNodeShears() const {
		return _nodeShears;
	}
	const vector<Vector3d> & getNodeNormals() const {
		return _nodeNormals;
	}
	 
	// non const
	vector<Vector3d> & getNodeShearsNonConst() {
		return _nodeShears;
	}
	vector<size_t> & getNodeContactNonConst() {
		return _nodeContact;
	}
	
private:
	void updateQuat(const Matrix3d & R) {
		double tr = R.trace();
		_quat(3) = sqrt(tr+1.)/2.;
		_quat(0) = (R(0,2)-R(2,0))/(4*_quat(3));
		_quat(1) = (R(1,2)-R(2,1))/(4*_quat(3));
		_quat(2) = (R(0,1)-R(1,0))/(4*_quat(3));
	}
	// Converts a quaternion into a rotation matrix which is stored in _rotMatrix
	void updateRotationMatrix(const Vector4d & quat) {
		_rotMatrix(0,0) = -quat(0)*quat(0) + quat(1)*quat(1) - quat(2)*quat(2) + quat(3)*quat(3);
		_rotMatrix(0,1) = -2*(quat(0)*quat(1) - quat(2)*quat(3));
		_rotMatrix(0,2) =  2*(quat(1)*quat(2) + quat(0)*quat(3));
		_rotMatrix(1,0) = -2*(quat(0)*quat(1) + quat(2)*quat(3));
		_rotMatrix(1,1) =  quat(0)*quat(0) - quat(1)*quat(1) - quat(2)*quat(2) + quat(3)*quat(3);
		_rotMatrix(1,2) = -2*(quat(0)*quat(2) - quat(1)*quat(3));
		_rotMatrix(2,0) =  2*(quat(1)*quat(2) - quat(0)*quat(3));
		_rotMatrix(2,1) = -2*(quat(0)*quat(2) + quat(1)*quat(3));
		_rotMatrix(2,2) = -quat(0)*quat(0) - quat(1)*quat(1) + quat(2)*quat(2) + quat(3)*quat(3);
	}
	double 		_mass;
	Vector3d 	_position; 		// location of center of mass in real space
	Vector3d 	_velocity;		
	Vector3d 	_momentInertia;// moment of inertia in principal frame (purely diagonal terms)
	Vector4d 	_quat; 			// quaternion representing the rotation in (X,Y,Z,W) form
	Matrix3d		_rotMatrix; 	// rotation matrix from principle frame to current frame
	Vector3d		_omega; 			// angular velocity IN THE PRINCIPLE FRAME
	Vector3d		_dOmega; 		// change in angular velocity between timesteps (needed for predictor-corrector algorithm) IN THE PRINCIPLE FRAME
	Vector3d		_cmLset; 		// center of mass wrt the level set reference configuration
	vector<Vector3d> _pointList; 		// list of points comprising the grain in real space (translated and rotated)
	vector<Vector3d> _refPointList; 	// list of points in reference config (center of mass is at (0,0,0) and I is diagonal)
	double 		_radius; // radius of bounding sphere
	double		_rsq;		// squared radius
	Levelset3d *_lset;	// level set of grain
	double		_kn; 		// normal stiffness
	double		_ks; 		// shear stiffness (currently not used)
	double		_mu; 		// interparticle friction
	double 		_density;// density of particle (default = 1)
	Vector3d _omegaGlobal;
	size_t _id;
//	map<size_t, ConData> _conMap; // contact map for node contacts
	vector<Vector3d> _nodeShears; // shear force at each node
	vector<size_t> _nodeContact;  // index of grain the node is contacting
	vector<Vector3d> _nodeNormals;// contact normals at each node
	
};

#endif /* Grain2D_H_ */
