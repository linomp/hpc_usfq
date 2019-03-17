/*
 * WorldStates.h
 *
 *  Created on: Dec 7, 2015
 *      Author: reid
 */

#ifndef WORLDSTATES_H_
#define WORLDSTATES_H_

struct GrainState3d {
	GrainState3d() {}
	GrainState3d(const Vector6d & svtemp, const vector<Vector3d> & gftemp, const vector<Vector3d> & gmtemp):
		_stressVoigt(svtemp), _grainForces(gftemp), _grainMoments(gmtemp) {}
	void resize(const int & ngrains) {
		_grainForces.resize(ngrains);
		_grainMoments.resize(ngrains);
	}
	void reset() {
		_stressVoigt << 0., 0., 0., 0., 0., 0.;
		for (size_t i = 0; i < _grainForces.size(); i++) {
			_grainForces[i] << 0., 0., 0.;
			_grainMoments[i] << 0., 0., 0.;
		}
	}
	void operator+=(const GrainState3d & w) {
		_stressVoigt += w._stressVoigt; 
		for (size_t i = 0; i < _grainForces.size(); i++) {
			_grainForces[i] += w._grainForces[i];
			_grainMoments[i] += w._grainMoments[i];
		}
		
	}
	// member variables
	Vector6d _stressVoigt;			// macroscopic stress of assembly
	vector<Vector3d> _grainForces;	// forces on grains
	vector<Vector3d> _grainMoments;	// moments on grains
};

struct WallState3d {
	WallState3d() {_cWallForce = 0; _cWallContacts = 0;}
	WallState3d(const vector<Vector3d> & wtemp):
		_wallForces(wtemp) {_cWallForce = 0; _cWallContacts = 0;}
	void reset() {
		for (size_t i = 0; i < _wallForces.size(); i++) {
			_wallForces[i] << 0., 0., 0.;
			_wallContacts[i] = 0;
		}
		for (size_t i = 0; i < _ballWallForces.size(); i++) {
			_ballWallForces[i] << 0,0, 0;
			_ballWallContacts[i] = 0;
		}
		_cWallForce = 0;
		_cWallContacts = 0;
	}
	void resize(const int & nwalls) {
		_wallForces.resize(nwalls);
		_wallContacts.resize(nwalls);
	}
	void resize(const int & nwalls, const int & nshellwalls) {
		_wallForces.resize(nwalls);
		_wallContacts.resize(nwalls);
		_ballWallForces.resize(nshellwalls);
		_ballWallContacts.resize(nshellwalls);
	}
	void operator+=(const WallState3d & w) {
		for (size_t i = 0; i < _wallForces.size(); i++) {
			_wallForces[i] += w._wallForces[i];
			_wallContacts[i] += w._wallContacts[i];
		}
		for (size_t i = 0; i < _ballWallForces.size(); i++) {
			_ballWallForces[i] += w._ballWallForces[i];
			_ballWallContacts[i] += w._ballWallContacts[i];
		}
		_cWallForce += w._cWallForce;
		_cWallContacts += w._cWallContacts;
	}
	
	
	// member variables
	// forces/contacts on regular walls
	vector<size_t> _wallContacts; // number of contacts on walls
	vector<Vector3d> _wallForces; // forces on walls
	// forces/contacts on cylinder wall
	double _cWallForce; // force on cylinder wall
	size_t _cWallContacts; // number of contacts on cylinder wall
	// forces/contacts on shell wall
	vector<size_t> _ballWallContacts;
	vector<Vector3d> _ballWallForces;
};

struct CData {
	// member variables
	vector<Vector2i> _cpairs;
	vector<size_t> _nodes; // nodes of _cpairs[i](0)
	vector<Vector3d> _forces;
	vector<Vector3d> _normals;
	vector<Vector3d> _clocs;
	
	void operator+=(const CData & c) {
		_cpairs.insert( _cpairs.end(),	c._cpairs.begin(),	c._cpairs.end());
		_nodes.insert(_nodes.end(),		c._nodes.begin(),		c._nodes.end());
		_forces.insert( _forces.end(),	c._forces.begin(),	c._forces.end());
		_normals.insert(_normals.end(),	c._normals.begin(),	c._normals.end());
		_clocs.insert(_clocs.end(),		c._clocs.begin(),		c._clocs.end());
	}
};



#endif /* WORLDSTATES_H_ */
