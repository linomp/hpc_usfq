/*
 * World3d.h
 *
 *  Created on: September 9, 2014
 *      Author: Reid
 */

#ifndef WORLD3D_H_
#define WORLD3D_H_

#include "definitions.h"
#include "WorldStates.h"
#include "Grain3d.h"
#include "Wall3d.h"
#include "WallCylinder.h"
//#include "WallShell.h"
#include "WallBalls.h"

// comment!
class World3d {
	
public:
	World3d() {_dt = 0; _gdamping = 0; _ngrains = 0; _nwalls = 0; _nshellwalls = 0;}
	
	// regular walls only
	World3d(const vector<Grain3d> & grains, const vector<Wall3d> & walls, const double & dt, const double & gdamping):
		_grains(grains), _walls(walls), _dt(dt), _gdamping(gdamping) {
		_ngrains = grains.size();
		_nwalls = walls.size();
		_nshellwalls = 0;
		_globalGrainState.resize(_ngrains);	_globalGrainState.reset();
		_globalWallState.resize(_nwalls);	_globalWallState.reset();
	}
	// 2 regular walls + wall cylinder
	World3d(const vector<Grain3d> & grains, const vector<Wall3d> & walls, const WallCylinder & wallCylinder, const double & dt, const double & gdamping):
		_grains(grains), _walls(walls), _wallCylinder(wallCylinder), _dt(dt), _gdamping(gdamping) {
//		_wallShell = WallShell();
		_ngrains = grains.size();
		_nwalls = walls.size();
		_nshellwalls = 0;
		_globalGrainState.resize(_ngrains);	_globalGrainState.reset();
		_globalWallState.resize(_nwalls,0);	_globalWallState.reset();
	}
	// 2 regular walls + wall shell
//	World3d(const vector<Grain3d> & grains, const vector<Wall3d> & walls, const WallShell & wallShell, const double & dt, const double & gdamping):
//		_grains(grains), _walls(walls), _wallShell(wallShell), _dt(dt), _gdamping(gdamping) {
//		_ngrains = grains.size();
//		_nwalls = walls.size();
//		_nshellwalls = wallShell.getNumWalls();
//		_globalGrainState.resize(_ngrains);					_globalGrainState.reset();
//		_globalWallState.resize(_nwalls, _nshellwalls);	_globalWallState.reset();
//	}
	// 2 regular walls + wall balls
	World3d(const vector<Grain3d> & grains, const vector<Wall3d> & walls, const WallBalls & wallBalls, const double & dt, const double & gdamping):
		_grains(grains), _walls(walls), _wallBalls(wallBalls), _dt(dt), _gdamping(gdamping) {
		_ngrains = grains.size();
		_nwalls = walls.size();
		_nshellwalls = wallBalls.getNumWalls();
		_globalGrainState.resize(_ngrains);					_globalGrainState.reset();
		_globalWallState.resize(_nwalls, _nshellwalls);	_globalWallState.reset();
	}
	
	// updates _globalGrainState and _globalWallState
	void computeWorldState() {
//		double s2 = omp_get_wtime();
		// zero out the global state and define temp variables
		_globalGrainState.reset();
		_globalWallState.reset();
		Vector3d force;
		Vector3d momenti;
		Vector3d momentj;
		Vector3d cmvec; 
		Vector6d stress;
		Vector3d fn, fs;
		int numprocessors, rank; 
		MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
//		GrainState3d threadGrainState;
//		WallState3d threadWallState;
//		#pragma omp declare reduction( + : WallState3d : omp_out += omp_in)
//		#pragma omp declare reduction( + : GrainState3d : omp_out += omp_in)
		#pragma omp parallel default(none) shared(numprocessors,rank,cout) firstprivate(force, momenti, momentj, cmvec, stress, fn, fs) // reduction(+:threadWallState,threadGrainState) //num_threads(1)
		{
			GrainState3d threadGrainState;
			WallState3d threadWallState;
//			for (size_t repeat = 0; repeat < 1000; repeat++) {
			threadGrainState.resize(_ngrains);
			threadGrainState.reset();
			if (_wallBalls.exists() ) {
				threadWallState.resize(_nwalls, _wallBalls.getNumWalls());
			}
			else {
				threadWallState.resize(_nwalls);
			}
			threadWallState.reset();
			size_t ncontacts = 0;
			// go through the grains and compute grain state
			
			#pragma omp for schedule(dynamic, 5) 
			for (size_t i = rank; i < _ngrains; i+=numprocessors) {
				// grain-grain contacts
				for (size_t j = i+1; j < _ngrains; j++) {
					if (_grains[i].bCircleCheck(_grains[j])) { 
						if (_grains[i].findInterparticleForceMoment(_grains[j], _dt, force, momentj, momenti)) {
							cmvec = _grains[i].getPosition() - _grains[j].getPosition();
							threadGrainState._grainForces[i] += force;
							threadGrainState._grainForces[j] -= force;
							threadGrainState._grainMoments[i] += momentj;
							threadGrainState._grainMoments[j] += momenti;
							threadGrainState._stressVoigt(0) += force(0)*cmvec(0);
							threadGrainState._stressVoigt(1) += force(1)*cmvec(1);
							threadGrainState._stressVoigt(2) += force(2)*cmvec(2);
							threadGrainState._stressVoigt(3) += 0.5*(force(1)*cmvec(2) + force(2)*cmvec(1));
							threadGrainState._stressVoigt(4) += 0.5*(force(2)*cmvec(0) + force(0)*cmvec(2));
							threadGrainState._stressVoigt(5) += 0.5*(force(1)*cmvec(0) + force(0)*cmvec(1));
						}
					}
				}
				// wall balls
				if (_wallBalls.exists()) {
					vector<Vector2i> possWalls = _wallBalls.bCircleContact(_grains[i]);
					if (possWalls.size() > 0) {
						if (_wallBalls.findWallForceMoment(_grains[i], possWalls, force, momenti, stress, threadWallState)) {
							threadGrainState._grainForces[i] += force;
							threadGrainState._grainMoments[i] += momenti;
							threadGrainState._stressVoigt += stress;
						}
					}
				}
				// cylinder wall
				if (_wallCylinder.getHeight() > 0) {
					if (_wallCylinder.bCircleContact(_grains[i])) {
						if (_wallCylinder.findWallForce(_grains[i], force, momenti, ncontacts, stress) ) {
							threadWallState._cWallContacts += ncontacts;
							threadWallState._cWallForce += force.norm();
							threadGrainState._grainForces[i] += force;
							threadGrainState._grainMoments[i] += momenti;
							threadGrainState._stressVoigt += stress;
						}
					}
				}
				// flat walls
				for (size_t j = 0; j < _nwalls; j++) {
					if (_walls[j].bCircleContact(_grains[i])) {
						if (_walls[j].findWallForceMomentFriction(_grains[i], force, momenti, ncontacts, stress, _dt) ) {
							threadGrainState._grainForces[i] += force;
							threadGrainState._grainMoments[i] += momenti;
							threadWallState._wallForces[j] -= force;
							threadWallState._wallContacts[j] += ncontacts;
							threadGrainState._stressVoigt += stress;
						}
					}
				} // end loop over wall j
			} // end loop over grains
			
			#pragma omp critical 
			{
				_globalWallState += threadWallState; // FIXME
				_globalGrainState += threadGrainState; // FIXME
			} // FIXME
		
		} // closes openmp parallel section
//		cout << "time to compute world state: rank: " << rank << ", time " << omp_get_wtime() - s2 << endl;
		// MPI calls  sendbuf       recvbuff                                   count       		type               op       comm
		MPI_Allreduce(MPI_IN_PLACE, _globalGrainState._grainForces[0].data(),  _ngrains*3,  	MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, _globalGrainState._grainMoments[0].data(), _ngrains*3,  	MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, _globalGrainState._stressVoigt.data(),     6,           	MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &_globalWallState._cWallContacts,          1,           	MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &_globalWallState._cWallForce,             1,           	MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, _globalWallState._wallContacts.data(),     _nwalls,     	MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, _globalWallState._wallForces[0].data(),    3*_nwalls,   	MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, _globalWallState._ballWallContacts.data(),_nshellwalls,		MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, _globalWallState._ballWallForces[0].data(),  3*_nshellwalls,	MPI_DOUBLE,        MPI_SUM, MPI_COMM_WORLD);
		
		// 2 walls + wall cylinder
		if (_walls.size() == 2 && _wallCylinder.getHeight() > 0 ) {
			_globalGrainState._stressVoigt /= M_PI*_wallCylinder.getRadius()*_wallCylinder.getRadius()*_wallCylinder.getHeight();
		}
		// 2 walls + wall shell
//		else if (_walls.size() == 2 && _wallShell.exists()) {
//			_globalGrainState._stressVoigt /= _wallShell.findVolume();
//		}
		// 2 walls + wall ball
		else if (_walls.size() == 2 && _wallBalls.exists()) {
			_globalGrainState._stressVoigt /= _wallBalls.findVolume();
		}
		// 6 wall cube
		else if (_walls.size() == 6) {
			_globalGrainState._stressVoigt /= (_walls[1].getPosition()(2) - _walls[0].getPosition()(2))*
													 (_walls[3].getPosition()(1) - _walls[2].getPosition()(1))*
													 (_walls[5].getPosition()(0) - _walls[4].getPosition()(0));
		}
		// print stuff for testing
//		char fname[100]; 
//		sprintf(fname, "worldState_%d.dat", numprocessors);
//		if (rank == 0) {
//			FILE * worldStateNonMPI = fopen (fname,"w");
////			for (size_t i = 0; i < _ngrains; i++) {
////				fprintf(worldStateNonMPI, "%f %f %f\n", _globalGrainState._grainForces[i](0), _globalGrainState._grainForces[i](1), _globalGrainState._grainForces[i](2));
////				fprintf(worldStateNonMPI, "%f %f %f\n", _globalGrainState._grainMoments[i](0), _globalGrainState._grainMoments[i](1), _globalGrainState._grainMoments[i](2));
////			}
//			fprintf(worldStateNonMPI, "%d\n", _globalWallState._cWallContacts);
//			fprintf(worldStateNonMPI, "%f\n", _globalWallState._cWallForce);
//			for (size_t i = 0; i < _nwalls; i++) {
//				fprintf(worldStateNonMPI, "%d\n", _globalWallState._wallContacts[i]);
//				fprintf(worldStateNonMPI, "%f\n %f\n %f\n", _globalWallState._wallForces[i](0), _globalWallState._wallForces[i](1), _globalWallState._wallForces[i](2) );
//			}
//			
//			fclose(worldStateNonMPI);
//		}
	} // end computeWorldState method
	
	CData computeCstate() const {
		CData cDataRank;
		int numprocessors, rank; 
		MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank==0) {
//			cout << "1" << endl;
		}
		#pragma omp parallel default(none) shared(numprocessors,rank,cout,cDataRank) // num_threads(1)
		{
			CData cDataThread;
			#pragma omp for schedule(dynamic, 5)
			for (size_t i = rank; i < _ngrains; i+=numprocessors) {
				for (size_t j = i+1; j < _ngrains; j++) {
					if (_grains[i].bCircleCheck(_grains[j])) { 
						CData cDataContact = _grains[i].findContactData(_grains[j]);
						if (cDataContact._clocs.size() > 0) {
							cDataThread += cDataContact;
						}
					}
				} // close grain subloop
			} // end loop over grains
			#pragma omp critical
			{
				cDataRank += cDataThread;
			}
		} // closes openmp parallel section
		return cDataRank;
	} // end computeCState method
	
	
	void applyBodyForce(Vector3d bodyForce) {
		for (size_t i = 0; i < _ngrains; i++) {
			_globalGrainState._grainForces[i] += bodyForce;
		}
	}
	
	void applyAcceleration(Vector3d acceleration) {
		for (size_t i = 0; i < _ngrains; i++) {
			_globalGrainState._grainForces[i] += _grains[i].getMass()*acceleration;
		}
	}
	 
	// take a timestep for each grain based off of the world's _globalGrainState
	void grainTimestep() {
		int numprocessors, rank; 
		MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
//		double s2 = omp_get_wtime();
//		for (size_t repeat = 0; repeat < 1000; repeat++) {
		// update center of mass, rotation, etc
		#pragma omp parallel for default(none) schedule(static,1) //num_threads(12)
		for (size_t i = 0; i < _ngrains; i++) {
			_grains[i].takeTimestep(_globalGrainState._grainForces[i], _globalGrainState._grainMoments[i], _gdamping, _dt);
//			_grains[i].updatePoints();
		}
//		 update points for grains on this rank
		#pragma omp parallel for default(none) shared(rank,numprocessors) schedule(static,1) //num_threads(12)
		for (size_t i = rank; i < _ngrains; i+=numprocessors) {
			_grains[i].updatePoints();
		}
	}
	
	void dragTopWall() {
		Vector3d force;
		Vector3d momenti;
		size_t ncontacts=0;
		Vector6d stress;
		Vector3d threadVelocity(0,0,0);
		Vector3d totalVelocity(0,0,0);
		size_t threadGcount = 0;
		size_t totalGcount = 0;
		#pragma omp parallel default(none) firstprivate(force,momenti,ncontacts,stress,threadVelocity,totalVelocity,threadGcount,totalGcount)
		{
			#pragma omp for schedule(dynamic,5)
			for (size_t i = 0; i < _ngrains; i+=1) {
				if (_walls[1].bCircleContact(_grains[i])) {
					if (_walls[1].findWallForceMomentFriction(_grains[i], force, momenti, ncontacts, stress, _dt) ) {
						threadVelocity += _grains[i].getVelocity();
						threadGcount++;
					}
				}
			}
		}
		#pragma omp critical
		{
			totalVelocity += threadVelocity;
			totalGcount += threadGcount;
		}
		double fac = 1.0;
		Vector3d projectedVelocity = totalVelocity - _walls[1].getNormal()*totalVelocity.dot(_walls[1].getNormal());
		_walls[1].changeVelocity( projectedVelocity*fac/(double)totalGcount );
	}
	 
	void wallCylinderTimestep() {
		_wallCylinder.takeTimestep(_globalWallState._cWallForce, _globalWallState._cWallContacts);
	}
	
	void wallBallsTimestep() {
		_wallBalls.takeTimestep(_globalWallState);
	}
	
	void topWallTimestep(const double & pressure) {
		_walls[1].takeTimestep(pressure*M_PI*_wallCylinder.getRadius()*_wallCylinder.getRadius() - _globalWallState._wallForces[1](2), _globalWallState._wallContacts[1]);
	}
	
	
//	void topWallTimestepShell(const double & pressure) {
//		_walls[1].takeTimestep(pressure*_wallShell.findTopArea() - _globalWallState._wallForces[1](2), _globalWallState._wallContacts[1]);
//	}
	
	void topWallTimestepBalls(const double & pressure) {
		_walls[1].takeTimestep(pressure*_wallBalls.findTopArea() - _globalWallState._wallForces[1](2), _globalWallState._wallContacts[1]);
	}
	
	void wallTimestep(const double & force, const size_t & wallId) {
		_walls[wallId].takeTimestep(force - _globalWallState._wallForces[wallId].norm(), _globalWallState._wallContacts[wallId]);
	}
	
	void topWallMoveWithGrains(const double & amount, const double & height) {
		int numprocessors, rank; 
		MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
		_walls[1].moveWall(Vector3d(0,0,amount)); 
		#pragma omp parallel for schedule(static,1) // num_threads(2)
		for (size_t i = rank; i < _ngrains; i+=numprocessors) {
			_grains[i].moveGrain(Vector3d(0,0,amount*_grains[i].getPosition()(2)/height ) );
		}
	}
	
	double findVolume() {
		double volume = 1;
		if (_walls.size() == 2 && _wallCylinder.getHeight() > 0 ) {
			volume = M_PI*_wallCylinder.getRadius()*_wallCylinder.getRadius()*_wallCylinder.getHeight();
		}
		else if (_walls.size() == 2 && _wallBalls.exists()) {
			volume = _wallBalls.findVolume();
		}
		// 6 walls
		else if (_walls.size() == 6) {
			volume = (_walls[1].getPosition()(2) - _walls[0].getPosition()(2))*
						(_walls[3].getPosition()(1) - _walls[2].getPosition()(1))*
						(_walls[5].getPosition()(0) - _walls[4].getPosition()(0));
		}
		return volume;
	}
	
	void changeDt(const double & dt) {
		_dt = dt;
	}
	void changeGdamping(const double & gdamping) {
		_gdamping = gdamping;
	}
	// changes kn of grains, walls, membrane
	void changeKnAll(const double & kn) {
		for (size_t i = 0; i < _grains.size(); i++) {
			_grains[i].changeKn(kn);
		}
		for (size_t i = 0; i < _walls.size(); i++) {
			_walls[i].changeKn(kn);
		}
		_wallBalls.changeKn(kn);
	}
	// changes ks of grains, walls
	void changeKsAll(const double & ks) {
		for (size_t i = 0; i < _grains.size(); i++) {
			_grains[i].changeKs(ks);
		}
		for (size_t i = 0; i < _walls.size(); i++) {
			_walls[i].changeKs(ks);
		}
	}
	// const get methods
	const vector<Grain3d> & getGrains() const {
		return _grains;
	}
	const vector<Wall3d> & getWalls() const {
		return _walls;
	}
	const GrainState3d & getGrainState() const {
		return _globalGrainState;
	}
	const WallState3d & getWallState() const {
		return _globalWallState;
	}
	const WallCylinder & getWallCylinder() const {
		return _wallCylinder;
	}
	const WallBalls & getWallBalls() const {
		return _wallBalls;
	}
	
//	 non const get methods (members can be changed such as wall positions)
	vector<Wall3d> & getWallsNonConst() {
		return _walls;
	}
	WallCylinder & getWallCylinderNonConst() {
		return _wallCylinder;
	}
//	WallShell & getWallShellNonConst() {
//		return _wallShell;
//	}
	WallBalls & getWallBallsNonConst() {
		return _wallBalls;
	}
	vector<Grain3d> & getGrainsNonConst() {
		return _grains;
	}
	
private:
	vector<Grain3d>	_grains;
	vector<Wall3d>		_walls;
//	WallShell			_wallShell;
	WallBalls			_wallBalls;
	WallCylinder		_wallCylinder;
	GrainState3d		_globalGrainState;	// grain state of entire assembly
	WallState3d			_globalWallState;		// wall state of entire assembly
	double _dt;			 // time increment
	double _gdamping; // global damping
	size_t _ngrains;		// number of grains
	size_t _nwalls;		// number of walls
	size_t _nshellwalls;
};


#endif
