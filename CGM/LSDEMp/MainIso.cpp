
#include "definitions.h"
#include "Levelset3d.h"
#include "WorldStates.h"
#include "Grain3d.h"
#include "readInputFile3d.h"
#include "generationMethods.h"
#include "Wall3d.h"
#include "World3d.h"
#include "WallCylinder.h"

void runIso(string, double);

int main(int argc, char * argv[]){ 
    // comunicaion entre los nodos   	 
	MPI_Init(&argc, &argv);
	int numprocessors, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
	char name[50]; int count;
	MPI_Get_processor_name(name, &count);

    // If I am master node...
	if (rank == 0) {
	   cout << "Using " << numprocessors << " nodes with " << omp_get_max_threads()  << " threads/node for a total of " << numprocessors*omp_get_max_threads() << " processes." << endl;
	}

    runIso("P=5_1_tramp=0", 5.1);
    
    MPI_Finalize();
}

void runIso(string testID, double pp){
	// MPI stuff
	int numprocessors, rank; 
	MPI_Comm_size(MPI_COMM_WORLD, &numprocessors);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

    // If I am master node...
	if (rank == 0) {
		cout << "Running MainCycIso" << endl;
		cout << "Using " << numprocessors << " node." <<  endl;
	}
    
    // load grains
    size_t ngrains = 100; 
    vector<Grain3d> grains = generateGrainsFromFile("100grains.dat");

    //vector<Grain3d> grains = generateGrainsFromFile("jaramijo23-01-19.dat");
    Vector3d translation(-245.4364,-245.6927,0.0);

    // could be parallelized!!
    for (size_t i = 0; i < ngrains; i++) {
        grains[i].moveGrain(translation);
    }
   
    
    /*#pragma omp parallel for schedule(dynamic,200) num_threads(4)
    for (size_t g = 1; g < ngrains+1; g++) {
        //grains[g-1].changeMu(0);
        //grains[g-1].changeKn(1e4);
    }*/ 
        
    // find maximum radius and height (offset to ensure the sample lies on the xy plane smoothly)
    double maxr = 0;
    double maxh = 0;
    double minh = 1e3; // ADDED VARIABLE
    double r_to_add = 0; // ADDED VARIABLE
    for (size_t g = 0; g < ngrains; g++) {
        double r = sqrt(double( grains[g].getPosition()(0)*grains[g].getPosition()(0) + grains[g].getPosition()(1)*grains[g].getPosition()(1)));
        if (r > maxr) {
            maxr = r;
        }
        double h = grains[g].getPosition()(2) + grains[g].getRadius();
        if (h > maxh) {
            maxh = h;
        }
        //*********** ADDED SECTION   *********
        h = grains[g].getPosition()(2) - grains[g].getRadius();
        if (h < minh) {
            minh = h;
            r_to_add = grains[g].getRadius();
        }
        //******** END OF ADDED SECTION ********
    }
    
    //********************* ADDED SECTION  *******************
    cout << "max height: "<< maxh << endl;
    cout << "min height: "<<minh << endl;
    for (size_t i = 0; i < ngrains; i++) {
        grains[i].moveGrain(Vector3d(0,0,-minh));
    }
    minh = 1e3;
    maxh = 0;
    for (size_t i = 0; i < ngrains; i++) {
        double h = grains[i].getPosition()(2)- grains[i].getRadius();
        if (h < minh) {
            minh = h;
        }
        h = grains[i].getPosition()(2) + grains[i].getRadius();
        if (h > maxh) {
            maxh = h;
        }
    }
    cout << "new min height: " << minh << endl;
    cout << "radius of grain at bottom: "<<r_to_add<< endl;
    cout << "new max height: "<< maxh << endl;
	double ho = maxh;
    //******************* END OF ADDED SECTION ********************
  
    
	double kn = 30000; //grains[0].getKn();
	double mu = 0.65;
	double density = grains[0].getDensity();
	double dt = 0.2*(2*sqrt(.4*1000*density/kn));
	double gdamping = 1e4;
	double pressure = pp;//0.1; // pressure of 20 makes everything nan
	double stepPressure = 0; //pressure that will ramp up slowlyup to desired value
	double k = 0; // kn that will ramp up slowly un to the actual kn (grains)
    
	if (rank == 0 ) {
		cout << "Number of grains: " << ngrains << endl;
	}
	// set up walls
	vector<Wall3d> walls(2);
	
    double radiusInc = 8; 
    double height = maxh; //1118.; // FIXME
    double radius = maxr + radiusInc; // FIXME
	
	Vector3d botCorner(0.,0.,0.);
	Vector3d topCorner(0.,0.,height);
	// bottom wall
	walls[0] = Wall3d( Vector3d(0.,0.,1.), botCorner, kn, 0);
	// top wall
	walls[1] = Wall3d( Vector3d(0.,0.,-1.), topCorner, kn, 0);
	double ballRadius = 15;
	double kmn = 100;
	double kms = 100;
	WallBalls wallBalls(height, radius, ballRadius, pressure, density, kn, kmn, kms, dt, gdamping); // 99999 = id
	
	
	// create world
	World3d world(grains, walls, wallBalls, dt, gdamping);
	// clear grains from memory since it takes up a shitload of memory and creating world makes a copy of it anyway
	grains.clear();
	
    // set up file IO    
    string fileName = "wallsBallsIso/positions_iso_" + testID + ".dat";
    FILE * positions = fopen (fileName.c_str(),"w");
    fileName = "wallsBallsIso/kinetic_iso_" + testID + ".dat";
    FILE * kinetic = fopen (fileName.c_str(),"w");
    fileName = "wallsBallsIso/rotations_iso_" + testID + ".dat";
    FILE * quaternions = fopen (fileName.c_str(),"w");
    fileName = "wallsBallsIso/wallPos_iso_" + testID + ".dat";
    FILE * wallPositions = fopen(fileName.c_str(),"w");
    fileName = "wallsBallsIso/stress_iso_" +  testID + ".dat";
    FILE * stress = fopen(fileName.c_str(),"w");
    fileName = "wallsBallsIso/wallBallPos_iso_" +  testID + ".dat";
    FILE * wallBallPositions = fopen(fileName.c_str(),"w");
    fileName = "wallsBallsIso/test_params" +  testID + ".txt"; 
    FILE * params = fopen(fileName.c_str(), "w");

    FILE * initialDims = fopen ("wallsBallsTemp/initialDims.dat","w");    
    fprintf(initialDims, "%f \n%f", ho, radius);
    fclose(initialDims);

    // Files to store only last-iteration data
    FILE * lastPositions = fopen("wallsBallsTemp/lastPositions_iso.dat","w");
    FILE * lastRotations = fopen("wallsBallsTemp/lastRotations_iso.dat","w");
    FILE * WBLastPositions = fopen("wallsBallsTemp/wallsBallslastPositions_iso.dat","w");
    FILE * lastWallPos = fopen("wallsBallsTemp/lastWallPos_iso.dat","w");
    
	// to backup last-iteration data
	fileName = "wallsBallsBackups/lastPositions_iso_" +  testID + ".dat";
    FILE * lastPositions_backup = fopen(fileName.c_str(),"w");
	fileName = "wallsBallsBackups/lastRotations_iso_" +  testID + ".dat";
    FILE * lastRotations_backup = fopen(fileName.c_str(),"w");
	fileName = "wallsBallsBackups/wallsBallslastPositions_iso_" +  testID + ".dat";
    FILE * WBLastPositions_backup = fopen(fileName.c_str(),"w");
	fileName = "wallsBallsBackups/lastWallPos_iso_" +  testID + ".dat";
    FILE * lastWallPos_backup = fopen(fileName.c_str(),"w");
	
    size_t tiso = 5.0e3; //150; //1.5e5; // FIXME
    size_t tcomp = 0;
    // ramp pressure slowly
    size_t tramp = 3e3;//2e5;
    size_t ttopMove = 1e3;//2.75e4;
    double disp = -0.001;//- 0.001;
    
    /*fprintf(params,"Simulation Time: %f \ntop wall displacement rate: %f \ntop wall displacement time = %f \nPressure: %f \nRadius Increment: %f", (double)tiso, disp, (double)ttopMove, pressure, radiusInc);
    fclose(params);*/
        
	if ( rank == 0 ) {
    	cout << "kn = " << kn << endl;
    	cout << "mu = " << mu << endl;
    	cout << "gdamping = " << gdamping << endl;
        cout << "r = " << radius << endl;
	}
	 
    if (tramp > 0){   
        for (size_t i = 0; i < ngrains; i++) {
            world.getGrainsNonConst()[i].changeKn(1.);
            world.getGrainsNonConst()[i].changeKs(0.9*1.);
        }
    }

    double v0 = world.getWallBallsNonConst().findVolume();
    cout<< "v0: " << v0 << endl;
    double start = omp_get_wtime();
	
    int consolePrintStep = 100;
    int filePrintStep = 1000;
    //int logFileSwitchStep = 50000;

    //std::ofstream logFile;
    //logFile.open("logs/log1.dat");

	for (size_t t = 0; t < tiso+tcomp; t++) {

        //string logFileName = "logs/log";
        //if(t > 0 && t % logFileSwitchStep == 0){
          //  logFile.close();
            
            //std::ostringstream num ;
            //num << (t/logFileSwitchStep) + 1;

            //logFileName += num.str() + ".dat";
            //logFile.open( Name.c_str() );
        //}
	 	
		world.computeWorldState();
		
		if (rank == 0) {
			if (t % consolePrintStep == 0 || t == tiso-1) {				
				double duration = omp_get_wtime() - start;
                cout << "\n-----------------------------------------------------------------------------------------\n";
                //logFile << "\n-----------------------------------------------------------------------------------------\n";
		printf("Timestep %lu of %lu (%4.2f%% complete, %4.2f minutes, ~%4.2f remaining)\n", t+1, tiso+tcomp, 100*double(t+1)/double(tiso+tcomp), duration/60., duration/60.*(tiso+tcomp)/(t+1) - duration/60.);
		//logFile << "Timestep " << (t+1) << " of " << tiso+tcomp << " (" << 100*double(t+1)/double(tiso+tcomp) << " complete, " << (duration/60.) << " minutes, ~" << duration/60.*(tiso+tcomp)/(t+1) - duration/60. << " remaining)" << endl;
                fflush(stdout);
                cout << "height = " << world.getWalls()[1].getPosition()(2) << endl;
                cout << "top wall position: " << world.getWalls()[1].getPosition().transpose() << endl;
                cout << "v/v0 = " << world.getWallBalls().findVolume()/v0 << endl;
                cout << "stress = " << (world.getGrainState()._stressVoigt) << endl;
                cout << "force on top wall = " << world.getWallState()._wallForces[1](2) << endl;

                //logFile << "height = " << world.getWalls()[1].getPosition()(2) << endl;
                //logFile << "top wall position: " << world.getWalls()[1].getPosition().transpose() << endl;
                //logFile << "v/v0 = " << world.getWallBalls().findVolume()/v0 << endl;
                //logFile << "stress = " << (world.getGrainState()._stressVoigt) << endl;
		        //logFile << "force on top wall = " << world.getWallState()._wallForces[1](2) << endl;
				double ke = 0;
				for (size_t i = 0; i < world.getGrains().size(); i++) {
					ke += world.getGrains()[i].computeKineticEnergy();
					//ke += world.getGrains()[i].getVelocity().norm();
				}
		cout << "total kinetic energy = " << ke << endl;
                fprintf(kinetic, "%4.12f\n", ke);
        }
		
		// file writing 
        if(t % filePrintStep == 0 || t == tiso){
            //if (t == titestIDso+tcomp-1) {  //
				for (size_t i = 0; i < world.getGrains().size(); i++) {
	
					fprintf(positions, "%.12f %.12f %.12f\n", world.getGrains()[i].getPosition()(0), 
																				world.getGrains()[i].getPosition()(1), 
																				world.getGrains()[i].getPosition()(2));
					fprintf(quaternions, "%.12f %.12f %.12f %.12f\n", world.getGrains()[i].getQuat()(0), 
																					  	   world.getGrains()[i].getQuat()(1), 
																					  	   world.getGrains()[i].getQuat()(2), 
																					  	   world.getGrains()[i].getQuat()(3));
				}
				for (size_t i = 0; i < world.getWallBallsNonConst().getPositions().size(); i++) {
					fprintf(wallBallPositions, "%.12f %.12f %.12f\n", world.getWallBallsNonConst().getPositions()[i](0), 
																						  world.getWallBallsNonConst().getPositions()[i](1), 
																						  world.getWallBallsNonConst().getPositions()[i](2));
				}
			//}
			 
			//if (t >= tiso && t % 1000 == 0) {
				fprintf(wallPositions, "%.12f\n", world.getWalls()[1].getPosition()(2));
				//fprintf(wallPositions, "%4.12f\n", world.getWallCylinder().getRadius());
				//fprintf(stress, "%4.6g\n", world.getWallState()._wallForces[1](2));
				fprintf(stress,"%4.6g %4.6g %4.6g %4.6g %4.6g %4.6g\n", world.getGrainState()._stressVoigt(0), 
																						  world.getGrainState()._stressVoigt(1), 
																						  world.getGrainState()._stressVoigt(2),
																						  world.getGrainState()._stressVoigt(3),
																						  world.getGrainState()._stressVoigt(4),
                                                                                          world.getGrainState()._stressVoigt(5) );
            }
            //}
		} 
		
        // LAST ITERATION DATA File writing
        if(t == tiso-1){            
            for (size_t i = 0; i < world.getGrains().size(); i++) {                
                fprintf(lastPositions, "%.12f %.12f %.12f\n", world.getGrains()[i].getPosition()(0),
                        world.getGrains()[i].getPosition()(1),
                        world.getGrains()[i].getPosition()(2));
				fprintf(lastPositions_backup, "%.12f %.12f %.12f\n", world.getGrains()[i].getPosition()(0),
                        world.getGrains()[i].getPosition()(1),
                        world.getGrains()[i].getPosition()(2));		
						
                fprintf(lastRotations, "%.12f %.12f %.12f %.12f\n", world.getGrains()[i].getQuat()(0),
                        world.getGrains()[i].getQuat()(1),
                        world.getGrains()[i].getQuat()(2),
                        world.getGrains()[i].getQuat()(3));
				fprintf(lastRotations_backup, "%.12f %.12f %.12f %.12f\n", world.getGrains()[i].getQuat()(0),
                        world.getGrains()[i].getQuat()(1),
                        world.getGrains()[i].getQuat()(2),
                        world.getGrains()[i].getQuat()(3));
            }
           
            for (size_t i = 0; i < world.getWallBallsNonConst().getPositions().size(); i++) {
                fprintf( WBLastPositions, "%.12f %.12f %.12f\n", world.getWallBallsNonConst().getPositions()[i](0),
																						  world.getWallBallsNonConst().getPositions()[i](1),
																						  world.getWallBallsNonConst().getPositions()[i](2));
				fprintf( WBLastPositions_backup, "%.12f %.12f %.12f\n", world.getWallBallsNonConst().getPositions()[i](0),
																						  world.getWallBallsNonConst().getPositions()[i](1),
																						  world.getWallBallsNonConst().getPositions()[i](2));
            
			}
            
            fprintf(lastWallPos, "%4.12f\n", world.getWalls()[1].getPosition()(2));
			fprintf(lastWallPos_backup, "%4.12f\n", world.getWalls()[1].getPosition()(2));            
        }
        
        
        if (t < ttopMove) {
            world.getWallsNonConst()[1].moveWall(Vector3d(0,0,disp));
			world.getWallBallsNonConst().moveHeight(disp);
		}
		
		if (t < tramp) {
			for (size_t i = 0; i < ngrains; i++) {
				k =  double(t+1)*kn/double(tramp);
				world.getGrainsNonConst()[i].changeKn(k);
				world.getGrainsNonConst()[i].changeKs(0.9*k);
			}
            world.getWallBallsNonConst().changePressure(double(t+1)/double(tramp)*pressure);
			world.wallBallsTimestep();
            stepPressure = double(t+1)/double(tramp)*pressure;
        }
        else if (t<tiso) {
			world.wallBallsTimestep();
		}
		else {
			world.wallBallsTimestep();
		}
        if (t % consolePrintStep == 0 || t == tiso-1){
            cout<<"Pressure Applied: "<< world.getWallBallsNonConst().getPressure() << endl;
            cout<<"Grains Kn: "<< world.getGrainsNonConst()[0].getKn() << endl;
            cout<<"Grains Ks: "<< world.getGrainsNonConst()[0].getKs() << endl;

            //logFile << "Pressure Applied: " << world.getWallBallsNonConst().getPressure() << endl;
            //logFile << "Grains Kn: " << world.getGrainsNonConst()[0].getKn() << endl;
            //logFile << "Grains Ks: " << world.getGrainsNonConst()[0].getKs() << endl;
        }
        // take timestep
		world.grainTimestep();
		 
	} // end time loop
	 	
	fclose(positions);
        //logFile.close();
	fclose(quaternions);
	fclose(wallPositions);
	fclose(stress);
	fclose(wallBallPositions);
    fclose(kinetic);
    fclose(lastPositions);
    fclose(lastRotations);
    fclose(WBLastPositions);
    fclose(lastWallPos);
	fclose(lastPositions_backup);
    fclose(lastRotations_backup);
    fclose(WBLastPositions_backup);
    fclose(lastWallPos_backup);
    
	if (rank ==0) {
		double duration = omp_get_wtime() - start;
		printf("Time taken: %dh, %dm, %2.2fs\n", int(floor(duration/3600.)), -int(floor(duration/3600.))+int(floor(fmod(duration/60., 60.))), -floor(fmod(duration/60., 60.))*60. + fmod(duration, 3600.) );
		cout << "program terminated" << endl;
	}
	
    fprintf(params,"Simulation Time: %f \ntop wall displacement rate: %f \ntop wall displacement time = %f \nPressure Reached: %f \nGrains Kn reached: %f", (double)tiso, disp, (double)ttopMove, stepPressure, k);
    
    // find maximum radius and height
    maxr = 0;
    maxh = 0;
    for (size_t g = 0; g < ngrains; g++) {
        double r = sqrt(world.getGrains()[g].getPosition()(0)*world.getGrains()[g].getPosition()(0) + world.getGrains()[g].getPosition()(1)*world.getGrains()[g].getPosition()(1));
        if (r > maxr) {
            maxr = r;
        }
        double h = world.getGrains()[g].getPosition()(2) + world.getGrains()[g].getRadius();
        if (h > maxh) {
            maxh = h;
        }
    }
    
    double vf = world.getWallBallsNonConst().findVolume();
    double ratio = maxh/(2*maxr);
    fprintf(params,"\nVf: %f \nro: %f \nho: %f \nrf: %f \nhf: %f \nfinal H/D: %f", vf, radius, ho, maxr,maxh,ratio);
    //logFile.close();
    fclose(params);
}