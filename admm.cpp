/*
 * admm.cpp
 *
 *  Created on: Feb 2, 2013
 *      Author: xiaocai
 */

#include <iostream>
#include "Grid.h"
#include "mpi.h"
#include "definition.h"

using namespace std;

#define PROFILE(x) {MPI_Barrier(MPI_COMM_WORLD); \
    (x) = MPI_Wtime() - (startTime);\
    MPI_Barrier(MPI_COMM_WORLD);\
    startTime = MPI_Wtime();}

#define MPI_GETMAX(x, y)  MPI_Reduce(&(x), &(y), 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
#define MPI_GETAVG(x, y)  MPI_Reduce(&(x), &(y), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#define MPI_GETMIN(x, y)  MPI_Reduce(&(x), &(y), 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

int main(int argc, char* argv[])
{
  // Init Grid
  MPI_Init(&argc, &argv);
  MPI_Barrier(MPI_COMM_WORLD);
  double startRun = MPI_Wtime();
  string fname("climate.nc");
  Grid grid;
  grid.SetPartitioner(ROW);

#ifdef GENERATE
  grid.SetGlobalInfo(10000, 1000, 3, 9);
  grid.Partition();
  grid.Generate();
  grid.Save(fname);
#else
  //grid.SetGlobalInfo(10000, 1000, 3, 9);
  int irun = 0, totalRuns = 0;
  
  double globalTime = 0; 
  double iterTime = 0;
  double startTime = 0;
  
  double allocTime;
  double loadTime;
  double neighborTime;
  double optTime;
  
  double commTime = 0;
  double minCommTime = 0;
  double maxCommTime = 0;
  
  double compTime = 0;
  double minCompTime = 0;
  double maxCompTime = 0;
 
  PROFILE(startTime);
  grid.SetGlobalInfo(fname);

  grid.Partition();
  PROFILE(allocTime);

  grid.Load(fname);
  PROFILE(loadTime);

  grid.InitOptimization();
  PROFILE(neighborTime);
 
  // Optimization
  while(irun++ < totalRuns) {
    iterTime = MPI_Wtime();
    grid.Compute();
    compTime += MPI_Wtime() - iterTime;
    iterTime = MPI_Wtime();
    grid.Communicate();
    MPI_Barrier(MPI_COMM_WORLD);
    commTime += MPI_Wtime() - iterTime;
  }
 
  grid.FinalizeOptimization();
  PROFILE(optTime);
  MPI_Barrier(MPI_COMM_WORLD);
  double runTime = MPI_Wtime() - startRun;

  MPI_GETMAX(commTime, maxCommTime);
  MPI_GETMIN(commTime, minCommTime);
  MPI_GETMAX(compTime, maxCompTime);
  MPI_GETMIN(compTime, minCompTime);

  if (grid.GetRank()== 0) {
    fprintf(stderr, "processes:       %d\n", grid.GetNprocs());
    fprintf(stderr, "globalTime:      %lf\n", globalTime);
    fprintf(stderr, "allocTime:       %lf\n", allocTime);
    fprintf(stderr, "loadTime:        %lf\n", loadTime);
    fprintf(stderr, "neighborTime:    %lf\n", neighborTime);
    fprintf(stderr, "optTime:         %lf\n", optTime);
    fprintf(stderr, "maxCompTime:     %lf\n", maxCompTime);
    fprintf(stderr, "minCompTime:     %lf\n", minCompTime);
    fprintf(stderr, "maxCommTime:     %lf\n", maxCommTime);
    fprintf(stderr, "minCommTime:     %lf\n", minCommTime);
    fprintf(stderr, "maxTotal:        %lf\n", runTime);
  }
#endif
  MPI_Finalize();
  return 0;
}
