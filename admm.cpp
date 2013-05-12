/*
 * admm.cpp
 *
 *  Created on: Feb 2, 2013
 *      Author: xiaocai
 */

#include <iostream>
#include "mpi.h"
#include "Grid.h"
#include "Definition.h"

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
#ifdef GENERATE
  if (grid.GetRank() == 0)
    cout << "generate grid\n";
  grid.SetPartitioner(ROW);
  grid.SetGlobalInfo(4, 3, 3, 9);
  grid.Partition();
  grid.Generate();
  grid.Save(fname);
#else
  grid.SetPartitioner(EDGE);
  //grid.SetGlobalInfo(10000, 1000, 3, 9);
  int irun = 0, totalRuns = 0;
  
  double iterTime = 0;
  double startTime = 0;
  
  double allocTime;
  double loadTime;
  double initTime;
  double optTime;
  
  double compTime = 0;
  double commTime = 0;
  double updateTime = 0;
   
  double maxCompTime = 0;
  double maxCommTime = 0;
  double maxUpdateTime = 0;
  
  double minCompTime = 0;
  double minCommTime = 0;
  double minUpdateTime = 0;
  
  PROFILE(startTime);
  grid.SetGlobalInfo(fname);
  PROFILE(allocTime);

  //grid.LoadExample();
  grid.Load(fname);
  PROFILE(loadTime);

  grid.InitOptimization();
  totalRuns = 10;
  PROFILE(initTime);
  // Optimization
  double obj, totalObj;
  while(irun++ < totalRuns) {
    iterTime = MPI_Wtime();
    grid.Compute();
    compTime += MPI_Wtime() - iterTime;
    
    iterTime = MPI_Wtime();
    grid.Communicate();
    MPI_Barrier(MPI_COMM_WORLD);
    commTime += MPI_Wtime() - iterTime;
    //grid.DisplayEdge();
    
    iterTime = MPI_Wtime();
    grid.Update();
    updateTime += MPI_Wtime() - iterTime;
    //grid.DisplayEdge();
    
    obj = grid.ComputeLinearObj();
    MPI_Reduce(&obj, &totalObj, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (grid.GetRank() == 0) {
      fprintf(stderr, "totalObj:       %lf\n", totalObj);
    }
  }
  grid.FinalizeOptimization();
  PROFILE(optTime);
  
  double runTime = MPI_Wtime() - startRun;
  MPI_GETMAX(commTime, maxCommTime);
  MPI_GETMIN(commTime, minCommTime);
  MPI_GETMAX(compTime, maxCompTime);
  MPI_GETMIN(compTime, minCompTime);
  MPI_GETMAX(updateTime, maxCompTime);
  MPI_GETMIN(updateTime, minCompTime);

  if (grid.GetRank()== 0) {
    fprintf(stderr, "processes:       %d\n", grid.GetNprocs());
    fprintf(stderr, "allocTime:       %lf\n", allocTime);
    fprintf(stderr, "loadTime:        %lf\n", loadTime);
    fprintf(stderr, "initTime:        %lf\n", initTime);
    fprintf(stderr, "optTime:         %lf\n", optTime);
    fprintf(stderr, "maxCompTime:     %lf\n", maxCompTime);
    fprintf(stderr, "maxCommTime:     %lf\n", maxCommTime);
    fprintf(stderr, "maxUpdateTime:   %lf\n", maxUpdateTime);
    fprintf(stderr, "minCompTime:     %lf\n", minCompTime);
    fprintf(stderr, "minCommTime:     %lf\n", minCommTime);
    fprintf(stderr, "minUpdateTime:   %lf\n", minUpdateTime);
    fprintf(stderr, "maxTotal:        %lf\n", runTime);
  }
#endif
  MPI_Finalize();
  return 0;
}
