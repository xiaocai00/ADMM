
#include <iostream>
#include "Grid.h"
#include "mpi.h"
#include "Definition.h"

using namespace std;
int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Barrier(MPI_COMM_WORLD);
  string fnameIn(argv[0]);
  string fnameOut(argv[1]);

  Grid grid;
  if (grid.GetRank() == 0)
    grid.Convert(fnameIn);
  grid.Save(fnameOut);
}
