
#include <iostream>
#include "Grid.h"
#include "mpi.h"
#include "Definition.h"

using namespace std;
int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Barrier(MPI_COMM_WORLD);
  if (argc < 3) {
    cerr << "./converter <infile> <outfile>\n";
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  string fnameIn(argv[1]);
  string fnameOut(argv[2]);

  Grid grid;
  if (grid.GetRank() == 0)
    grid.Convert(fnameIn);
  grid.Save(fnameOut);
}
