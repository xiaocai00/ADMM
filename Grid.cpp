/*
 * GraphGenerator.cpp
 *
 *  Created on: Feb 1, 2013
 *      Author: xiaocai
 */

#include <sstream>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <climits>
#include <cstring>
#include <iostream>
#include <iostream>
#include <algorithm>
#include "mpi.h"
#include "Grid.h"
#include "pnetcdf.h"
#include "definition.h"

#define NDIMS 2
#define HANDLE_ERROR {                                \
  if (err != NC_NOERR)                             	 	\
  printf("Error at line %d (%s) at %d\n", __LINE__,   			\
      ncmpi_strerror(err), rank);                  					\
}

//#define DEBUG
//#define GENERATE 

/**
 *
 */
Grid::Grid(int globalRowSize, int globalColSize,
    int verticePotentialSize, int edgePotentialSize)
: GraphBase(verticePotentialSize, edgePotentialSize),
  numLocalRows(0), numLocalCols(0),
  numGlobalRows(globalRowSize), numGlobalCols(globalColSize),
  numLocalVerts(0), numLocalEdges(0), numLocalEdgesGhost(0), nprocs(0), rank(0),
  edgeList(NULL), request(NULL), commBuffer(NULL), rhoEdge(1.0)
{
  numGlobalVertices = globalRowSize * globalColSize;
  numGlobalEdges = globalRowSize * (globalColSize - 1) +
    globalColSize * (globalRowSize - 1);
  edgeVal.resize(numLocalEdges);
  part = ROW; 
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}

/**
 *
 */
Grid::Grid(int localRowSize, int localColSize,
    int globalRowSize, int globalColSize,
    int verticePotentialSize, int edgePotentialSize): 
  GraphBase(verticePotentialSize, edgePotentialSize),
  numLocalRows(localRowSize), numLocalCols(localColSize),
  numGlobalRows(globalRowSize), numGlobalCols(globalColSize),
  numLocalVerts(0), numLocalEdges(0), numLocalEdgesGhost(0), nprocs(0), rank(0),
  edgeList(NULL), request(NULL), commBuffer(NULL), rhoEdge(1.0)
{
  numGlobalVertices = globalRowSize * globalColSize;
  numGlobalEdges = globalRowSize * (globalColSize - 1) +
    globalColSize * (globalRowSize - 1);

  numLocalVerts = numLocalRows * numLocalCols;
  numLocalEdges = localRowSize * (localColSize - 1) +
    localColSize * (localRowSize - 1);
  numLocalEdgesGhost = numLocalEdges;
  edgeVal.resize(numLocalEdges);
  part = ROW; 
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}

Grid::Grid() : GraphBase(), edgeList(NULL), request(NULL), 
  commBuffer(NULL), rhoEdge(1.0)
{
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#if DEBUG
  std::ostringstream ss;
  ss << "log_" << rank << ".txt";
  logFile.open(ss.str().c_str()); 
#endif
}

void Grid::SetGlobalInfo(int globalRowSize, int globalColSize,
    int verticePotentialSize, int edgePotentialSize)
{
  numGlobalRows = globalRowSize;
  numGlobalCols = globalColSize;
  numVertPotentials = verticePotentialSize;
  numEdgePotentials = edgePotentialSize;

  numGlobalVertices = globalRowSize * globalColSize;
  numGlobalEdges = globalRowSize * (globalColSize - 1) +
    globalColSize * (globalRowSize - 1);
#if 0 
  int numLocalRows, int numLocalCols,
  numLocalVerts = numLocalRows * numLocalCols;
  numLocalEdges = localRowSize * (localColSize - 1) +
    localColSize * (localRowSize - 1);
  numLocalEdgesGhost = numLocalEdges;
  edgeVal.resize(numLocalEdges);
#endif
}

/**
 *
 */
Grid::~Grid()
{
#if DEBUG
  logFile.close();
#endif
  if(edgeList)
    delete edgeList;
  if (vertPotentials)
    delete vertPotentials;
  if (edgePotentials)
    delete edgePotentials;
  for (int i = 0; i < edgeVal.size(); i++)
    delete edgeVal[i];
}
  
/**
 *
 */
void Grid::Partition()
{
  int boundingBoxPerProc[nprocs * BBOXSIZE];
  int boundingBox[BBOXSIZE];

  /*
   * prepare the partition boundary for each processor
   */
  if (part == ROW) {
    int rowOffset = 0;
    int rowSize = numGlobalRows / nprocs;
    int extraRows = numGlobalRows % nprocs; 
    /*
     * initialize the bounding box by row partition 
     */
    if (rank < extraRows) {
      rowSize = numGlobalRows / nprocs + 1;
      rowOffset = rank * rowSize;
    } else {
      rowSize = numGlobalRows / nprocs;
      rowOffset = extraRows * ((numGlobalRows + nprocs - 1)/nprocs)
                  + (rank-extraRows)*rowSize; 
    }

    if (numLocalCols != numGlobalCols)
      numLocalCols = numGlobalCols;
    
    bbox.setBBox(rowOffset, 0, rowSize, numLocalCols);

    /**
     * set the local grid
     * ______
     * | | |
     * ______
     * | | |
     * ______
     * 3 * 4
     */
    numLocalRows = bbox.getBoxRowSize();
    numLocalCols = bbox.getBoxColSize();

    numLocalVerts = numLocalRows * numLocalCols;
    numLocalEdges = 2 * numLocalRows * numLocalCols - numLocalRows ;

    if ((numLocalRows + bbox.getUpperLeftRow()) == numGlobalRows)
      numLocalEdges -= numLocalCols;
    
    numLocalEdgesGhost = numLocalEdges;
#ifndef GENERATE
    if (rank != 0)
      numLocalEdgesGhost += numGlobalCols;
#endif

    // TODO: verify the number of edges and the number of vertices
    if (numLocalVerts > 0) {

      edgeList = new int[numLocalEdgesGhost * NUM_NODE];
      
      vertPotentials = new float[numVertPotentials * numLocalVerts];
      edgePotentials = new float[numEdgePotentials * numLocalEdges];
            
      edgeVal.resize(numLocalEdges);
      for (int i = 0; i < numLocalEdges; i++)
        edgeVal[i] = new Edge(numVertPotentials);
    }

#ifdef DEBUG
    logFile << "\n-------------- PARTITION RESULTS " 
            << " at Rank " << rank << "---------------\n";
    logFile << "bbox:               " << bbox << std::endl;
    logFile << "numLocalRows:       " << numLocalRows<< "\n";
    logFile << "numLocalCols:       " << numLocalCols<< "\n";
    logFile << "numLocalVerts:      " << numLocalVerts<< "\n";
    logFile << "numLocalEdges:      " << numLocalEdges<< "\n";
    logFile << "numLocalEdgesGhost: " << numLocalEdgesGhost << "\n";
#endif
  } else {
    logFile << "only supported for the case: " << numGlobalCols << " = "
      << numLocalCols << std::endl;
    MPI_Abort(MPI_COMM_WORLD, -1);	
  }

}

/**
 * Generate the graph
 */
void Grid::Generate() 
{
  int startID, endID, ie = 0;
  int rowStart = bbox.getUpperLeft().first;
  int colStart = bbox.getUpperLeft().second;
  int rowEnd = rowStart + bbox.getBoxSize().first;
  int colEnd = colStart + bbox.getBoxSize().second;
#ifdef DEBUG 
  logFile << "rank " << rank << " has bbox: " << bbox << std::endl;
#endif 
  std::vector<std::pair<int, int> > edgeVec(numLocalEdges);
  
  for (int i = rowStart ; i < rowEnd; i++) {
    /**
     * generate row edges
     */
    startID = i * numGlobalCols + colStart;
    for (int j = colStart; j < colEnd-1; j++) {
      endID = startID + 1;
      try {
        edgeVec[ie++] = (std::pair<int, int>(startID, endID));
      } catch (std::exception &e) {
        logFile << e.what() << std::endl;
      }
      startID++;
    }

    if (i == (numGlobalRows - 1)) {
      break;
    }
    
    /**
     * generate column edges
     */
    startID = i * numGlobalCols + colStart;
    for (int j = colStart; j < colEnd; j++) {
      endID = startID + numGlobalCols;
      try {
        edgeVec[ie++] = (std::pair<int, int>(startID, endID));
      } catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
      }
      startID++;
    }
  }
  
  edgeVec.resize(ie);
  numLocalEdges = ie;
  numLocalVerts = numLocalRows * numLocalCols;
  /*
   * copy vector to the contiguous linear buffer
   * maybe there is a better way to implement this
   */
  ie = 0;
  std::vector<std::pair<int, int> >::iterator iter;
  for(iter = edgeVec.begin(); iter != edgeVec.end(); iter++) {
    edgeList[ie++] = iter->first;
    edgeList[ie++] = iter->second;
  }

  srand(time(NULL));
  vertPotentials = new float[numLocalVerts * numVertPotentials];
  edgePotentials = new float[numLocalEdges * numEdgePotentials];

  for (int i = 0; i < numLocalVerts; i++) {
    for (int j = 0; j < numVertPotentials; j++) {
      vertPotentials[i * numVertPotentials + j] = ((float)rand())/RAND_MAX;
    }
  }

  for (int i = 0; i < numLocalEdges; i++) {
    for (int j = 0; j < numEdgePotentials; j++) {
      edgePotentials[i * numEdgePotentials + j] = ((float) rand())/RAND_MAX;
    }
  }

  // Initialize Neighbors
  InitNeighbors();
}

/**
 *
 */
void Grid::SetGlobalInfo(const std::string fname)
{
  int err, ncid;
  int type, gRows, gCols, gVertices, gEdges, nvPotentials, nePotentials;
  
  // open a file
  err = ncmpi_open (MPI_COMM_WORLD, fname.c_str(),
        NC_CLOBBER | NC_64BIT_DATA, MPI_INFO_NULL, &ncid);
  HANDLE_ERROR

  ncmpi_get_att_int(ncid, NC_GLOBAL, "graphType", &type);
  ncmpi_get_att_int(ncid, NC_GLOBAL, "numRows", &gRows);
  ncmpi_get_att_int(ncid, NC_GLOBAL, "numCols", &gCols);
  ncmpi_get_att_int(ncid, NC_GLOBAL, "numVertices", &gVertices);
  ncmpi_get_att_int(ncid, NC_GLOBAL, "numEdges", &gEdges);
  ncmpi_get_att_int(ncid, NC_GLOBAL, "numVerticePotentials", &nvPotentials);
  ncmpi_get_att_int(ncid, NC_GLOBAL, "numEdgePotentials", &nePotentials);
  HANDLE_ERROR

  err = ncmpi_close(ncid);
  HANDLE_ERROR
  
  numGlobalRows = gRows;
  numGlobalCols = gCols;
  numVertPotentials = nvPotentials;
  numEdgePotentials = nePotentials;
  
  numGlobalVertices = numGlobalRows * numGlobalCols;
  numGlobalEdges = numGlobalRows * (numGlobalCols - 1) +
                   numGlobalCols * (numGlobalRows - 1);

#ifdef DEBUG
  logFile << "\n-------------- GET GLOBAL INFO --------------\n";
  logFile << "numRows:            " << gRows << std::endl;
  logFile << "numCols:            " << gCols<< std::endl;
  logFile << "numVertices:        " << gVertices<< std::endl;
  logFile << "numEdges:           " << gEdges<< std::endl;
  logFile << "numVertPotentials:  " << nvPotentials << std::endl;
  logFile << "numEdgePotentials:  " << nePotentials << std::endl;
  logFile << "\n";
  logFile << "numGlobalRows:      " << numGlobalRows << std::endl;
  logFile << "numGlobalCols:      " << numGlobalCols << std::endl;
  logFile << "numGlobalVertices:  " << numGlobalVertices << std::endl;
  logFile << "numGlobalEdges:     " << numGlobalEdges << std::endl;
  logFile << "numVertPotentials:  " << numVertPotentials << std::endl;
  logFile << "numEdgePotentials:  " << numEdgePotentials << std::endl;
#endif
}

void Grid::Load(const std::string fname)
{
  int err, ncid, vid;
  err = ncmpi_open(MPI_COMM_WORLD, fname.c_str(),
      NC_CLOBBER | NC_64BIT_DATA, MPI_INFO_NULL, &ncid);
  long long starts[NDIMS];
  long long counts[NDIMS];
  assert(NDIMS>1);
  // load the vertPotentials
  starts[0] = bbox.getUpperLeftRow() * numGlobalCols;
  starts[1] = 0;
  counts[0] = numLocalVerts;
  counts[1] = numVertPotentials;
  
  if (starts[0] == numGlobalVertices) {
    starts[0] = 0;
    counts[0] = 0;
    counts[1] = 0;
  }

  std::string varname = "verticePotentials";

#ifdef DEBUG
  logFile << "\n-------------- LOAD DATA --------------\n";
  logFile << varname
          << "\n\tstart(" << starts[0] << ", " << starts[1] << "), " 
          << " end(" << counts[0] << ", " << counts[1] << ") " << std::endl;
#endif
  err = ncmpi_inq_varid(ncid, varname.c_str(), &vid);
  err = ncmpi_get_vara_float_all(ncid, vid, starts, counts, vertPotentials);
  HANDLE_ERROR
  
  // load the edgePotentials
  starts[0] = bbox.getUpperLeftRow() * (2 * numGlobalCols - 1);
  starts[1] = 0;
  counts[0] = numLocalEdges;
  counts[1] = numEdgePotentials;
  varname = "edgePotentials";

#ifdef DEBUG 
  logFile << varname
          << "\n\tstart(" << starts[0] << ", " << starts[1] << "), " 
          << " end(" << counts[0] << ", " << counts[1] << ") " << std::endl;
#endif

  err = ncmpi_inq_varid(ncid, varname.c_str(), &vid);
  err = ncmpi_get_vara_float_all(ncid, vid, starts, counts, edgePotentials);
  HANDLE_ERROR
  
  starts[0] = bbox.getUpperLeftRow() * (2 * numGlobalCols - 1);
  if (rank != 0)
    starts[0] -= numGlobalCols;
  starts[1] = 0;
  counts[0] = numLocalEdgesGhost;
  counts[1] = 2;
  if (starts[0] == numGlobalEdges) {
    counts[1] = 0;
    starts[0] = numGlobalRows-1;
  }
  
  // load the edgeList
  varname = "edgeList";
#ifdef DEBUG 
  logFile << varname <<" with ghostEdge " << numLocalEdgesGhost;
  logFile << "\n\tstart(" << starts[0] << ", " << starts[1] << "), "
          << "count(" << counts[0] << ", " << counts[1] << ") " << std::endl; 
  logFile << "-------------- END LOAD DATA --------------\n";
#endif
  err = ncmpi_inq_varid(ncid, varname.c_str(), &vid);
  err = ncmpi_get_vara_int_all(ncid, vid, starts, counts, edgeList);
  HANDLE_ERROR

  err = ncmpi_close(ncid);
  HANDLE_ERROR
}

void Grid::InitTree()
{
  if(rank == -1) {
    DisplaySelf();
    DisplayNeighbors();
  }
  
  // assign potentials to edge vector
  // ROW PARTITION
  int gvid1, gvid2, lvid1, lvid2, nedges1, nedges2;
  int offset = (rank == 0) ? 0 : numGlobalCols; 
  for(int i = 0; i < numLocalEdges; i++) {
    
    gvid1 = edgeList[2*i + offset];
    lvid1 = vSelfMapG2L[gvid1];
    gvid2 = edgeList[2*i+1 + offset];
    lvid2 = vSelfMapG2L[gvid2];
    nedges1 = internalNeighborMap[lvid1].size();
    nedges2 = internalNeighborMap[lvid2].size();
  
    // need normalization
    for (int j = 0; j< numVertPotentials; j++) {
      edgeVal[i]->x_potential_node[j] =
        vertPotentials[numVertPotentials*lvid1 + j] / (rhoEdge * nedges1);

      edgeVal[i]->x_potential_node[numVertPotentials + j] =
        vertPotentials[numVertPotentials*lvid2 + j] / (rhoEdge * nedges2);
    }

    for (int j = 0; j< numEdgePotentials; j++) {
      edgeVal[i]->x_potential_edge[j] =
        edgePotentials[numEdgePotentials*i+j]/ rhoEdge;
    }
  }
}

/**
 * store the graph
 */
void Grid::Save(const std::string fname)
{
  int err, ncid;

  err = ncmpi_create(MPI_COMM_WORLD, fname.c_str(),
      NC_NOWRITE, MPI_INFO_NULL, &ncid);
  HANDLE_ERROR

  int nodeDims[2], edgeDims[2], edgeListDims[2];
  ncmpi_def_dim(ncid, "numVertices", numGlobalVertices, &nodeDims[0]);
  ncmpi_def_dim(ncid, "numVerticePotentials", numVertPotentials, &nodeDims[1]);
  ncmpi_def_dim(ncid, "numEdges", numGlobalEdges, &edgeDims[0]);
  ncmpi_def_dim(ncid, "numEdgePotentials", numEdgePotentials, &edgeDims[1]);
  ncmpi_def_dim(ncid, "numEndPointsPerEdge", 2, &edgeListDims[1]);
  edgeListDims[0] = edgeDims[0];

  int vidNode, vidEdge, vidEdgeList;
  std::string varname = "edgeList";
  ncmpi_def_var(ncid, varname.c_str(), NC_INT, NDIMS, edgeListDims, &vidEdgeList);
  varname = "edgePotentials";
  ncmpi_def_var(ncid, varname.c_str(), NC_FLOAT, NDIMS, edgeDims, &vidEdge);
  varname = "verticePotentials";
  ncmpi_def_var(ncid, varname.c_str(), NC_FLOAT, NDIMS, nodeDims, &vidNode);

  ncmpi_put_att_text(ncid, NC_GLOBAL, "graphType", 5, "Grid");
  ncmpi_put_att_int(ncid, NC_GLOBAL, "numRows", NC_INT, 1, &numGlobalRows);
  ncmpi_put_att_int(ncid, NC_GLOBAL, "numCols", NC_INT, 1, &numGlobalCols);
  ncmpi_put_att_int(ncid, NC_GLOBAL, "numEdges", NC_INT, 1, &numGlobalEdges);
  ncmpi_put_att_int(ncid, NC_GLOBAL, "numVertices", NC_INT, 1, 
      &numGlobalVertices);
  ncmpi_put_att_int(ncid, NC_GLOBAL, "numVerticePotentials", NC_INT, 1,
      &numVertPotentials);
  ncmpi_put_att_int(ncid, NC_GLOBAL, "numEdgePotentials", NC_INT, 1,
      &numEdgePotentials);

  err = ncmpi_enddef(ncid);
  HANDLE_ERROR
    /**
     * write process
     */
  long long starts[NDIMS];
  long long counts[NDIMS];

  starts[0] = bbox.getUpperLeft().first * (2 * numGlobalCols - 1);
  starts[1] = 0;
  counts[0] = numLocalEdges;
  counts[1] = 2;
  err = ncmpi_put_vara_int_all(ncid, vidEdgeList, starts, counts, edgeList);
  HANDLE_ERROR

  starts[0] = bbox.getUpperLeft().first * numGlobalCols;
  starts[1] = 0;
  counts[0] = numLocalVerts;
  counts[1] = numVertPotentials;
  err = ncmpi_put_vara_float_all(ncid, vidNode, starts, counts, vertPotentials);
  HANDLE_ERROR

  starts[0] = bbox.getUpperLeft().first * (2 * numGlobalCols - 1);
  starts[1] = 0;
  counts[0] = numLocalEdges;
  counts[1] = numEdgePotentials;
  err = ncmpi_put_vara_float_all(ncid, vidEdge, starts, counts, edgePotentials);
  HANDLE_ERROR

  err = ncmpi_close(ncid);
  HANDLE_ERROR
}

/**
 *
 */
int Grid::GetProcForGVertByRow(int gvid)
{
  int nrows = (numGlobalRows + nprocs - 1) / nprocs;
  int vidPivot = nrows * (numGlobalRows % nprocs) * numGlobalCols;
  if (gvid < vidPivot){
    return gvid / numGlobalCols / nrows;
  } else {
    return (gvid-vidPivot) / numGlobalCols / (numGlobalRows/nprocs) + (numGlobalRows % nprocs);
  }
}

void Grid::InitNeighbors()
{
  switch(part) {
    case ROW:
      InitNeighborsByRow();
      break;
    case COL:
      std::cerr << "no implementation for COL\n";
      exit(-1);
    case BLOCK:
      std::cerr << "no implementation for BLOCK\n";
      exit(-1);
    default:
      break;
  }
}

/**
 *
 */
void Grid::AddToNeighborMap(int id, std::map<int, std::set<int> > &neighorMap)
{
  int pid, vidSelf, vidOther = edgeList[id];
  // starting vertex
  if (id & 0x1 == 0) {
    vidSelf = edgeList[id-1];
  } else {
    vidSelf = edgeList[id];
  }
  pid = GetProcForGVertByRow(vidOther);
  if (neighorMap.find(pid) != neighorMap.end()) {
    neighorMap[pid].insert(vidSelf);
  } else {
    std::set<int> newset;
    newset.insert(vidSelf);
    neighorMap.insert(std::make_pair(pid, newset));
  }
}

void Grid::DisplayNeighbors()
{
  std::set<int>::iterator siter;
  std::map<int, std::set<int> >::iterator miter;
  std::cout << "DisplayOutNeighbors " << std::endl;
  for (miter = outNeighborMap.begin(); miter != outNeighborMap.end(); miter++) {
    std::cout << miter->first << ": ";
    for (siter = miter->second.begin(); siter != miter->second.end(); siter++) {
      std::cout << (*siter) << ", "; 
    }
    std::cout << std::endl;
  }

  std::cout << "DisplayInNeighbors " << std::endl;
  for (miter = inNeighborMap.begin(); miter != inNeighborMap.end(); miter++) {
    std::cout << miter->first << ": ";
    for (siter = miter->second.begin(); siter != miter->second.end(); siter++) {
      std::cout << (*siter) << ", "; 
    }
    std::cout << std::endl;
  }
  
  std::cout << "DisplayInternalNeighbors " << std::endl;
  for (miter = internalNeighborMap.begin(); miter != internalNeighborMap.end(); miter++) {
    
    std::cout << miter->first << ": " ;
    for (siter = miter->second.begin(); siter != miter->second.end(); siter++) {
      std::cout << (*siter) << ", "; 
    }
    std::cout << std::endl;
  }
  std::cout << "EndofNeighbors " << std::endl;
}

void Grid::DisplaySelf()
{
  std::cout << "self vertices ID\n";
  for (int i = 0; i < vSelfVec.size(); i++) 
    std::cout << vSelfVec[i] << ", ";
  std::cout << std::endl;

  std::cout << "self edgelist\n";
  for (int i = rank==0?0:numGlobalCols; i < numLocalEdges; i++) 
    std::cout << edgeList[2*i] << ", " << edgeList[2*i+1] << std::endl;
 
  std::cout << "self G2L\n";
  std::map<int, int>::const_iterator miter;
  for (miter = vSelfMapG2L.begin(); miter != vSelfMapG2L.end(); miter++) 
    std::cout << miter->first << ", " << miter->second << std::endl;
  
}
/**
 *
 */
void Grid::InitNeighborsByRow()
{
  /**
   * build vertex mapping: local id --> global id
   */
  vSelfVec.resize(numLocalVerts); 
  int startVId = bbox.getUpperLeftRow() * numGlobalCols;
  int endVId = startVId + bbox.getBoxRowSize();
  int gid = startVId;
  int lid = 0;

  for (int i = startVId; i < endVId; i++) {
    for (int j = 0; j < numGlobalCols; j++) {
      vSelfVec[lid] = gid;
      vSelfMapG2L[gid++] = lid++;
    }
  }
 
  for (int i = 0; i < (numLocalEdgesGhost << 1); i++) {
    if (vSelfVec[0] > edgeList[i]) {
      AddToNeighborMap(i, inNeighborMap);
    } else if (vSelfVec[numLocalVerts-1] < edgeList[i]) {
      AddToNeighborMap(i, outNeighborMap);
    }
    if (i%2==0) {
      if (edgeList[i] <= vSelfVec[numLocalVerts - 1] && edgeList[i] >= vSelfVec[0])
        internalNeighborMap[edgeList[i]].insert(edgeList[i+1]);
      if (edgeList[i+1] <= vSelfVec[numLocalVerts - 1] && edgeList[i+1] >= vSelfVec[0])
        internalNeighborMap[edgeList[i+1]].insert(edgeList[i]);
    }
  }
}

void Grid::VerifyNeighbors() 
{
  std::set<int>::iterator siter;
  std::map<int, std::set<int> >::iterator miter;
  int prevVid;
  for (miter = outNeighborMap.begin(); miter != outNeighborMap.end(); miter++) {
    siter = miter->second.begin();
    prevVid = *siter;
    siter++;
    for (; siter != miter->second.end(); siter++) {

      if ((*siter) - prevVid != 1) {
        std::cerr << "boundary vertex id not contiguous for outNeighborMap:" 
          << *siter << ", " << prevVid << std::endl;
        break;
      }
      prevVid = *siter;
    }
  }

  for (miter = inNeighborMap.begin(); miter != inNeighborMap.end(); miter++) {
    siter = miter->second.begin();
    prevVid = *siter;
    siter++;
    for (; siter != miter->second.end(); siter++) {
      if ((*siter) - prevVid != 1) {
        std::cerr << "boundary vertex id not contiguous for inNeighborMap!\n";
        break;
      }
      prevVid = *siter;
    }
  }

#if 0
  //iter = std::upper_bound(startArr.begin(), startArr.end(), vid);
  //pid = iter - startArr.begin();
  if (rank == me) {
    std::set<int>::iterator siter;
    std::map<int, std::set<int> >::iterator miter;
    std::cout << "my neighbor : \n";
    for (miter = outNeighborMap.begin(); miter != outNeighborMap.end(); miter++) {
      std::cout << miter->first << ", ";
      for (siter = miter->second.begin(); siter != miter->second.end(); siter++)
        std::cout << *siter << " ";
      std::cout << std::endl;
    }
  }
#endif
}

/**
 *
 */
void Grid::Compute()
{
  //assert(edgeVal.size() == numLocalEdges);
  //std::cout << "edgeVal: " << edgeVal.size () << std::endl;
  for (int i = 0; i < edgeVal.size(); i++) {
    // update the lambda
    edgeVal[i]->updateLambda();
    // adjust the potential
    edgeVal[i]->adjustPotentials();
    // sum product
    edgeVal[i]->optimizeEdge();
  }
}

void Grid::FinalizeOptimization()
{
  FinalizeCommunication();
}

void Grid::InitOptimization()
{
  if (numLocalVerts == 0)
    return;
  InitNeighbors();
  InitTree();
  InitCommunication();
}

void Grid::InitCommunication()
{
  int totalRequests = inNeighborMap.size() + outNeighborMap.size();
  if (totalRequests == 0) {
    //std::cerr << "rank  " << rank << " has no neighbor!\n";
    return;
  }
  commBuffer = new float*[totalRequests];
  memset(commBuffer, 0, sizeof(float*) * totalRequests);
  request = (MPI_Request *) malloc(sizeof(MPI_Request) * totalRequests);
}

void Grid::FinalizeCommunication()
{
  int totalRequests = inNeighborMap.size() + outNeighborMap.size();
  if (totalRequests == 0)
    return;
  if (commBuffer) {
    delete commBuffer;
  }
  if (request) delete request;
}

/**
 * exchange the boundary vertcies' potentials
 * in order to get convergence. 
 */
void Grid::Communicate()
{
  if (numLocalVerts == 0)
    return;
  int idxSend = 0, idxRecv = 0, sidx = 0, vsize;
  // send the lower boudnary for row partition
  int edgeIdx = numLocalEdges - numGlobalCols;
  
  std::set<int>::iterator siter;
  std::map<int, std::set<int> > :: const_iterator citer;

  for (citer = outNeighborMap.begin(); citer != outNeighborMap.end(); citer++) {
    sidx = 0;
    vsize = citer->second.size() * numVertPotentials;
#if 0
    std::cout << "send rank: " << rank << "-->" << citer->first << ", " 
              << vsize << ", " << idxSend << ", " << edgeVal.size() << std::endl;
#endif
    commBuffer[idxSend] = new float[vsize];
    for (siter = citer->second.begin(); siter != citer->second.end(); siter++) {
      memcpy(commBuffer[idxSend] + (sidx++)*numVertPotentials, 
          edgeVal[edgeIdx++]->edge_mu_node, numVertPotentials);
    }
#if 1
    MPI_Isend(commBuffer[idxSend], vsize, MPI_FLOAT,
        citer->first, idxSend, MPI_COMM_WORLD, request + idxSend);
#endif
    idxSend++;
  }
  
  // send the uppper boudnary for row partition
  edgeIdx = -1;
  for (citer = inNeighborMap.begin(); citer != inNeighborMap.end(); citer++) {
    sidx = 0;
    vsize = citer->second.size() * numVertPotentials;
#if 0
    std::cout << "recv rank: " << rank << "-->" << citer->first << ", " 
              << vsize << ", " << idxSend + idxRecv << std::endl;
    std::cout << "edgeVal " << edgeVal.size() << ", " << numLocalEdges << std::endl;
#endif
    commBuffer[idxSend+idxRecv] = new float[vsize];
    for (siter = citer->second.begin(); siter != citer->second.end(); siter++) {
      try {
        if (sidx == citer->second.size() - 1)
          break;
        memcpy(edgeVal[++edgeIdx]->edge_mu_node, 
              commBuffer[idxSend + idxRecv] + (sidx++)*numVertPotentials, 
              numVertPotentials);
      } catch (std::exception &e) {
        std::cout << e.what() << std::endl;
      }
    }
    //the last node
    memcpy(edgeVal[edgeIdx]->edge_mu_node + numVertPotentials, 
        commBuffer[idxSend + idxRecv] + (sidx)*numVertPotentials, 
        numVertPotentials);
#if 1
    MPI_Irecv(commBuffer[idxRecv+idxSend], vsize, MPI_FLOAT,
        citer->first, idxRecv, MPI_COMM_WORLD, request+idxSend+idxRecv);
#endif
    idxRecv++;
  }

  MPI_Status status;
#if 1
  for (int i = 0; i < idxSend + idxRecv; i++) {
    MPI_Wait(request+i, &status);
  }
  for (int i = 0; i < idxSend + idxRecv; i++)
    if (commBuffer[i]) delete commBuffer[i];
#endif
}

std::ostream& operator<<(std::ostream& strm, const Grid &g)
{
  return strm << "G (" << g.numGlobalRows << ", " << g.numGlobalCols << ") "
    << "(" << g.numGlobalVertices << ", " << g.numGlobalEdges << ")";
}
