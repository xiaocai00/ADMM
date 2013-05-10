/*
 * GraphGenerator.cpp
 *
 *  Created on: Feb 1, 2013
 *      Author: xiaocai
 */
#include <iomanip>
#include <utility>      // std::pair
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
#include <fstream>
#include "pnetcdf.h"
#include "Definition.h"
#include <limits>       // std::numeric_limits
#define NDIMS 2
#define HANDLE_ERROR {                                  \
  if (err != NC_NOERR)                             	 	\
  printf("Error at line %d (%s) at %d\n", __LINE__,   	\
      ncmpi_strerror(err), rank);                  	    \
}

//#define DEBUG
//#define GENERATE 

/**
 *
 */
Grid::Grid(int globalRowSize, int globalColSize, int verticePotentialSize,
    int edgePotentialSize) :
  GraphBase(verticePotentialSize, edgePotentialSize), numLocalRows(0), numLocalCols(
      0), numGlobalRows(globalRowSize), numGlobalCols(globalColSize), numLocalVerts(
        0), numLocalEdges(0), numLocalEdgesGhost(0), numProcs(0), rank(0), edgeList(
          NULL), request(NULL), commBuffer(NULL), aux(NULL), edgeOp(NULL), rhoEdge(1.0),
          penalty(0.1), penaltyPrime(0.5) {
          numGlobalVerts = globalRowSize * globalColSize;
          numGlobalEdges = globalRowSize * (globalColSize - 1)
            + globalColSize * (globalRowSize - 1);
          edgeVal.resize(numLocalEdges);
          part = ROW;
          MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
          MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        }

/**
 *
 */
Grid::Grid(int localRowSize, int localColSize, int globalRowSize,
    int globalColSize, int verticePotentialSize, int edgePotentialSize) :
  GraphBase(verticePotentialSize, edgePotentialSize), numLocalRows(
      localRowSize), numLocalCols(localColSize), numGlobalRows(globalRowSize), 
  numGlobalCols(globalColSize), numLocalVerts(0), 
  numLocalEdges(0), numLocalEdgesGhost(0), 
  numProcs(0), rank(0), edgeList(NULL), request(NULL), 
  commBuffer(NULL), aux(NULL), edgeOp(NULL), rhoEdge(1.0),
  penalty(0.1), penaltyPrime(0.5) {
    numGlobalVerts = globalRowSize * globalColSize;
    numGlobalEdges = globalRowSize * (globalColSize - 1)
      + globalColSize * (globalRowSize - 1);

    numLocalVerts = numLocalRows * numLocalCols;
    numLocalEdges = localRowSize * (localColSize - 1)
      + localColSize * (localRowSize - 1);
    numLocalEdgesGhost = numLocalEdges;
    edgeVal.resize(numLocalEdges);
    part = ROW;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }

Grid::Grid() :
  GraphBase(), edgeList(NULL), request(NULL), commBuffer(NULL), 
  aux(NULL), edgeOp(NULL), rhoEdge(1.0),
  penalty(0.1), penaltyPrime(0.5) {
    numLocalEdges = 0;
    numLocalVerts = 0;
    numVertPotentials = 0;
    numEdgePotentials = 0;
    numGlobalRows = -1;
    numGlobalCols = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#if DEBUG
    std::ostringstream ss;
    ss << "log_" << rank << ".txt";
    logFile.open(ss.str().c_str());
#endif
  }

void Grid::SetGlobalInfo(int globalRowSize, int globalColSize,
    int verticePotentialSize, int edgePotentialSize) {

  numGlobalRows = globalRowSize;
  numGlobalCols = globalColSize;
  numVertPotentials = verticePotentialSize;
  numEdgePotentials = edgePotentialSize;

  numGlobalVerts = globalRowSize * globalColSize;
  numGlobalEdges = globalRowSize * (globalColSize - 1)
    + globalColSize * (globalRowSize - 1);

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
Grid::~Grid() {
#if DEBUG
  logFile.close();
#endif
  if (edgeList)
    delete edgeList;
  if (vertPotentials)
    delete vertPotentials;
  if (edgePotentials)
    delete edgePotentials;
  for (int i = 0; i < edgeVal.size(); i++)
    delete edgeVal[i];
  // release the memory
  for (int i = 0; i < vertexMuArr.size(); i++) {
    if (vertexMuArr[i])
      delete vertexMuArr[i];
  }
  if (aux)
    delete aux;
  if (edgeOp)
    delete edgeOp;
}

void Grid::ResetVertexMuArr() {
  for (int i = 0; i < numLocalVerts; i++) {
    memset(vertexMuArr[i], 0, sizeof(float) * numVertPotentials);
  }
}

void Grid::UpdateMuArrBeforeComm() {
  ResetVertexMuArr();
#if 1
  // update all the local information
  for (int i = 0; i < (numLocalEdges << 1); i++) {
    int nodeid = vSelfMapG2L[edgeList[i]];
    for (int j = 0; j < numVertPotentials; j++) {
      vertexMuArr[nodeid][j] += edgeVal[i / 2]->x_mu_node[j
        + numVertPotentials * i % 2];
    }
  }
#endif
}

void Grid::UpdateMuArrAfterComm() {
  // add the neighbors' information
  for (int i = 0; i < neighbors.size(); i++) {
    for (int j = 0; j < outNeighborMapEdge[i].size(); j++) {
      int nodeid = vSelfMapG2L[outNeighborMapEdge[i][j]];
      for (int k = 0; k < numVertPotentials; k++) {
        vertexMuArr[nodeid][k] += commBuffer[i][j*numVertPotentials+k] * cntNeighborMapEdge[i][j];
      }
    }
  }
  // this is a hack implementation, since each edge might have different
  // average the value
  for (int i = 0; i < numLocalVerts; i++) {
    for (int j = 0; j < numVertPotentials; j++) {
      vertexMuArr[i][j] /= vertCnt[i];
    }
  }
}

/**
 *
 */
void Grid::Update() {
  UpdateMuArrAfterComm();  
  for (int i = 0; i < (numLocalEdges << 1); i++) {
    int nodeid = vSelfMapG2L[edgeList[i]];
    memcpy(edgeVal[i/2]->edge_mu_node + (i % 2) * numVertPotentials,
        vertexMuArr[nodeid], numVertPotentials * sizeof(float));
  }
}

/**
 *
 */
void Grid::InitNeighborsByEdge() {
  int sbuf[2];
  sbuf[0] = vertexStart;
  sbuf[1] = vertexEnd;

  int *rbuf = new int[2 * numProcs];
  // broadcast my offset range
  MPI_Allgather(sbuf, 2, MPI_INT, rbuf, 2, MPI_INT, MPI_COMM_WORLD );

  for (int i = 0; i < numProcs; i++) {
    if (i == rank || sbuf[0] > rbuf[2 * i + 1] || sbuf[1] < rbuf[2 * i])
      continue;
    neighbors.push_back(i);
  }
  
  // overlap based on edges
  int totalCnt = 0;
  outNeighborMapEdge.resize(neighbors.size());
  
  for (int i = 0; i < neighbors.size(); i++) {
    int overlapStart, overlapLen;
    // left
    int nid = neighbors[i];
    if (sbuf[0] < rbuf[nid * 2]) {
      overlapStart = rbuf[nid * 2];
      overlapLen = sbuf[1] - rbuf[nid * 2] + 1;
    } else {
      overlapStart = sbuf[0];
      overlapLen = rbuf[nid * 2 + 1] - sbuf[0] + 1;
    }
    for (int j = overlapStart; j < overlapLen + overlapStart; j++) {
      outNeighborMapEdge[i].push_back(j);
    }
    totalCnt += outNeighborMapEdge[i].size();
  }
  MPI_Request request;
  std::vector<MPI_Request> requests(neighbors.size() * 2);
  requests.resize(2*neighbors.size());
  
  int* recvBuf = new int[totalCnt];
  if (2*numProcs < totalCnt) { 
    delete rbuf;
    rbuf = new int[totalCnt];
  }
  totalCnt = 0;
  int reqidx = 0;
  for (int i = 0; i < neighbors.size(); i++) {
    int nid = neighbors[i];
    int cnt = outNeighborMapEdge[i].size();
    for (int idx = 0; idx < outNeighborMapEdge[i].size(); idx++) {
      int vid = vSelfMapG2L[outNeighborMapEdge[i][idx]];
      rbuf[totalCnt + idx] = vertCnt[vid];
    }
    logFile << std::endl;
    
    MPI_Isend(rbuf + totalCnt, cnt, 
        MPI_INT, nid, rank, MPI_COMM_WORLD, &requests[reqidx++]);
    MPI_Irecv(recvBuf + totalCnt, cnt, MPI_INT, nid, nid, MPI_COMM_WORLD, 
        &requests[reqidx++]);
    
    totalCnt += cnt;
  }
  MPI_Status status;
  
#if 1
  for (int i = 0; i < reqidx; i++)
    MPI_Wait(&requests[i], &status);
  
  totalCnt = 0;
  cntNeighborMapEdge.resize(neighbors.size());
  for (int i = 0; i < neighbors.size(); i++) {
    int nid = neighbors[i];
    int cnt = outNeighborMapEdge[i].size();
    for (int j = 0; j < cnt; j++) {
      int vcnt = recvBuf[j + totalCnt];
      int vid = outNeighborMapEdge[i][j];
      vertCnt[vSelfMapG2L[vid]] += vcnt;
      cntNeighborMapEdge[i].push_back(vcnt);
    }
    totalCnt += cnt;
  }

  totalCnt = 0;
  for (int i = 0; i < numLocalVerts; i++) {
      totalCnt += vertCnt[i];
  }

  delete rbuf;
  delete recvBuf;
#endif
}

void Grid::DisplayNeighborsByEdge() {
  logFile << "---------------------\n";
  logFile << "DisplayOutNeighbors " << neighbors.size() << std::endl;
  for (int i = 0; i < neighbors.size(); i++) {
   logFile << neighbors[i]<< "(" << outNeighborMapEdge[i].size() << "): ";
    for (int j = 0; j < outNeighborMapEdge[i].size(); j++) {
      logFile << "(" << outNeighborMapEdge[i][j] << ","
        << cntNeighborMapEdge[i][j] << ")\t";
    }
    logFile << std::endl;
  }
}

void Grid::Load(std::string fname) {
  if (part == ROW) {
    PartitionByRow();
    LoadByRow(fname);
  } else if (part == EDGE) {
    PartitionByEdge(fname);
  } else {
    std::cerr << "only supported for the case: " << numGlobalCols << " = "
      << numLocalCols << "\n";
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

#ifdef DEBUG
  logFile << "\n-------------- PARTITION RESULTS " << "---------------\n";
  //logFile << "bbox:               " << bbox << std::endl;
  logFile << "numLocalRows:       " << numLocalRows << "\n";
  logFile << "numLocalCols:       " << numLocalCols << "\n";
  logFile << "numLocalVerts:      " << numLocalVerts << "\n";
  logFile << "numLocalEdges:      " << numLocalEdges << "\n";
  logFile << "numLocalEdgesGhost: " << numLocalEdgesGhost << "\n";
#endif
}

void Grid::PartitionByRow() {
  /*
   * prepare the partition boundary for each processor
   */
  int boundingBoxPerProc[numProcs * BBOXSIZE];
  int boundingBox[BBOXSIZE];
  int rowOffset = 0;
  int rowSize = numGlobalRows / numProcs;
  int extraRows = numGlobalRows % numProcs;
  /*
   * initialize the bounding box by row partition
   */
  if (rank < extraRows) {
    rowSize = numGlobalRows / numProcs + 1;
    rowOffset = rank * rowSize;
  } else {
    rowSize = numGlobalRows / numProcs;
    rowOffset = extraRows * ((numGlobalRows + numProcs - 1) / numProcs)
      + (rank - extraRows) * rowSize;
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
  numLocalEdges = 2 * numLocalRows * numLocalCols - numLocalRows;

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
}

void getVerticesRange(int *edgeList, int count, int &minId, int &maxId) {
  minId = std::numeric_limits<int>::max();
  maxId = 0;
  for (int i = 0; i < count; i++) {
    if (edgeList[i] < minId) {
      minId = edgeList[i];
    }
    if (edgeList[i] > maxId) {
      maxId = edgeList[i];
    }
  }
}

void Grid::PartitionByEdge(std::string fname) {
  numLocalEdges = numGlobalEdges / numProcs;
  int divid = numGlobalEdges % numProcs;
  int edgeOffset;
  if (rank < divid) {
    numLocalEdges += 1;
    edgeOffset = numLocalEdges * rank;
  } else {
    edgeOffset = (numLocalEdges + 1) * divid + numLocalEdges * (rank - divid);
  }
  numLocalEdgesGhost = numLocalEdges;
  edgeList = new int[numLocalEdges * NUM_NODE];
  edgePotentials = new float[numEdgePotentials * numLocalEdges];
  edgeVal.resize(numLocalEdges);
  for (int i = 0; i < numLocalEdges; i++)
    edgeVal[i] = new Edge(numVertPotentials);

  int err, ncid, vid;

  err = ncmpi_open(MPI_COMM_WORLD, fname.c_str(), NC_CLOBBER | NC_64BIT_DATA,
      MPI_INFO_NULL, &ncid);
  long long starts[NDIMS];
  long long counts[NDIMS];
  assert(NDIMS>1);

  //##################################
  // load the edgeList
  //##################################
  std::string varname = "edgeList";
  starts[0] = edgeOffset;
  starts[1] = 0;
  counts[0] = numLocalEdges;
  counts[1] = 2;
  if (starts[0] == numGlobalEdges) {
    counts[1] = 0;
    starts[0] = numGlobalRows - 1;
  }

#ifdef DEBUG
  logFile << "\n-------------- LOAD DATA --------------\n";
  logFile << varname << " with Edge " << numLocalEdges;
  logFile << "\n\tstart(" << starts[0] << ", " << starts[1] << "),  count("
    << counts[0] << ", " << counts[1] << ") " << std::endl;
#endif
  err = ncmpi_inq_varid(ncid, varname.c_str(), &vid);
  err = ncmpi_get_vara_int_all(ncid, vid, starts, counts, edgeList);
  HANDLE_ERROR

  //##################################
  //load the vertPotentials
  //##################################
  
  GetVerts();
  //getVerticesRange(edgeList, 2 * numLocalEdges, vertexStart, vertexEnd);
  //numLocalVerts = vertexEnd - vertexStart + 1;
#ifdef DEBUG
  logFile << "\n-------------- GET LOCAL INFO --------------\n";
  logFile << "vertexRange: " << vertexStart << ", " << vertexEnd << std::endl;
  logFile << "numLocalEdges:     " << numLocalEdges<< std::endl;
  logFile << "numLocalVertices:  " << numLocalVerts << std::endl;
#endif

  starts[0] = vertexStart;
  starts[1] = 0;
  counts[0] = numLocalVerts;
  counts[1] = numVertPotentials;
  varname = "verticePotentials";
  vertPotentials = new float[numVertPotentials * numLocalVerts];

#ifdef DEBUG
  logFile << varname << "\n\tstart(" << starts[0] << ", " << starts[1]
    << "), " << " count(" << counts[0] << ", " << counts[1] << ") "
    << std::endl;
#endif
  err = ncmpi_inq_varid(ncid, varname.c_str(), &vid);
  err = ncmpi_get_vara_float_all(ncid, vid, starts, counts, vertPotentials);
  HANDLE_ERROR

  //##################################
  // load the edgePotentials
  //##################################
  starts[0] = edgeOffset;
  starts[1] = 0;
  counts[0] = numLocalEdges;
  counts[1] = numEdgePotentials;
  varname = "edgePotentials";

#ifdef DEBUG
  logFile << varname << "\n\tstart(" << starts[0] << ", " << starts[1]
    << "), " << " end(" << counts[0] << ", " << counts[1] << ") "
    << std::endl;
#endif

  err = ncmpi_inq_varid(ncid, varname.c_str(), &vid);
  err = ncmpi_get_vara_float_all(ncid, vid, starts, counts, edgePotentials);
  HANDLE_ERROR
  err = ncmpi_close(ncid);
  HANDLE_ERROR
 
  int k = 0;
  for (int i = 0; i < numLocalVerts; i++) {
    for (int j = 0; j < numVertPotentials; j++) {
      vertPotentials[k] = -vertPotentials[k];
      k++;
    }
  } 
  k = 0; 
  for (int i = 0; i < numLocalEdges; i++) {
    for (int j = 0; j < numEdgePotentials; j++) {
      edgePotentials[k] = -edgePotentials[k];
      k++;
    }
  } 
  
  numLocalRows = -1;
  numLocalCols = -1; 
  numLocalEdgesGhost = -1;
}

/**
 * Generate the graph
 */
void Grid::Generate() {
  int startID, endID, ie = 0;
  int rowStart = bbox.getUpperLeft().first;
  int colStart = bbox.getUpperLeft().second;
  int rowEnd = rowStart + bbox.getBoxSize().first;
  int colEnd = colStart + bbox.getBoxSize().second;
#ifdef DEBUG 
  logFile << "rank " << rank << " has bbox: " << bbox << std::endl;
#endif 
  std::vector<std::pair<int, int> > edgeVec(numLocalEdges);

  for (int i = rowStart; i < rowEnd; i++) {
    /**
     * generate row edges
     */
    startID = i * numGlobalCols + colStart;
    for (int j = colStart; j < colEnd - 1; j++) {
      endID = startID + 1;
      try {
        edgeVec[ie++] = (std::pair<int, int>(startID, endID));
      } catch (std::exception &e) {
#ifdef DEBUG
        logFile << e.what() << std::endl;
#endif
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
  for (iter = edgeVec.begin(); iter != edgeVec.end(); iter++) {
    edgeList[ie++] = iter->first;
    edgeList[ie++] = iter->second;
  }

  srand(time(NULL));
  vertPotentials = new float[numLocalVerts * numVertPotentials];
  edgePotentials = new float[numLocalEdges * numEdgePotentials];

  for (int i = 0; i < numLocalVerts; i++) {
    for (int j = 0; j < numVertPotentials; j++) {
      vertPotentials[i * numVertPotentials + j] = ((float) rand()) / RAND_MAX;
    }
  }

  for (int i = 0; i < numLocalEdges; i++) {
    for (int j = 0; j < numEdgePotentials; j++) {
      edgePotentials[i * numEdgePotentials + j] = ((float) rand()) / RAND_MAX;
    }
  }

  // Initialize Neighbors
  InitNeighbors();
}

/**
 *
 */
void Grid::SetGlobalInfo(const std::string fname) {
  int err, ncid;
  int type, gRows, gCols, gVertices, gEdges, nvPotentials, nePotentials;

  // open a file
  err = ncmpi_open(MPI_COMM_WORLD, fname.c_str(), NC_CLOBBER | NC_64BIT_DATA,
      MPI_INFO_NULL, &ncid);
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

  numGlobalVerts = numGlobalRows * numGlobalCols;
  numGlobalEdges = numGlobalRows * (numGlobalCols - 1)
    + numGlobalCols * (numGlobalRows - 1);

#ifdef DEBUG
  logFile << "\n-------------- GET GLOBAL INFO --------------\n";
  logFile << "numRows:            " << gRows << std::endl;
  logFile << "numCols:            " << gCols << std::endl;
  logFile << "numVertices:        " << gVertices << std::endl;
  logFile << "numEdges:           " << gEdges << std::endl;
  logFile << "numVertPotentials:  " << nvPotentials << std::endl;
  logFile << "numEdgePotentials:  " << nePotentials << std::endl;
  logFile << "\n";
  logFile << "numGlobalRows:      " << numGlobalRows << std::endl;
  logFile << "numGlobalCols:      " << numGlobalCols << std::endl;
  logFile << "numGlobalVertices:  " << numGlobalVerts << std::endl;
  logFile << "numGlobalEdges:     " << numGlobalEdges << std::endl;
  logFile << "numVertPotentials:  " << numVertPotentials << std::endl;
  logFile << "numEdgePotentials:  " << numEdgePotentials << std::endl;
#endif
}

void Grid::LoadExample()
{
  numGlobalRows = 4;
  numGlobalCols = 3;
  numVertPotentials = 2;
  numEdgePotentials = 4;

  numGlobalVerts = numGlobalRows * numGlobalCols;
  numGlobalEdges = numGlobalRows * (numGlobalCols - 1)
    + numGlobalCols * (numGlobalRows - 1);

  numLocalEdges = numGlobalEdges / numProcs;
  if (rank < numGlobalEdges%numProcs) {
    numLocalEdges++;
  }

  edgeList = new int[2*numLocalEdges];
  if (rank == 0) {
    numLocalEdges = 4;
    edgeList[0] = 0; edgeList[1] = 1;
    edgeList[2] = 1; edgeList[3] = 2;
    edgeList[4] = 0; edgeList[5] = 3;
    edgeList[6] = 1; edgeList[7] = 4;
  }
  if (rank == 1) {
    numLocalEdges = 4;
    edgeList[0] = 2; edgeList[1] = 5;
    edgeList[2] = 3; edgeList[3] = 4;
    edgeList[4] = 4; edgeList[5] = 5;
    edgeList[6] = 3; edgeList[7] = 6;
  }

  if (rank == 2) {
    numLocalEdges = 4;
    edgeList[0] = 4; edgeList[1] = 7;
    edgeList[2] = 5; edgeList[3] = 8;
    edgeList[4] = 6; edgeList[5] = 7;
    edgeList[6] = 7; edgeList[7] = 8;
  }

  if (rank == 3) {
    numLocalEdges = 5;
    edgeList[0] = 6; edgeList[1] = 9;
    edgeList[2] = 7; edgeList[3] = 10;
    edgeList[4] = 8; edgeList[5] = 11;
    edgeList[6] = 9; edgeList[7] = 10;
    edgeList[8] = 10; edgeList[9] = 11;
  }

  std::vector<std::pair<int,int> > offlenList = GetVerts();
  logFile << "***************************************\n";
  logFile << "***************************************\n";
  logFile << "numLocalEdges/numLocalVerts: " << numLocalEdges << ", " << numLocalVerts << std::endl;
  logFile << "numVert/EdgePotentials: " << numVertPotentials << ", " 
          << numEdgePotentials << std::endl << std::endl; 
  vertPotentials = new float[(vertexEnd - vertexStart + 1) * numVertPotentials];
  edgePotentials = new float[numLocalEdges * numEdgePotentials];
  
  for (int i = 0; i < (vertexEnd-vertexStart+1); i++) {
    for (int j = 0; j < numVertPotentials; j++) {
      vertPotentials[i*numVertPotentials + j] = -(i + vertexStart);
    }
  }

  for (int i = 0; i < numLocalEdges; i++) {
    for (int j = 0; j < numEdgePotentials; j++) {
      edgePotentials[i*numEdgePotentials + j] = -(i + 4*rank);
    }
  }
  
  edgeVal.resize(numLocalEdges);
  for (int i = 0; i < numLocalEdges; i++)
    edgeVal[i] = new Edge(numVertPotentials);

}

std::vector<std::pair<int,int> > Grid::GetVerts() {
  std::vector<std::pair<int, int> > offlenList;
  int *verts = new int[numLocalEdges << 1];
  memcpy(verts, edgeList, (numLocalEdges << 1) * sizeof(int));
  std::sort(verts, verts + (numLocalEdges<<1));
  vertexStart = verts[0];
  vertexEnd = verts[(numLocalEdges << 1) - 1];
  vSelfVec.push_back(verts[0]);
  vertCnt.push_back(1);

  int idx = 0;
  int len = 1; 
 
  int vidx = 0;
  for (int i = 1; i < (numLocalEdges << 1); i++) {
    if (verts[i] == verts[idx]) {
      vertCnt[vidx]++; 
    } else {
      if (verts[i] == verts[idx] + 1) {
        len++;
      } else {
        offlenList.push_back(std::make_pair<int, int> (verts[idx], len));
        len = 1;
      }
      vSelfVec.push_back(verts[i]);
      vertCnt.push_back(1);
      idx = i;
      vidx++;
    }
  }

  offlenList.push_back(std::make_pair<int, int> (verts[idx], len));
  numLocalVerts = vSelfVec.size();
  for (int i = 0; i < numLocalVerts; i++) {
    vSelfMapG2L[vSelfVec[i]] = i; 
  }
  
  // constrcuct vertex offset and length list
  delete verts;
  return offlenList;
}

void Grid::LoadByRow(const std::string fname) {
  int err, ncid, vid;
  err = ncmpi_open(MPI_COMM_WORLD, fname.c_str(), NC_CLOBBER | NC_64BIT_DATA,
      MPI_INFO_NULL, &ncid);
  long long starts[NDIMS];
  long long counts[NDIMS];
  assert(NDIMS > 1);
  // load the vertPotentials
  starts[0] = bbox.getUpperLeftRow() * numGlobalCols;
  starts[1] = 0;
  counts[0] = numLocalVerts;
  counts[1] = numVertPotentials;

  if (starts[0] == numGlobalVerts) {
    starts[0] = 0;
    counts[0] = 0;
    counts[1] = 0;
  }

  std::string varname = "verticePotentials";

#ifdef DEBUG
  logFile << "\n-------------- LOAD DATA --------------\n";
  logFile << varname << "\n\tstart(" << starts[0] << ", " << starts[1]
    << "), " << " end(" << counts[0] << ", " << counts[1] << ") "
    << std::endl;
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
  logFile << varname << "\n\tstart(" << starts[0] << ", " << starts[1]
    << "), " << " end(" << counts[0] << ", " << counts[1] << ") "
    << std::endl;
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
    starts[0] = numGlobalRows - 1;
  }

  // load the edgeList
  varname = "edgeList";
#ifdef DEBUG 
  logFile << varname << " with ghostEdge " << numLocalEdgesGhost;
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

// remodified by partition by edge
void Grid::InitTree() {
  /**
   * This is not a good practice since the memory is allocated in the middle of the problem
   * need have finalizeTree come with it
   */
  
  // assign potentials to edge vector
  // ROW PARTITION
  int lvid1, lvid2;
  
  for (int i = 0; i < numLocalEdges; i++) {
    lvid1 = vSelfMapG2L[edgeList[2*i]];
    lvid2 = vSelfMapG2L[edgeList[2*i+1]];
#if 0
    logFile << "lr: " << lvid1 << ", " << lvid2 << "\t";
    logFile << "counter= " << vertCnt[lvid1] << ", " << vertCnt[lvid2] << std::endl;
#endif
    // need normalization
    for (int j = 0; j < numVertPotentials; j++) {
      edgeVal[i]->x_potential_node[j] = 
        vertexMuArr[lvid1][j] / (rhoEdge * vertCnt[lvid1]);
      edgeVal[i]->x_potential_node[numVertPotentials + j] = 
        vertexMuArr[lvid2][j] / (rhoEdge * vertCnt[lvid2]);
      edgeVal[i]->x_mu_node[j] = 1.0/(float)numVertPotentials;
      edgeVal[i]->edge_mu_node[j] = 1.0/(float)numVertPotentials;
#if 0      
      if (i == 2) {
        logFile << "x_potential_node=" << edgeVal[i]->x_potential_node[j] << std::endl; 
        logFile << "numerator= " << vertexMuArr[lvid1][j] *rhoEdge << std::endl; 
        logFile << "denumerator= " << rhoEdge * vertCnt[lvid2] << std::endl;
        
        logFile << "x_potential_node=" << edgeVal[i]->x_potential_node[numVertPotentials + j] << std::endl; 
        logFile << "numerator= " << vertexMuArr[lvid2][j] *rhoEdge << std::endl; 
        logFile << "denumerator= " << rhoEdge * vertCnt[lvid2] << std::endl;
      }
#endif
    }
#if 0
    logFile << std::endl;
#endif
    for (int j = 0; j < numEdgePotentials; j++) {
      edgeVal[i]->x_potential_edge[j] = 
        edgePotentials[numEdgePotentials * i + j] / rhoEdge;
      edgeVal[i]->x_mu_edge[j] = 1.0/(float)numEdgePotentials;
    }
  }
}

/**
 * convert a binary file to netcdf file
 */
void Grid::Convert(const std::string fname)
{
  std::string line;
  std::string infoFile;
  std::string edgelistFile;
  std::string nodePoFile;
  std::string edgePoFile;

  std::ifstream inFile(fname.c_str());
  if (inFile.is_open()) {
    if (inFile.good()) {
      getline(inFile, infoFile);
      logFile << "info file:\t\t" << infoFile << std::endl;
      getline(inFile, edgelistFile);
      logFile << "edge list file:\t\t" << edgelistFile << std::endl;
      getline(inFile, nodePoFile);
      logFile << "node potential file:\t" << nodePoFile << std::endl;
      getline(inFile, edgePoFile);
      logFile << "edge potential file:\t" << edgePoFile << std::endl;
    }
  }
  inFile.close();

  inFile.open(infoFile.c_str());
  inFile >> numLocalVerts;
  inFile >> numLocalEdges;
  inFile >> numVertPotentials;
  inFile.close();

  numEdgePotentials = numVertPotentials * numVertPotentials;
  edgeList = new int[numLocalEdges * 2];
  memset(edgeList, 0, numLocalEdges * 2 * sizeof(int));
  vertPotentials = new float[numLocalVerts * numVertPotentials];
  memset(vertPotentials, 0, numLocalVerts * numVertPotentials * sizeof(float));
  edgePotentials = new float[numLocalEdges * numEdgePotentials];
  memset(edgePotentials, 0, numLocalEdges * numEdgePotentials * sizeof(float));
  
  numGlobalVerts = numLocalVerts;
  numGlobalEdges = numLocalEdges;
  return;
  inFile.open(edgelistFile.c_str(), std::ios::in | std::ios::binary);
  for (int i = 0; i < 2*numLocalEdges; i++) {
    inFile.read((char*)&edgeList[i++], sizeof edgeList[i++]);
  }
  inFile.close();

  inFile.open(nodePoFile.c_str(), std::ios::in | std::ios::binary);
  int k = 0;
  for (int i = 0; i < numLocalVerts * numVertPotentials; i++) {
    inFile.read((char*)&vertPotentials[k++], sizeof vertPotentials[k++]);
  }
  inFile.close();

  inFile.open(edgePoFile.c_str(), std::ios::in | std::ios::binary);
  k = 0;
  for (int i = 0; i < numLocalEdges * numEdgePotentials; i++) {
    inFile.read((char*)&edgePotentials[k++], sizeof edgePotentials[k++]);
  }
  inFile.close();
}
/**
 * store the graph
 */
void Grid::Save(const std::string fname) {
  MPI_Bcast(&numGlobalVerts, 1, MPI_INT, 0, MPI_COMM_WORLD);  
  MPI_Bcast(&numGlobalEdges, 1, MPI_INT, 0, MPI_COMM_WORLD);  
  MPI_Bcast(&numVertPotentials, 1, MPI_INT, 0, MPI_COMM_WORLD);  
  MPI_Bcast(&numEdgePotentials, 1, MPI_INT, 0, MPI_COMM_WORLD);  
  MPI_Barrier(MPI_COMM_WORLD);
  logFile << "---------------------------\n";
  logFile << "numLocalVerts=" << numLocalVerts << std::endl;
  logFile << "numLocalEdges=" << numLocalEdges << std::endl;
  logFile << "numVertPotentials=" << numVertPotentials << std::endl;
  logFile << "numEdgePotentials=" << numEdgePotentials << std::endl;
  logFile << "numGlobalVerts=" << numGlobalVerts << std::endl;
  logFile << "numGlobalEdges=" << numGlobalEdges << std::endl;
  logFile << bbox << std::endl;
  logFile << "---------------------------\n";
  int err, ncid;
  err = ncmpi_create(MPI_COMM_WORLD, fname.c_str(), NC_NOWRITE, MPI_INFO_NULL, &ncid);
  HANDLE_ERROR
  int nodeDims[2], edgeDims[2], edgeListDims[2];
  ncmpi_def_dim(ncid, "numVertices", numGlobalVerts, &nodeDims[0]);
  ncmpi_def_dim(ncid, "numVerticePotentials", numVertPotentials,
      &nodeDims[1]);
  ncmpi_def_dim(ncid, "numEdges", numGlobalEdges, &edgeDims[0]);
  ncmpi_def_dim(ncid, "numEdgePotentials", numEdgePotentials, &edgeDims[1]);
  ncmpi_def_dim(ncid, "numEndPointsPerEdge", 2, &edgeListDims[1]);
  edgeListDims[0] = edgeDims[0];

  int vidNode, vidEdge, vidEdgeList;
  std::string varname = "edgeList";
  ncmpi_def_var(ncid, varname.c_str(), NC_INT, NDIMS, edgeListDims,
      &vidEdgeList);
  varname = "edgePotentials";
  ncmpi_def_var(ncid, varname.c_str(), NC_FLOAT, NDIMS, edgeDims, &vidEdge);
  varname = "verticePotentials";
  ncmpi_def_var(ncid, varname.c_str(), NC_FLOAT, NDIMS, nodeDims, &vidNode);

  ncmpi_put_att_text(ncid, NC_GLOBAL, "graphType", 5, "Grid");
  ncmpi_put_att_int(ncid, NC_GLOBAL, "numRows", NC_INT, 1, &numGlobalRows);
  ncmpi_put_att_int(ncid, NC_GLOBAL, "numCols", NC_INT, 1, &numGlobalCols);
  ncmpi_put_att_int(ncid, NC_GLOBAL, "numEdges", NC_INT, 1, &numGlobalEdges);
  ncmpi_put_att_int(ncid, NC_GLOBAL, "numVertices", NC_INT, 1,
      &numGlobalVerts);
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
  err = ncmpi_put_vara_float_all(ncid, vidNode, starts, counts,
      vertPotentials);
  HANDLE_ERROR

  starts[0] = bbox.getUpperLeft().first * (2 * numGlobalCols - 1);
  starts[1] = 0;
  counts[0] = numLocalEdges;
  counts[1] = numEdgePotentials;
  err = ncmpi_put_vara_float_all(ncid, vidEdge, starts, counts,
      edgePotentials);
  HANDLE_ERROR

  err = ncmpi_close(ncid);
  HANDLE_ERROR
}

/**
 *
 */
int Grid::GetProcForGVertByRow(int gvid) {
  int nrows = (numGlobalRows + numProcs - 1) / numProcs;
  int vidPivot = nrows * (numGlobalRows % numProcs) * numGlobalCols;
  if (gvid < vidPivot) {
    return gvid / numGlobalCols / nrows;
  } else {
    return (gvid - vidPivot) / numGlobalCols / (numGlobalRows / numProcs)
      + (numGlobalRows % numProcs);
  }
}

void Grid::Partition()
{
  switch (part) {
    case ROW:
      PartitionByRow();
      break;
    case COL:
      std::cerr << "no implementation for COL\n";
      break;
    case BLOCK:
      std::cerr << "no implementation for BLOCK\n";
      break;
    default:
      break;
  }

}

void Grid::InitNeighbors() {
  switch (part) {
    case ROW:
      InitNeighborsByRow();
      break;
    case COL:
      std::cerr << "no implementation for COL\n";
      break;
    case BLOCK:
      std::cerr << "no implementation for BLOCK\n";
      break;
    case EDGE:
      InitNeighborsByEdge();
      break;
    default:
      break;
  }
}

/**
 *
 */
void Grid::AddToNeighborMap(int id,
    std::map<int, std::set<int> > &neighorMap) {
  int pid, vidSelf, vidOther = edgeList[id];
  // starting vertex
  if (id & 0x1 == 0) {
    vidSelf = edgeList[id - 1];
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

void Grid::DisplayNeighbors() {
  std::set<int>::iterator siter;
  std::map<int, std::set<int> >::iterator miter;
  std::cout << "DisplayOutNeighbors " << std::endl;
  for (miter = outNeighborMap.begin(); miter != outNeighborMap.end();
      miter++) {
    std::cout << miter->first << ": ";
    for (siter = miter->second.begin(); siter != miter->second.end();
        siter++) {
      std::cout << (*siter) << ", ";
    }
    std::cout << std::endl;
  }

  std::cout << "DisplayInNeighbors " << std::endl;
  for (miter = inNeighborMap.begin(); miter != inNeighborMap.end(); miter++) {
    std::cout << miter->first << ": ";
    for (siter = miter->second.begin(); siter != miter->second.end();
        siter++) {
      std::cout << (*siter) << ", ";
    }
    std::cout << std::endl;
  }

  std::cout << "DisplayInternalNeighbors " << std::endl;
  for (miter = internalNeighborMap.begin();
      miter != internalNeighborMap.end(); miter++) {

    std::cout << miter->first << ": ";
    for (siter = miter->second.begin(); siter != miter->second.end();
        siter++) {
      std::cout << (*siter) << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << "EndofNeighbors " << std::endl;
}

void Grid::DisplaySelf() {
  logFile << "-----------SELF NODE/EDGE Potentials----------\n";
  logFile << "self vID:\t";
  for (int i = 0; i < vSelfVec.size(); i++)
    logFile << "(" << vSelfVec[i] << ", " << vertCnt[i] <<")\t";
  logFile << std::endl;

  logFile << "self elist:\t";
  for (int i = 0 ; i < numLocalEdges; i++)
    logFile << "(" << edgeList[2 * i] << ", " << edgeList[2 * i + 1] << ")\t"; 
  logFile << std::endl;

  logFile << "self G2L:\t";
  std::map<int, int>::const_iterator miter;
  for (miter = vSelfMapG2L.begin(); miter != vSelfMapG2L.end(); miter++)
    logFile << "(" << miter->first << ", " << miter->second << ")\t"; 
  logFile << std::endl;

  logFile << "vertex potential:\n";
  for (int i = 0; i < numLocalVerts; i++) {
    for (int j = 0; j < numVertPotentials; j++) {
      logFile << vertPotentials[j + i*numVertPotentials] << " "; 
    }
    logFile << std::endl;
  }

  logFile << "edge potential:\n";
  for (int i = 0; i < numLocalEdges; i++) {
    for (int j = 0; j < numEdgePotentials; j++) {
      logFile << edgePotentials[j + i *numEdgePotentials] << " "; 
    }
    logFile << std::endl;
  }
}

/**
 *
 */
void Grid::InitNeighborsByRow() {
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
    } else if (vSelfVec[numLocalVerts - 1] < edgeList[i]) {
      AddToNeighborMap(i, outNeighborMap);
    }

    if (i % 2 == 0) {

      if (edgeList[i] <= vSelfVec[numLocalVerts - 1]
          && edgeList[i] >= vSelfVec[0])
        internalNeighborMap[edgeList[i]].insert(edgeList[i + 1]);

      if (edgeList[i + 1] <= vSelfVec[numLocalVerts - 1]
          && edgeList[i + 1] >= vSelfVec[0])
        internalNeighborMap[edgeList[i + 1]].insert(edgeList[i]);
    }
  }
}

void Grid::VerifyNeighbors() {
  std::set<int>::iterator siter;
  std::map<int, std::set<int> >::iterator miter;
  int prevVid;
  for (miter = outNeighborMap.begin(); miter != outNeighborMap.end();
      miter++) {
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
  if (rank == me)
  {
    std::set<int>::iterator siter;
    std::map<int, std::set<int> >::iterator miter;
    std::cout << "my neighbor : \n";
    for (miter = outNeighborMap.begin(); miter != outNeighborMap.end(); miter++)
    {
      std::cout << miter->first << ", ";
      for (siter = miter->second.begin(); siter != miter->second.end(); siter++)
        std::cout << *siter << " ";
      std::cout << std::endl;
    }
  }
#endif
}

void Grid::DisplayMuArr()
{
  logFile << "------------------MuArr----------------\n";
  for (int i = 0; i < vertexMuArr.size(); i++) {
    logFile << i << ": ";
    for (int j = 0; j < numVertPotentials; j++) {
      logFile << vertexMuArr[i][j] << " ";
    }
    logFile << "\n";
  }
}
/**
 *
 */
void Grid::Compute() {
  //logFile << "numVertPotentials=" << numVertPotentials << std::endl;
  logFile <<" before compute.....\n";
  DisplayEdge();
  int lvid, rvid;
  // reset vertexMuArr b/c we will use it to hold the partial sum for shared nodes
  ResetVertexMuArr();
#if 1
  for (int i = 0; i < edgeVal.size(); i++) {
    lvid = vSelfMapG2L[edgeList[i*2]];
    rvid = vSelfMapG2L[edgeList[i*2+1]];
    // update the lambda
    edgeVal[i]->UpdateLambda(penalty);
    // adjust the potential
    //if (rank == 0)
    edgeVal[i]->AdjustPotentials(edgeOp, penalty, rhoEdge);
    // sum product
    edgeVal[i]->OptimizeEdge(edgeOp, aux, vertexMuArr, lvid, rvid, penaltyPrime);
  }
  logFile <<" after compute.....\n";
  //DisplayEdge();
#endif
}

void Grid::FinalizeOptimization() {
  FinalizeCommunication();
}

void Grid::DisplayEdge()
{
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
  //std::cout << etiosflags(ios::fixed) << setprecision(2);
  //std::cout.unsetf(std::ios::floatfield);
  //std::cout.precision(4);
  logFile << std::fixed << std::setprecision(2); 
  logFile << "--------------Each Edge Tree----------------\n";
  for (int i = 0; i < numLocalEdges; i++) {
    logFile << "edge \t1: (mu)"; 
    for (int j = 0; j < numVertPotentials; j++) {
      logFile << edgeVal[i]->x_mu_node[j] << " ";
    }
    logFile<<", (po)";
    for (int j = 0; j < numVertPotentials; j++) {
      logFile << edgeVal[i]->x_potential_node[j] << " ";
    }
    logFile<<", (lambda)";
    for (int j = 0; j < numVertPotentials; j++) {
      logFile << edgeVal[i]->lambda_node[j] << " ";
    }
    logFile<<", (edge_mu)";
    for (int j = 0; j < numVertPotentials; j++) {
      logFile << edgeVal[i]->edge_mu_node[j] << " ";
    }
    logFile << "\n\t2: (mu)"; 
    
    for (int j = 0; j < numVertPotentials; j++) {
      logFile << edgeVal[i]->x_mu_node[j + numVertPotentials] << " ";
    }
    logFile<<", (po)";
    for (int j = 0; j < numVertPotentials; j++) {
      logFile << edgeVal[i]->x_potential_node[j + numVertPotentials] << " ";
    }
    logFile<<", (lambda)";
    for (int j = 0; j < numVertPotentials; j++) {
      logFile << edgeVal[i]->lambda_node[j + numVertPotentials] << " ";
    }
    logFile<<", (edge_mu)";
    for (int j = 0; j < numVertPotentials; j++) {
      logFile << edgeVal[i]->edge_mu_node[j + numVertPotentials] << " ";
    }
    logFile << "\n\te: (mu)"; 
    for (int j = 0; j < numEdgePotentials; j++) {
      logFile << edgeVal[i]->x_mu_edge[j] << " ";
    }
    logFile<<", (po)";
    for (int j = 0; j < numEdgePotentials; j++) {
      logFile << edgeVal[i]->x_potential_edge[j] << " ";
    }
    logFile << std::endl << std::endl;
  }
}

void Grid::InitOptimization() {
  logFile << "------------InitOptimization-------------\n";
  InitNeighbors();
  DisplaySelf();
  InitCommunication();
  ResetVertexMuArr();
  for (int i = 0; i < numLocalEdges; i++) {
    int lvid = vSelfMapG2L[edgeList[2*i]];
    int rvid = vSelfMapG2L[edgeList[2*i+1]];
    for (int j = 0; j < numVertPotentials; j++) {
      vertexMuArr[lvid][j] += vertPotentials[lvid*numVertPotentials + j];
      vertexMuArr[rvid][j] += vertPotentials[rvid*numVertPotentials + j];
    }
  }
  Communicate();
  UpdateMuArrAfterComm(); 
  InitTree();
  //DisplayEdge();
}

void Grid::InitCommunication() {
  //logFile << "k=" << edgeVal[0]->k << ", ksquare=" << edgeVal[0]->ksquare << std::endl;
  aux = new APGMem(edgeVal[0]->k);
  edgeOp = new EdgeOp(edgeVal[0]->k);
 
  logFile << "Init vertexMuArr\n";
  vertexMuArr.resize(numLocalVerts);
  for (int i = 0; i < numLocalVerts; i++) {
    vertexMuArr[i] = new float[numVertPotentials];
    memset(vertexMuArr[i], 0, sizeof(float) * numVertPotentials);
  }
  
  int totalRequests = outNeighborMapEdge.size() * 2;
  if (totalRequests == 0) {
#if 0
    logFile << "rank  " << rank << " has no neighbor!\n";
#endif
    return;
  }
  commBuffer = new float*[totalRequests];
  for (int i = 0; i < neighbors.size(); i++) {
    commBuffer[2*i] = new float[outNeighborMapEdge[i].size() * numVertPotentials];
    commBuffer[2*i+1] = new float[outNeighborMapEdge[i].size() * numVertPotentials];
  }
  request = (MPI_Request *) malloc(sizeof(MPI_Request) * totalRequests);
}

void Grid::FinalizeCommunication() {
  int totalRequests = outNeighborMap.size() * 2;
  if (totalRequests == 0)
    return;
  for (int i = 0; i < totalRequests; i++) {
    if (commBuffer[i]) {
      delete commBuffer[i]; 
    }
  }
  if (commBuffer) {
    delete commBuffer;
  }
  if (request)
    delete request;
}

void Grid::Communicate() {
  switch (part) {
    case ROW:
      CommunicateByRow();
      break;
    case COL:
      std::cerr << "no implementation for COL\n";
      break;
    case BLOCK:
      std::cerr << "no implementation for BLOCK\n";
      break;
    case EDGE:
      CommunicateByEdge();
      break;
    default:
      break;
  }
}

void Grid::CommunicateByRow() {
  if (numLocalVerts == 0)
    return;

  int idxSend = 0, idxRecv = 0, sidx = 0, vsize;
  // send the lower boudnary for row partition
  int edgeIdx = numLocalEdges - numGlobalCols;

  std::set<int>::iterator siter;
  std::map<int, std::set<int> >::const_iterator citer;

  for (citer = outNeighborMap.begin(); citer != outNeighborMap.end();
      citer++) {
    sidx = 0;
    vsize = citer->second.size() * numVertPotentials;
#if 0
    std::cout << "send rank: " << rank << "-->" << citer->first << ", "
      << vsize << ", " << idxSend << ", " << edgeVal.size() << std::endl;
#endif
    commBuffer[idxSend] = new float[vsize];
    for (siter = citer->second.begin(); siter != citer->second.end();
        siter++) {
      memcpy(commBuffer[idxSend] + (sidx++) * numVertPotentials,
          edgeVal[edgeIdx++]->edge_mu_node, numVertPotentials);
    }
#if 1
    MPI_Isend(commBuffer[idxSend], vsize, MPI_FLOAT, citer->first, idxSend,
        MPI_COMM_WORLD, request + idxSend);
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
    commBuffer[idxSend + idxRecv] = new float[vsize];
    for (siter = citer->second.begin(); siter != citer->second.end();
        siter++) {
      try {
        if (sidx == citer->second.size() - 1)
          break;
        memcpy(edgeVal[++edgeIdx]->edge_mu_node,
            commBuffer[idxSend + idxRecv] + (sidx++) * numVertPotentials,
            numVertPotentials);
      } catch (std::exception &e) {
        std::cout << e.what() << std::endl;
      }
    }
    //the last node
    memcpy(edgeVal[edgeIdx]->edge_mu_node + numVertPotentials,
        commBuffer[idxSend + idxRecv] + (sidx) * numVertPotentials,
        numVertPotentials);
#if 1
    MPI_Irecv(commBuffer[idxRecv + idxSend], vsize, MPI_FLOAT, citer->first,
        idxRecv, MPI_COMM_WORLD, request + idxSend + idxRecv);
#endif
    idxRecv++;
  }

  MPI_Status status;
#if 1
  for (int i = 0; i < idxSend + idxRecv; i++) {
    MPI_Wait(request + i, &status);
  }
  for (int i = 0; i < idxSend + idxRecv; i++)
    if (commBuffer[i])
      delete commBuffer[i];
#endif
}

void Grid::CommunicateByEdge() {
  int requestIdx = 0;
  int cpySize = sizeof(float) * numVertPotentials;
  int numNeighbors = neighbors.size(); 
  // vertexMuArr is the updated copy of node potential after compute  
  for (int i = 0; i < numNeighbors; i++) {
    int payloadSize = outNeighborMapEdge[i].size() * numVertPotentials;
    int nid = neighbors[i];
    for (int j = 0; j < outNeighborMapEdge[i].size(); j++) {
      memcpy(commBuffer[i] + j * numVertPotentials,
          vertexMuArr[vSelfMapG2L[outNeighborMapEdge[i][j]]], cpySize);
    }
       
    MPI_Isend(commBuffer[requestIdx], payloadSize, MPI_FLOAT, nid,
        rank, MPI_COMM_WORLD, request + requestIdx);
    requestIdx++;
    MPI_Irecv(commBuffer[requestIdx], payloadSize, MPI_FLOAT, nid, 
        nid, MPI_COMM_WORLD, request + requestIdx);
    requestIdx++;
  }

  MPI_Status status;
  for (int i = 0; i < requestIdx; i++) {
    MPI_Wait(request + i, &status);
  }
}

double Grid::ComputeLinearObj()
{
  double obj = 0.0;
  int lvid, rvid;
  for (int i = 0; i < numLocalEdges; i++) {
    lvid = vSelfMapG2L[edgeList[2*i]];
    rvid = vSelfMapG2L[edgeList[2*i+1]];
    for (int j = 0; j < numVertPotentials * 2; j++) {
      obj += edgeVal[i]->x_potential_node[j] * edgeVal[i]->x_mu_node[j]; 
    }
    for(int j = 0; j < numEdgePotentials; j++) {
      obj += edgeVal[i]->x_potential_edge[j] * edgeVal[i]->x_mu_edge[j]; 
    }
  }
  return obj*rhoEdge;
}
/**
 * exchange the boundary vertcies' potentials
 * in order to get convergence.
 */
std::ostream&
operator<<(std::ostream& strm, const Grid &g) {
  return strm << "G (" << g.numGlobalRows << ", " << g.numGlobalCols << ") "
    << "(" << g.numGlobalVerts << ", " << g.numGlobalEdges << ")";
}
