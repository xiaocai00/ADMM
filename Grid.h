/*
 * GraphGenerator.h
 *
 *  Created on: Feb 1, 2013
 *      Author: xiaocai
 */

#ifndef GRAPHGENERATOR_H_
#define GRAPHGENERATOR_H_
#include <fstream>
#include <iostream>
#include <exception>
#include <vector>
#include <map>
#include <set>
#include "mpi.h"
#include "GraphBase.h"
#include "BBox.h"
#include "Tree.h"

enum Partitioner{ROW, COL, BLOCK, EDGE};

class Grid: public GraphBase {
  private:
    // Global Graph
    int numLocalRows;
    int numLocalCols;
    int numGlobalRows;
    int numGlobalCols;
    APGMem *aux;
    EdgeOp *edgeOp;
    // Partitioner type
    BBox bbox;
    float rhoEdge;
    Partitioner part;

    // Local Graph
    int numGlobalVerts;
    int numGlobalEdges;
    int numLocalVerts;
    int numLocalEdges;
    int numLocalEdgesGhost;
    // input data structure
    int *edgeList;
    // inferred data structure
    std::map<int, std::set<int> > outNeighborMap;
    std::vector<std::vector<int> > outNeighborMapEdge;
    std::vector<std::vector<int> > cntNeighborMapEdge;
    std::map<int, std::set<int> > inNeighborMap;
    std::map<int, std::set<int> > internalNeighborMap;
    std::map<int, int> vSelfMapG2L;
    std::vector<std::pair<int, int> > offlenList;
    std::vector<int> vSelfVec;
    std::vector<float*> vertexMuArr;
    std::vector<int> neighbors;
    // MPI Info
    int rank;
    int numProcs;
    int vertexStart;
    int vertexEnd;
    float penalty;
    float penaltyPrime;
#if DEBUG
    std::ofstream logFile;
#endif
    std::vector<Edge *> edgeVal;
    MPI_Request *request;
    float** commBuffer;

    int GetProcForGVertByRow(int gvid);
    void ResetVertexMuArr();
    
    void InitTree();
    void VerifyNeighbors();
    void AddToNeighborMap(int vid, std::map<int, std::set<int> > &neighorMap);
    void InitNeighborsByRow();
    void CommunicateByRow();
    void CommunicateByEdge();
    void InitCommunication();
    void FinalizeCommunication();

    void InitNeighbors();
    void InitNeighborsByEdge();
    
    void DisplaySelf();
    void DisplayNeighbors();
    void DisplayNeighborsByEdge();
    
    void UpdateMuArrBeforeComm();
    void UpdateMuArrAfterComm();
    std::vector<int> vertArr;
    std::vector<int> vertCnt;
  public:
    void SetPenalties(float p, float pp) {penalty = p; penaltyPrime = pp;}
    Grid();
    Grid(int globalRowSize, int globalColSize, int verticePotentialSize,
        int edgePotentialSize);
    Grid(int localRowSize, int localColSize, int glocalRowSize, int globalColSize,
        int verticePotentialSize, int edgePotentialSize);
    
    void GetVerts();
    void SetGlobalInfo(int globalRowSize, int globalColSize, 
    int verticePotentialSize, int edgePotentialSize);
    void SetGlobalInfo(const std::string fname);
    
    void SetPartitioner(Partitioner aPart) {part = aPart;}
    void SetRhoEdge(float aRhoEdge) {rhoEdge = aRhoEdge; }
    void Convert(const std::string fname);
    
    void Partition();
    void Generate();
    void Compute();
    void Communicate();
    void Update();
    double ComputeLinearObj();
    
    void InitOptimization();
    void FinalizeOptimization();
    
    void Load(std::string fname);
    void LoadByRow(const std::string fname);
    void LoadExample();
    void Save(const std::string fname);
    
    void DisplayEdge();
    
    int GetRank() {
      return rank;
    }
    int GetNprocs() {
      return numProcs;
    }
    ~Grid();
  private:
    void DisplayMuArr();
    void PartitionByEdge(std::string fname);
    void PartitionByRow();
    friend std::ostream& operator<<(std::ostream& strm, const Grid &g);
};

#endif /* GRAPHGENERATOR_H_ */
