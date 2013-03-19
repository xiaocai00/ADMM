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

enum Partitioner{ROW, COL, BLOCK};

class Grid: public GraphBase {
private:
	// Global Graph
	int numLocalRows;
	int numLocalCols;
	int numGlobalRows;
	int numGlobalCols;

	// Partitioner type
	BBox bbox;
	float rhoEdge;
	Partitioner part;

	// Local Graph
	int numGlobalVertices;
	int numGlobalEdges;
	int numLocalVerts;
	int numLocalEdges;
	int numLocalEdgesGhost;
	// input data structure
	int *edgeList;
	// inferred data structure
	std::map<int, std::set<int> > outNeighborMap;
	std::map<int, std::set<int> > inNeighborMap;
	std::map<int, std::set<int> > internalNeighborMap;
	std::map<int, int> vSelfMapG2L;
	std::vector<int> vSelfVec;

	// MPI Info
	int rank;
	int nprocs;
#if DEBUG
        std::ofstream logFile;
#endif
	std::vector<Edge *>edgeVal;
	MPI_Request *request;
        float** commBuffer;
	void InitNeighborsByRow();
	int GetProcForGVertByRow(int gvid);
	void AddToNeighborMap(int vid, std::map<int, std::set<int> > &neighorMap);
        void DisplaySelf();
        void VerifyNeighbors(); 
        void InitTree();
        void InitCommunication();
        void FinalizeCommunication();
	void InitNeighbors();
	void DisplayNeighbors();
public:
	Grid();
	Grid(int globalRowSize, int globalColSize, int verticePotentialSize,
			int edgePotentialSize);
	Grid(int localRowSize, int localColSize, int glocalRowSize, int globalColSize,
			int verticePotentialSize, int edgePotentialSize);
	void SetGlobalInfo(int globalRowSize, int globalColSize, 
                int verticePotentialSize, int edgePotentialSize);
        void SetGlobalInfo(const std::string fname);
        void SetPartitioner(Partitioner aPart) {part = aPart;}
        void setRhoEdge(float aRhoEdge) {rhoEdge = aRhoEdge; }
	void Partition();
	void Generate();
	void Save(const std::string fname);
	void Load(const std::string fname);
	void Compute();
	void Communicate();
        void InitOptimization();
        void FinalizeOptimization();
        int GetRank() {return rank;}
        int GetNprocs() {return nprocs;}
	~Grid();
private:
	friend std::ostream& operator<<(std::ostream& strm, const Grid &g);
};

#endif /* GRAPHGENERATOR_H_ */
