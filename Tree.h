/*
 * Edge.h
 *
 *  Created on: Feb 26, 2013
 *      Author: xiaocai
 */

#ifndef EDGE_H_
#define EDGE_H_

#define NUM_NODE 2
#define NUM_EDGE 1
#define EPSILON 0.00000000000000001


class APGMem;
class Edge {
public:
	int k;
	int ksquare;
	int nedges;
	float rho;
	// node
	float *x_mu_node;
	float *x_potential_node;
	float *edge_mu_node;

	//edge
	float *x_mu_edge;
	float *x_potential_edge;

	// lambda
	float *lambda_node;

	float *x_node;
	float *x_edge;
	APGMem *aux;

	float penalty;
public:
	Edge(int kval);
	virtual ~Edge();
	void setPenalty(float p) {penalty = p;}
	void setRho(float arho) {rho = arho;}
	void setNedge(float anedges) {nedges = anedges;}
	void updateLambda();
	void adjustPotentials();
	void optimizeEdge();
	void sumProductEdge();
};

class APGMem
{
public:
	float *nodeDerivative;
	float *edgeDerivative;
	float *nodePotential;
	float *edgePotential;
	float *msg_12;
	float *msg_21;
	float *bufferSquare;
	float *buffer;
	int k;
	int ksquare;
public:
	APGMem(int kval);
	~APGMem();
	void updateMu(float *node, float *edge);
};



#endif /* TREE_H_ */
