/*
 * Edge.h
 *
 *  Created on: Feb 26, 2013
 *      Author: xiaocai
 */

#ifndef EDGE_H_
#define EDGE_H_

#include <iostream>
#include <vector>
#define NUM_NODE 2
#define NUM_EDGE 1
#define EPSILON 0.00000000000000001

class EdgeOp;
class APGMem;

class Edge {
  public:
    int k;
    int ksquare;
    int nedges;
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

  private:
    void UpdateMu();
    void sumProductEdge(APGMem *aux, std::vector<float*>&, int, int);
  public:
    Edge(int kval);
    virtual ~Edge();
    void SetNEdge(float anedges) {nedges = anedges;}
    void UpdateLambda(float penalty);
    void AdjustPotentials(EdgeOp *edgeOp, float penalty, float rho);
    void OptimizeEdge(EdgeOp *edgeOp, APGMem *aux, 
        std::vector<float*>&, int, int, float penaltyPrime); 
};

class EdgeOp
{
public:
  int k;
  int ksquare;
  float *x_node;
  float *x_edge;
public:
  EdgeOp(int kval);
  virtual ~EdgeOp();
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
};

#endif /* TREE_H_ */
