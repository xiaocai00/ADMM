/*
 * Edge.cpp
 *
 *  Created on: Feb 26, 2013
 *      Author: xiaocai
 */

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <exception>
#include "Tree.h"

Edge::Edge(int kval) : k(kval), ksquare(kval*kval){
  x_mu_node = new float[2*k];
  x_potential_node = new float[2*k];
  edge_mu_node = new float[2*k];
  lambda_node = new float[2*k];

  x_mu_edge = new float[ksquare * NUM_EDGE];
  x_potential_edge = new float[ksquare * NUM_EDGE];
  
  for (int i = 0; i < k; i++) {
    x_mu_node[i] = 1.0/k;
    edge_mu_node[i] = 1.0/k;
    lambda_node[i] = 0.0;
  }
  for (int i = 0; i < k; i++) {
    x_mu_node[i+k] = 1.0/k;
    edge_mu_node[i+k] = 1.0/k;
    lambda_node[i+k] = 0.0;
  }
  nedges = 1;
  rho = 1;
}

Edge::~Edge() {
  delete[] x_mu_node;
  delete[] x_potential_node;
  delete[] edge_mu_node;
  delete[] lambda_node;
  delete[] x_mu_edge;
  delete[] x_potential_edge;
}

void Edge::UpdateMu()
{
  for (int i = 0; i < k * 2; i++) {
    if (x_mu_node[i] <= 0.0)
      x_mu_node[i] = EPSILON;
  }

  for (int i = 0; i < ksquare; i++) {
    if (x_mu_edge[i] <= 0.0)
      x_mu_edge[i] = EPSILON;
  }
}

/**
 * ADMM step 3 update Lagrangian multiplier
 */
//float penalty;
//float penaltyPrime;
void Edge::updateLambda(float penalty)
{
  int j, len;
  for (int i = 0; i < nedges; i++) {
    len = 2*k;
    for (j = 0; j < len; j++) {
      lambda_node[j] += penalty*(x_mu_node[j] - edge_mu_node[j]);
    }
  }
}

/**
 * ADMM update node and edge potential
 */
void Edge::adjustPotentials(EdgeOp *edgeOp, float penalty)
{
  int len = k << 1;
#if 0
  std::cout << "-------------adjustPotentials BEGIN------------" << std::endl;
#endif
  for (int i = 0; i < len; i++) {
    edgeOp->x_node[i] = rho * (x_potential_node[i] + lambda_node[i])
                - penalty * (edge_mu_node[i] - x_mu_node[i]);
#if 0
    printf("%.2f, %.2f %.2f p=%.2f d=%.2f\n",
        edgeOp->x_node[i],
        x_potential_node[i],
        lambda_node[i], penalty, 
        (edge_mu_node[i] - x_mu_node[i]));
#endif
  }
#if 0
  std::cout << std::endl;
#endif
  for (int i = 0; i < ksquare; i++) {
    edgeOp->x_edge[i] = rho * x_potential_edge[i];
#if 0
    printf("%.2f ", edgeOp->x_edge[i]);
#endif
  }
#if 0
  std::cout << "\n-------------adjustPotentials END------------" << std::endl;
#endif
}

/**
 *  ADMM update the edge
 */
void Edge::optimizeEdge(EdgeOp *edgeOp, APGMem *aux, std::vector<float*>& vertexMuArr,
    int lvid, int rvid, float penaltyPrime) 
{
  // update mu
  UpdateMu();
  // derivatives
  int j = 0;
  float value;
  for (int i = 0; i < k; i++) {
    value = log(x_mu_node[i]);
    aux->nodeDerivative[i] = -1.0 - value;
    for (int ik = 0; ik < k; ik++) {
      try {
        aux->edgeDerivative[j] = value - log(x_mu_edge[j]) + 1;
#if 0
        printf("%.2f, ", aux->edgeDerivative[j]);
#endif
        j++;
      } catch (std::exception &e) {
        std::cerr << e.what() << ", " << j << std::endl;
      }
    }
  }
#if 0
  std::cout << std::endl;
#endif

  // derivatives
  int start;
  for (int i = 0; i < k; i++) {
    value = log(x_mu_node[i+k]);
    aux->nodeDerivative[i+k] = -1.0 - value ;
    start = i;
    for (j = start; j < ksquare; j += k) {
      aux->edgeDerivative[j] += value;
#if 0
      printf("%.2f, ", aux->edgeDerivative[j]);
#endif
    }
  }
#if 0
  std::cout << std::endl;
  std::cout << std::endl;
#endif
  // potentials
  for (int i = 0; i < (k<<1); i++) {
    value = edgeOp->x_node[i] + penaltyPrime * aux->nodeDerivative[i];
    value /= penaltyPrime;
    aux->nodePotential[i] = 0.0 - value;
#if 0
    printf("%.2f, ", aux->nodePotential[i]);
#endif
  }
#if 0
  std::cout << std::endl;
#endif
#if 1
  for (int i = 0; i < ksquare; i++) {
    value = edgeOp->x_edge[i] + penaltyPrime * aux->edgeDerivative[i];
    value /= penaltyPrime;
    aux->edgePotential[i] = 0.0 - value;
#if 0
    printf("%.2f ", aux->edgePotential[i]);
#endif
  }

#if 0
  std::cout << std::endl;
  std::cout << std::endl;
#endif
  sumProductEdge(aux, vertexMuArr, lvid, rvid);
#endif
}

void Edge::sumProductEdge(APGMem *aux, std::vector<float*> &vertexMuArr, int lvid, int rvid)
{
  int j, p, start, end;
  float sum = 0.0;
  memcpy(aux->bufferSquare, aux->edgePotential, ksquare*sizeof(float));
  for (int i = 0; i < k; i++) {
    p = 0;
    start = i;
    aux->buffer[i] = 0.0;
    for (j = start; j < ksquare; j += k) {
      aux->buffer[i] += exp(aux->edgePotential[j] + aux->nodePotential[p]);
      aux->bufferSquare[j] += aux->nodePotential[p];
      p++;
    }
    sum += aux->buffer[i];
  }

  for (int i = 0; i < k; i++) {
    aux->msg_12[i] = aux->buffer[i] / sum;
#if 0
    printf("%.2f ", aux->msg_12[i]);
#endif
  }
#if 0
  printf("\n");
#endif
  start = 0;
  sum = 0.0;
  for (int i = 0; i < k; i++) {
    p = k;
    end = start + k;
    aux->buffer[i] = 0.0;
    for (j = start; j < end; j++) {
      aux->buffer[i] += exp(aux->edgePotential[j] + aux->nodePotential[p]);
      aux->bufferSquare[j] += aux->nodePotential[p];
      p++;
    }
    start = end;
    sum += aux->buffer[i];
  }

  for (int i = 0; i < k; i++) {
    aux->msg_21[i] = aux->buffer[i] / sum;
#if 0
    printf("%.2f ", aux->msg_21[i]);
#endif
  }
#if 0
  printf("\n");
#endif
  sum = 0.0;
  for (int i = 0; i < k; i++) {
    aux->buffer[i] = exp(aux->nodePotential[i]) * aux->msg_21[i];
    sum += aux->buffer[i];
  }

  for (int i = 0; i < k; i++) {
    x_mu_node[i] = aux->buffer[i]/sum;
    vertexMuArr[lvid][i] += x_mu_node[i];
  }
  sum = 0.0;
  for (int i = 0; i < k; i++) {
    aux->buffer[i] = exp(aux->nodePotential[i+k]) * aux->msg_12[i];
    sum += aux->buffer[i];
  }
  for (int i = 0; i < k; i++) {
    x_mu_node[i+k] = aux->buffer[i]/sum;
    vertexMuArr[rvid][i] += x_mu_node[i+k];
  }
  
  start = 0;
  for (int i = 0; i < k; i++) {
    end = start + k;
    aux->buffer[i] = 0.0;
    for (j = start; j < end; j++) {
      aux->bufferSquare[j] = exp(aux->bufferSquare[j]);
      aux->buffer[i] += aux->bufferSquare[j];
    }

    for (j = start; j < end; j++) {
      if (aux->buffer[i] != 0.0) {
        x_mu_edge[j] = x_mu_node[i] * aux->bufferSquare[j] / aux->buffer[i];
      } else {
        x_mu_edge[j] = 0.0;
      }
    }
    start = end;
  }
}

EdgeOp::EdgeOp(int kval) : k(kval), ksquare(kval*kval)
{
  x_node = new float[k * NUM_NODE];
  x_edge = new float[ksquare * NUM_EDGE];
}

EdgeOp::~EdgeOp()
{
  delete x_node;
  delete x_edge;
}

APGMem::APGMem(int kval) : k(kval), ksquare(kval*kval)
{
  nodeDerivative = new float[NUM_NODE * k];
  edgeDerivative = new float[NUM_EDGE * ksquare];
  nodePotential = new float[NUM_NODE * k];
  edgePotential = new float[NUM_EDGE * ksquare];
  msg_12 = new float[k];
  msg_21 = new float[k];
  buffer = new float[k];
  bufferSquare = new float[ksquare];
}

APGMem::~APGMem()
{
  delete[] nodeDerivative;
  delete[] edgeDerivative;
  delete[] nodePotential;
  delete[] edgePotential;
  delete[] msg_12;
  delete[] msg_21;
  delete[] buffer;
  delete[] bufferSquare;
}

