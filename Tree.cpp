/*
 * Edge.cpp
 *
 *  Created on: Feb 26, 2013
 *      Author: xiaocai
 */

#include <cmath>
#include <iostream>
#include <cstring>
#include <exception>
#include "Tree.h"
Edge::Edge(int kval) : k(kval), ksquare(kval*kval)
{
	x_mu_node = new float[2*k];
	x_potential_node = new float[2*k];
	edge_mu_node = new float[2*k];
	x_mu_edge = new float[ksquare * NUM_EDGE];
	x_potential_edge = new float[ksquare * NUM_EDGE];
	lambda_node = new float[2*k];

	for (int i = 0; i < k; i++) {
		x_mu_node[i] = 1.0/k;
		edge_mu_node[i] = 1.0/k;
		lambda_node[i] = 1.0/k;
	}
	for (int i = 0; i < k; i++) {
		x_mu_node[i+k] = 1.0/k;
		edge_mu_node[i+k] = 1.0/k;
		lambda_node[i+k] = 1.0/k;
	}
	
        x_node = new float[k * NUM_NODE];
	x_edge = new float[ksquare * NUM_EDGE];
	aux = new APGMem(k);
	penalty = 0.0;
        nedges = 1;
}

Edge::~Edge() {
	delete[] x_mu_node;
	delete[] x_potential_node;
	delete[] edge_mu_node;
	delete[] x_mu_edge;
	delete[] x_potential_edge;
	delete[] lambda_node;
	delete[] x_node;
	delete[] x_edge;
}

void Edge::updateLambda()
{
	int j, len;
	for (int i = 0; i < nedges; i++) {
		len = 2*k;
		for (j = 0; j < len; j++) {
			lambda_node[j] += penalty*(x_mu_node[j] - edge_mu_node[j]);
		}
	}
}

void Edge::adjustPotentials()
{
	int len = k << 1;
	for (int i = 0; i < len; i++) {
		x_node[i] = rho * (x_potential_node[i] + lambda_node[i])
						- penalty * (edge_mu_node[i] - x_mu_node[i]);
	}
	for (int i = 0; i < ksquare; i++) {
		x_edge[i] = rho * x_potential_edge[i];
	}
}

void Edge::optimizeEdge()
{
	float value;
	int j = 0;

	// update mu
	aux->updateMu(x_mu_node, x_mu_edge);

	// derivatives
	for (int i = 0; i < k; i++) {
		value = log(x_mu_node[i]);
		aux->nodeDerivative[i] = 0.0 - value - 1;
		for (int ik = 0; ik < k; ik++) {
                    try {
			aux->edgeDerivative[j] = value - log(x_mu_edge[j]) + 1;
                        j++;
                    } catch(std::exception &e) {
                      std::cerr << e.what() << ", " << j << std::endl;
                    }
		}
	}
        
#if 1
        // derivatives
	int start;
	for (int i = 0; i < k; i++) {
		value = log(x_mu_node[i+k]);
		aux->nodeDerivative[i+k] = 0.0 - value -1;
		start = i;
		for (j = start; j < ksquare; j = j+k) {
			aux->edgeDerivative[j] += value;
		}
	}

	// potentials
	for (int i = 0; i < (k<<1); i++) {
		value = x_node[i] + penalty * aux->nodeDerivative[i];
		value /= penalty;
		aux->nodePotential[i] = 0.0 - value;
	}

	for (int i = 0; i < ksquare; i++) {
		value = x_edge[i] + penalty * aux->edgeDerivative[i];
		value /= penalty;
		aux->edgePotential[i] = 0.0 - value;
	}
	sumProductEdge();
#endif
}

void Edge::sumProductEdge()
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
	}

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
	}

	sum = 0.0;
	for (int i = 0; i < k; i++) {
		aux->buffer[i] = exp(aux->nodePotential[i]) * aux->msg_21[i];
		sum += aux->buffer[i];
	}

	for (int i = 0; i < k; i++) {
		x_mu_node[i] = aux->buffer[i]/sum;
	}

	sum = 0.0;
	for (int i = 0; i < k; i++) {
		aux->buffer[i] = exp(aux->nodePotential[i+k]) * aux->msg_12[i];
		sum += aux->buffer[i];
	}

	for (int i = 0; i < k; i++) {
		x_mu_node[i+k] = aux->buffer[i]/sum;
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
				x_mu_edge[i] = x_mu_node[i] * aux->bufferSquare[j] / aux->buffer[i];
			} else {
				x_mu_edge[i] = 0.0;
			}
		}
		start = end;
	}
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

void APGMem::updateMu(float *node, float *edge)
{
	for (int i = 0; i < k * 2; i++) {
		if (node[i] <= 0.0)
			node[i] = EPSILON;
	}

	for (int i = 0; i < ksquare; i++) {
		if (edge[i] <= 0.0)
			edge[i] = EPSILON;
	}
}
