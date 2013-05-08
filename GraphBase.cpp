/*
 * GraphBase.cpp
 *
 *  Created on: Feb 1, 2013
 *      Author: xiaocai
 */

#include "GraphBase.h"
#include <iostream>

GraphBase::GraphBase(int verticePotentialSize, int edgePotentialSize)
: vertPotentials(NULL), edgePotentials(NULL)
{
  numVertPotentials = verticePotentialSize;
  numEdgePotentials = edgePotentialSize;
}

GraphBase::~GraphBase() {
}

void GraphBase::Generate(){
}

void GraphBase::LoadByRow(const std::string fname){
}

void GraphBase::Save(const std::string fname){
}

