/*
 * GraphBase.h
 *
 *  Created on: Feb 1, 2013
 *      Author: xiaocai
 */

#ifndef GRAPHBASE_H_
#define GRAPHBASE_H_
#include <iostream>
class GraphBase {
  protected:
    int numVertPotentials;
    int numEdgePotentials;
    float *vertPotentials;
    float *edgePotentials;
  public:
    GraphBase():vertPotentials(NULL), edgePotentials(NULL), numVertPotentials(0), numEdgePotentials(0){}
    GraphBase(int verticePotentialSize, int edgePotentialSize);
    virtual ~GraphBase();
    virtual void Generate();
    virtual void Save(const std::string fname);
    virtual void LoadByRow(const std::string fname);
};

#endif /* GRAPHBASE_H_ */
