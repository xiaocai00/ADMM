/*
 * BBox.h
 *
 *  Created on: Feb 2, 2013
 *      Author: xiaocai
 */

#ifndef BBOX_H_
#define BBOX_H_
#define BBOXSIZE 4
#include <iostream>
class BBox {
  private:
    std::pair<int, int> upperLeft;
    std::pair<int, int> boxSize;

  public:
    BBox();
    BBox(int rowStart, int colStart, int rowSize, int colSize);
    virtual ~BBox();
    std::pair<int, int>& getUpperLeft();
    std::pair<int, int>& getBoxSize();
    int getUpperLeftRow();
    int getUpperLeftCol();
    int getBoxRowSize();
    int getBoxColSize();

    void setUpperLeft(int rowStart, int colStart);
    void setBoxSize(int rowSize, int colSize);
    void setBBox(int rowStart, int colStart, int rowSize, int colSize);
    void setBBox(int box[]);
  private:
    friend std::ostream& operator<<(std::ostream& strm, const BBox &b);
};

#endif /* BBOX_H_ */
