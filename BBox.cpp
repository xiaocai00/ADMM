/*
 * BBox.cpp
 *
 *  Created on: Feb 2, 2013
 *      Author: xiaocai
 */

#include "BBox.h"

BBox::BBox()
{}

BBox::BBox(int rowStart, int colStart,
    int rowSize, int colSize) {
  upperLeft = std::pair<int, int>(rowStart, colStart);
  boxSize = std::pair<int, int>(rowSize, colSize);
}

BBox::~BBox() {
}

std::pair<int, int>& BBox::getUpperLeft()
{
  return upperLeft;
}

int BBox::getUpperLeftRow()
{
	return upperLeft.first;
}

int BBox::getUpperLeftCol()
{
	return upperLeft.second;
}

std::pair<int, int>& BBox::getBoxSize()
{
  return boxSize;
}

int BBox::getBoxRowSize()
{
  return boxSize.first;
}

int BBox::getBoxColSize()
{
  return boxSize.second;
}
void BBox::setUpperLeft(int rowStart, int colStart)
{
  upperLeft = std::pair<int, int>(rowStart, colStart);
}

void BBox::setBoxSize(int rowSize, int colSize)
{
  boxSize = std::pair<int, int>(rowSize, colSize);
}

void BBox::setBBox(int rowStart, int colStart,
    int rowSize, int colSize)
{
  setUpperLeft(rowStart, colStart);
  setBoxSize(rowSize, colSize);
}

void BBox::setBBox(int box[])
{
  try {
    setUpperLeft(box[0], box[1]);
    setBoxSize(box[2], box[3]);
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
  }
}

std::ostream& operator<<(std::ostream& strm, const BBox &b)
{
  return strm << "B: (" 
              << b.upperLeft.first<< "," << b.upperLeft.second<<") "
              << ", (" 
              << b.boxSize.first << "," << b.boxSize.second<<") ";
}
