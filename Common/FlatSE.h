/*
 * This file is part of libTIM.
 *
 * Copyright (©) 2005-2013  Benoit Naegel
 * Copyright (©) 2013 Theo de Carpentier
 * Copyright (©) 2022-2023  Cyril Meyer
 */

#ifndef FlatSE_h
#define FlatSE_h

#include <iostream>
#include <limits>
#include <vector>

#include "Image.h"
#include "Point.h"

namespace LibTIM {
/// Container base class for flat structuring elements (or binary masks)
/**
    WARNING: some algorithms require a connexity rather than a structuring
element in parameters. To this end, use for example make2DN8() to compute a
8-neighborhood (now the center is not included in the structuring element)
**/

class FlatSE {
 protected:
  // The flat SE is stored as a vector of points and offsets
  std::vector<Point<TCoord> > points;
  std::vector<TOffset> offsets;

  // Negative and positive offsets of SE
  // Useful to add borders
  TCoord negativeOffsets[3];
  TCoord positiveOffsets[3];

 public:
  FlatSE() {}
  FlatSE(const Image<U8> &im);

  FlatSE &operator=(const FlatSE &se);
  FlatSE(const FlatSE &se) { operator=(se); }

  /// returns the number of points contained in the structuring element
  /// (cardinal of the set)
  size_t getNbPoints() const;

  /// computes the offset of each point, according to "size"
  void setContext(const TSize *size);

  void setNegPosOffsets();

  Point<TCoord> getPoint(int i) const { return points[i]; }
  void addPoint(Point<TCoord> p) { points.push_back(p); }

  /// returns the offset of "point"
  TOffset getOffset(int point) { return offsets[point]; }

  const TCoord *getNegativeOffsets() const;
  const TCoord *getPositiveOffsets() const;

  void makeSymmetric();

  Image<U8> toImage();

  FlatSE &operator+=(FlatSE &b) {
    for (unsigned int i = 0; i < b.points.size(); i++)
      points.push_back(b.points[i]);
    return *this;
  }

  // Iterators

  typedef std::vector<Point<TCoord> >::iterator iterator_point;
  typedef std::vector<TOffset>::iterator iterator_offset;
  typedef iterator_offset iterator;

  iterator begin() { return offsets.begin(); }
  iterator end() { return offsets.end(); }

  iterator_point begin_point() { return points.begin(); }
  iterator_point end_point() { return points.end(); }

  // Various neighborhoods

  /// Basic neighborhoods in 2D N4 and N8
  void make2DN4();
  void make2DN5();
  void make2DN8();
  void make2DN9();

  void make2DEuclidianBall(int r);

  /// In 3D N6,7,18,19,26,27
  void make3DN6();
  void make3DN7();
  void make3DN18();
  void make3DN19();
  void make3DN26();
  void make3DN27();

  void make3DAxialSegment(int l);

  // Methods to create various structuring elements

  template <class VoxelType>
  void makeBallEuclidian2D(Image<VoxelType> &img, double r);
  template <class VoxelType>
  void makeBallChessboard2D(Image<VoxelType> &img, double rx, double ry);

  template <class VoxelType>
  void makeBallEuclidian3D(Image<VoxelType> &img, double r);

  template <class VoxelType>
  void makeCircle2D(Image<VoxelType> &img, double r, double t);

  void print() {
    std::cout << "FlatSE \n";
    for (unsigned int i = 0; i < points.size(); i++) points[i].print();
    std::cout << " ";
    std::cout << "\n";
  }

  void reserve(size_t size) {
    points.reserve(size);
    offsets.reserve(size);
  }
  void clear() {
    points.clear();
    offsets.clear();
  }
};

}  // namespace LibTIM

#include "FlatSE.hxx"

#endif
