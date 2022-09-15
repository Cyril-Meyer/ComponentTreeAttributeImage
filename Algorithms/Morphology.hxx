/*
 * This file is part of libTIM.
 *
 * Copyright (©) 2005-2013  Benoit Naegel
 * Copyright (©) 2013 Theo de Carpentier
 * Copyright (©) 2022-2023  Cyril Meyer
 */

#include <algorithm>
#include <cmath>
#include <queue>

#include "Common/FlatSE.h"
#include "Common/Image.h"

namespace LibTIM {

template <class T>
void addBorders(Image<T> &im, const TCoord *preWidth, const TCoord *postWidth,
                T value) {
  TSize newSize[3];
  const TSize *oriSize = im.getSize();

  for (int i = 0; i < 3; i++) {
    newSize[i] = oriSize[i] + preWidth[i] + postWidth[i];
  }

  Image<T> temp(newSize);

  typename Image<T>::iterator it;
  typename Image<T>::iterator end = temp.end();

  std::fill(temp.begin(), end, value);
  temp.copy(im, preWidth[0], preWidth[1], preWidth[2]);

  im = temp;
}

template <class T>
void addBorders(Image<T> &im, FlatSE &se, T value) {
  TSize newSize[3];
  const TSize *oriSize = im.getSize();

  TCoord backOffsets[3];
  TCoord frontOffsets[3];

  FlatSE::iterator_point itSe;
  FlatSE::iterator_point endSe = se.end_point();

  for (int i = 0; i < 3; i++) backOffsets[i] = 0;
  for (itSe = se.begin_point(); itSe != endSe; ++itSe) {
    backOffsets[0] = std::min(backOffsets[0], (itSe->x));
    backOffsets[1] = std::min(backOffsets[1], (itSe->y));
    backOffsets[2] = std::min(backOffsets[2], (itSe->z));
  }

  for (int i = 0; i < 3; i++) frontOffsets[i] = 0;
  for (itSe = se.begin_point(); itSe != endSe; ++itSe) {
    frontOffsets[0] = std::max(frontOffsets[0], (itSe->x));
    frontOffsets[1] = std::max(frontOffsets[1], (itSe->y));
    frontOffsets[2] = std::max(frontOffsets[2], (itSe->z));
  }

  for (int i = 0; i < 3; i++) {
    newSize[i] =
        oriSize[i] + (TSize)abs(backOffsets[i]) + (TSize)abs(frontOffsets[i]);
  }

  Image<T> temp(newSize);

  typename Image<T>::iterator it;
  typename Image<T>::iterator end = temp.end();

  std::fill(temp.begin(), end, value);
  temp.copy(im, abs(backOffsets[0]), abs(backOffsets[1]), abs(backOffsets[2]));

  im = temp;
}

template <class T>
Image<T> dilation(Image<T> im, FlatSE se) {
  Image<T> res = im;

  // Symmetric of structuring element, according to Heijman's
  // definition of dilation.
  se.makeSymmetric();

  const TCoord *back = se.getNegativeOffsets();
  const TCoord *front = se.getPositiveOffsets();

  T minValue = std::numeric_limits<T>::min();
  addBorders(im, back, front, minValue);
  se.setContext(im.getSize());

  typename Image<T>::iteratorXYZ it;
  typename Image<T>::iteratorXYZ end = res.end();

  FlatSE::iterator itSe;
  FlatSE::iterator endSe = se.end();

  T min = std::numeric_limits<T>::min();

  T currentMax;

  for (it = res.begin(); it != end; ++it) {
    currentMax = min;
    TOffset offsetRes =
        im.getOffset(it.x + back[0], it.y + back[1], it.z + back[2]);
    for (itSe = se.begin(); itSe != endSe; ++itSe) {
      TOffset offsetIm = offsetRes + *itSe;
      currentMax = std::max(currentMax, im(offsetIm));
    }
    *it = currentMax;
  }

  return res;
}

template <class T>
Image<T> erosion(Image<T> im, FlatSE se) {
  Image<T> res = im;

  const TCoord *back = se.getNegativeOffsets();
  const TCoord *front = se.getPositiveOffsets();

  T maxValue = std::numeric_limits<T>::max();
  addBorders(im, back, front, maxValue);
  se.setContext(im.getSize());

  typename Image<T>::iteratorXYZ it;
  typename Image<T>::iteratorXYZ end = res.end();

  FlatSE::iterator itSe;
  FlatSE::iterator endSe = se.end();

  T max = std::numeric_limits<T>::max();

  T currentMin;

  for (it = res.begin(); it != end; ++it) {
    currentMin = max;
    TOffset offsetRes =
        im.getOffset(it.x + back[0], it.y + back[1], it.z + back[2]);
    for (itSe = se.begin(); itSe != endSe; ++itSe) {
      TOffset offsetIm = offsetRes + *itSe;
      currentMin = std::min(currentMin, im(offsetIm));
    }
    *it = currentMin;
  }

  return res;
}

template <class T>
Image<T> opening(Image<T> im, FlatSE se) {
  return dilation(erosion(im, se), se);
}

template <class T>
Image<T> closing(Image<T> im, FlatSE se) {
  return erosion(dilation(im, se), se);
}

template <class T>
Image<T> morphologicalGradient(Image<T> im, FlatSE se) {
  Image<T> tmp = erosion(im, se);
  im = dilation(im, se);
  im -= tmp;
  return im;
}

template <class T>
Image<T> internalMorphologicalGradient(Image<T> im, FlatSE se) {
  Image<T> tmp = erosion(im, se);
  im -= tmp;
  return im;
}

template <class T>
Image<T> externalMorphologicalGradient(Image<T> im, FlatSE se) {
  Image<T> tmp = dilation(im, se);
  tmp -= im;
  return tmp;
}

}  // namespace LibTIM
