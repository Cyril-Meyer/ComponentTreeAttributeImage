/*
 * This file is part of libTIM.
 *
 * Copyright (©) 2005-2013  Benoit Naegel
 * Copyright (©) 2013 Theo de Carpentier
 * Copyright (©) 2022-2023  Cyril Meyer
 */

#ifndef Types_h
#define Types_h

#include <cstdint>

namespace LibTIM {
typedef unsigned char U8;
typedef signed char S8;

typedef unsigned short U16;
typedef signed short S16;

typedef unsigned long U32;
typedef signed long S32;

// Table
template <class T, int N>
struct Table {
  Table(){};
  Table(const Table& v) {
    for (int i = 0; i < N; i++) el[i] = v.el[i];
  }
  Table(int p) {
    for (int i = 0; i < N; i++) el[i] = p;
  }
  Table(int* vect) {
    for (int i = 0; i < N; i++) el[i] = vect[i];
  }

  T el[N];
  T& operator[](int i) { return el[i]; }
};

// Type of RGB point
typedef Table<U8, 3> RGB;

// Type of image size
typedef long TSize;

// Type of point spacing
typedef double TSpacing;

// Type of points coordinates
typedef long TCoord;

// Type of label
typedef unsigned long TLabel;

// Type of offset
typedef long TOffset;

const float FLOAT_EPSILON = 0.0000000001f;
}  // namespace LibTIM

#endif
