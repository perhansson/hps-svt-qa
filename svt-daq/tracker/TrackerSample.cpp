//-----------------------------------------------------------------------------
// File          : TrackerSample.cpp
// Author        : Ryan Herbst  <rherbst@slac.stanford.edu>
// Created       : 08/26/2011
// Project       : Heavy Photon API
//-----------------------------------------------------------------------------
// Description :
// Sample Container
//-----------------------------------------------------------------------------
// Copyright (c) 2011 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 08/26/2011: created
//-----------------------------------------------------------------------------
#include <string.h>
#include <stdlib.h>
#include "TrackerSample.h"
using namespace std;

// Constructor for static pointer
TrackerSample::TrackerSample (uint size) {
   size_ = size;
   ldata_ = (uint *)malloc(sizeof(uint) * size);
   data_ = ldata_;
}

// Constructor with copy
TrackerSample::TrackerSample (uint size, uint *data ) {
   size_ = size;
   ldata_ = (uint *)malloc(sizeof(uint) * size);
   memcpy(ldata_, data, sizeof(uint)*size);
   data_ = ldata_;
}

TrackerSample::~TrackerSample () {
   free(ldata_);
}

// Set data pointer.
void TrackerSample::setData ( uint *data ) {
   data_ = data;
}

uint TrackerSample::size() {
   return size_;
}


