//-----------------------------------------------------------------------------
// File          : TriggerSample.cpp
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
#include "TriggerSample.h"
using namespace std;

// Constructor for static pointer
TriggerSample::TriggerSample () : TrackerSample(kSampleSize) {}

// Constructor with copy
TriggerSample::TriggerSample ( uint *data ) : TrackerSample(kSampleSize, data) {}


bool TriggerSample::head() {
   return ((data_[3]>>30)&0x1);
}

bool TriggerSample::tail() {
   return ((data_[3]>>29)&0x1);
}

// Get error flag
bool TriggerSample::error ( ) {
   return ((data_[3]>>28)&0x1);
}

// Get hybrid index.
uint TriggerSample::hybrid ( ) {
   return ((data_[3]>>26)&0x3);
}

// Get apv index.
uint TriggerSample::apv ( ) {
   return ((data_[3]>>23)&0x7);
}

// Get channel index.
uint TriggerSample::channel ( ) {
   return ((data_[3]>>16)&0x7F);
}

// Get FpgaAddress value from header.
uint TriggerSample::febAddress ( ) {
   return ((data_[3] >> 8) & 0xFF);
}

uint TriggerSample::rceAddress () {
   return (data_[3] & 0xFF);
}

uint TriggerSample::noSync () {
   return ((data_[3] >> 31) & 0x1);
}

// Get adc value at index.
uint TriggerSample::value ( uint index ) {
   switch(index) {
      case 0: return(data_[0]&0xFFFF);
      case 1: return((data_[0]>>16)&0xFFFF);
      case 2: return(data_[1]&0xFFFF);
      case 3: return((data_[1]>>16)&0xFFFF);
      case 4: return(data_[2]&0xFFFF);
      case 5: return((data_[2]>>16)&0xFFFF);
      default: return(0);
   }
}

