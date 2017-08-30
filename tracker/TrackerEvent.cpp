//-----------------------------------------------------------------------------
// File          : TrackerEvent.cpp
// Author        : Ryan Herbst  <rherbst@slac.stanford.edu>
// Created       : 08/26/2011
// Project       : Heavy Photon API
//-----------------------------------------------------------------------------
// Description :
// Event Container
//-----------------------------------------------------------------------------
// Copyright (c) 2011 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 08/26/2011: created
// 02/14/2012: Updates to match FPGA. Added hooks for future TI frames.
//-----------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "TrackerEvent.h"
#include "TrackerSample.h"
using namespace std;

void TrackerEvent::update() { }

// Constructor
TrackerEvent::TrackerEvent (uint eventCode, uint sampleSize) : Data() {
   eventCode_ = eventCode;
   sampleSize_ = sampleSize;
}

// Deconstructor
TrackerEvent::~TrackerEvent () {}

// Get FpgaAddress value from header.
uint TrackerEvent::dataEventCode ( ) {
   return((data_[0] >> 24) & 0xFF);
}

bool TrackerEvent::eventCodeMatch() {
   return (dataEventCode() == eventCode_);
}

// Get sequence count from header.
uint TrackerEvent::sequence ( ) {
   return(data_[0] & 0xFFFFFF);
}

// Get sample count
uint TrackerEvent::count ( ) {
   if (eventCodeMatch()) {
      return((size_- (kHeadSize + kTailSize) ) / sampleSize_);
   } else {
      return 0;
   }
}

// Get sample at index
void TrackerEvent::sample (uint index, TrackerSample* sample) {
   if ( index < count() ) {
      sample->setData(&(data_[kHeadSize + (index*sampleSize_)]));
   }
}

// Get sample at index
// TrackerSample *TrackerEvent::sampleCopy (uint index) {
//    TrackerSample *tmp;

//    if ( index >= count() ) return(NULL);
//    else {
//       tmp = new TrackerSample (&(data_[kHeadSize+(index*sampleSize_)]));
//       return(tmp);
//    }
// }

