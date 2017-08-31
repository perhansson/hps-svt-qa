//-----------------------------------------------------------------------------
// File          : TrackerEvioEvent.cpp
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
#include "TrackerEvioEvent.h"
using namespace std;

void TrackerEvioEvent::update() { }

// Constructor
TrackerEvioEvent::TrackerEvioEvent ()  {
  nbanks=0;
}

TrackerEvioEvent::~TrackerEvioEvent ()  {
}

void TrackerEvioEvent::restart ()  {
  nbanks=0;
}

void TrackerEvioEvent::addFPGAData (TrackerBank* newfpga)  {
  if(nbanks>7) cout<<"Why am I trying to add more than 8 banks???"<<endl;
  banks_[nbanks]=newfpga;
  nbanks++;
}

// Get sample at index
TrackerBank *TrackerEvioEvent::getFPGAData (uint index) {
  return banks_[index];
}

// Get sample at index
//TrackerBank *TrackerEvioEvent::getFPGACopy (uint index) {
//   TrackerBank *tmp;
//   tmp = new TrackerBank (banks_[index]);
//   return(tmp);
//}

