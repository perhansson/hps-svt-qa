//-----------------------------------------------------------------------------
// File          : TriggerEvent.cpp
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
#include "TriggerEvent.h"
#include "TrackerEvent.h"
using namespace std;

// Constructor
TriggerEvent::TriggerEvent () : TrackerEvent(1, 4) {}

// Deconstructor
TriggerEvent::~TriggerEvent () {}

void TriggerEvent::update() {}
