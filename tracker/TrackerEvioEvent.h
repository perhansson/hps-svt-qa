//-----------------------------------------------------------------------------
// File          : TrackerEvioEvent.h
// Author        : Matt Graham <mgraham@slac.stanford.edu>
// Created       : 04/05/2012
// Project       : Heavy Photon API
//-----------------------------------------------------------------------------
// Description :
// Evio Event data container.
//-----------------------------------------------------------------------------
// Copyright (c) 2011 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 4/05/2012: created
//-----------------------------------------------------------------------------
#ifndef __TRACKER_EVIO_EVENT_H__
#define __TRACKER_EVIO_EVENT_H__

#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <sys/types.h>
#include "TrackerBank.h"
using namespace std;

//! Tracker Event Container Class for an EVIO event with 1 bank-per-FPGA
class TrackerEvioEvent  {

  int nbanks;
  // Internal bank contrainer
  TrackerBank* banks_[8];
  
  // Update
  void update();
  
 public:
  
  //! Constructor
  TrackerEvioEvent ();
  
  //! Deconstructor
  ~TrackerEvioEvent ();
  
  //! Get sample count
  /*!
   * Returns sample count
   */
  int count ( ){return nbanks;};
  
  void addFPGAData(TrackerBank*);
  void restart();
  
  //! Get sample at index
  /*!
   * Returns pointer to static sample object without memory allocation.
   * Contents of returned object will change next time sample() is called.
   * \param index Sample index. 0 - count()-1.
   */
  TrackerBank *getFPGAData (uint index);
  
  //! Get sample at index
  /*!
   * Returns pointer to copy of sample object. A newly allocated sample object
   * is created and must be deleted after use.
   * \param index Sample index. 0 - count()-1.
   */
  //  TrackerBank *getFPGACopy (uint index);
  
};

#endif
