//-----------------------------------------------------------------------------
// File          : TrackerEvent.h
// Author        : Ryan Herbst  <rherbst@slac.stanford.edu>
// Created       : 08/26/2011
// Project       : Heavy Photon API
//-----------------------------------------------------------------------------
// Description :
// Event data container.
//-----------------------------------------------------------------------------
// Copyright (c) 2011 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 08/26/2011: created
// 02/14/2012: Updates to match FPGA. Added hooks for future TI frames.
//----------------------------------------------------------------------------
// Description :
// Event Container
// Event Data consists of the following: Z[xx:xx] = Zeros
//    Frame Size = 1 x 32-bits (32-bit dwords)
//    Header = 8 x 32-bits
//       Header[0] = T[0], Z[14:0], FpgaAddress[15:0] - T = 1 For TI FPGA
//       Header[1] = Sequence[31:0]
//
//    The rest of the event header depends on the T flag, For T = 0:
//
//       Header[2] = TempB[15:0], TempA[15:0] -- Hybrid 0
//       Header[3] = TempD[15:0], TempC[15:0] -- Hybrid 0
//       Header[4] = TempF[15:0], TempE[15:0] -- Hybrid 1
//       Header[5] = TempH[15:0], TempG[15:0] -- Hybrid 1
//       Header[6] = TempJ[15:0], TempI[15:0] -- Hybrid 2
//       Header[7] = TempL[15:0], TempK[15:0] -- Hybrid 2
//
//       Samples... (See TrackerSample.h)
//
//    For T = 1:
//
//       Header[2] = TBD, Waiting for clarification from JLAB
//       Header[3] = TBD, Waiting for clarification from JLAB
//       Header[4] = TBD, Waiting for clarification from JLAB
//       Header[5] = TBD, Waiting for clarification from JLAB
//       Header[6] = TBD, Waiting for clarification from JLAB
//       Header[7] = TBD, Waiting for clarification from JLAB
//
//    Tail = 1 x 32-bits
//       Should be zero
//-----------------------------------------------------------------------------
// Copyright (c) 2011 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 08/26/2011: created
//-----------------------------------------------------------------------------
#ifndef __TRACKER_EVENT_H__
#define __TRACKER_EVENT_H__

#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <sys/types.h>
#include "TrackerSample.h"
#include <Data.h>
using namespace std;

//! Tracker Event Container Class
class TrackerEvent : public Data {


 private:
  // Update
  void update();


 protected:
  // Frame Constants
  static const uint kHeadSize   = 1;
  static const uint kTailSize   = 1;
  
  uint eventCode_;
  uint sampleSize_;
  

   public:

      //! Constructor
      TrackerEvent (uint eventCode, uint sampleSize);

      //! Deconstructor
      ~TrackerEvent ();

      uint dataEventCode();

      bool eventCodeMatch();

      //! Get sequence count from header.
      /*!
       * Returns sequence count
      */
      uint sequence ( );

      //! Get sample count
      /*!
       * Returns sample count
      */
      uint count ( );

      //! Get sample at index
      /*!
       * Returns pointer to static sample object without memory allocation.
       * Contents of returned object will change next time sample() is called.
       * \param index Sample index. 0 - count()-1.
      */
      void sample (uint index, TrackerSample* sample);

      //! Get sample at index
      /*!
       * Returns pointer to copy of sample object. A newly allocated sample object
       * is created and must be deleted after use.
       * \param index Sample index. 0 - count()-1.
      */
/*       TrackerSample *sampleCopy (uint index); */

};

#endif
