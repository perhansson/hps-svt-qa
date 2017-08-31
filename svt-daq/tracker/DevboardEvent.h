//-----------------------------------------------------------------------------
// File          : DevboardEvent.h
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
//       Samples... (See DevboardSample.h)
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
#ifndef __DEVBOARD_EVENT_H__
#define __DEVBOARD_EVENT_H__

#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <sys/types.h>
#include "DevboardSample.h"
#include <Data.h>
using namespace std;

//! Devboard Event Container Class
class DevboardEvent : public Data {

      // Temperature Constants
      static const double coeffA_     = -1.4141963e1;
      static const double coeffB_     =  4.4307830e3;
      static const double coeffC_     = -3.4078983e4;
      static const double coeffD_     = -8.8941929e6;
      static const double beta_       = 3750;
      static const double constA_     = 0.03448533;
      static const double t25_        = 10000.0;
      static const double k0_         = 273.15;
      static const double vmax_       = 2.5;
      static const double vref_       = 2.5;
      static const double vrefNew_    = 2.048;
      static const double rdiv_       = 10000;
      static const double minTemp_    = -50;
      static const double maxTemp_    = 150;
      static const double incTemp_    = 0.01;
      static const uint   adcCnt_     = 4096;

      // Temperature lookup table
      double tempTable_[adcCnt_];
      double tempTableNew_[adcCnt_];

      // Frame Constants
      static const unsigned int headSize_   = 8;
      static const unsigned int tailSize_   = 1;
      static const unsigned int sampleSize_ = 4;

      // Internal sample contrainer
      DevboardSample sample_;

      // Update
      void update();

   public:

      //! Constructor
      DevboardEvent ();

      //! Deconstructor
      ~DevboardEvent ();

      //! Is frame TI frame?
      bool isTiFrame ( );

      //! Get FpgaAddress value from header.
      /*!
       * Returns fpgaAddress
      */
      uint fpgaAddress ( );

      //! Get sequence count from header.
      /*!
       * Returns sequence count
      */
      uint sequence ( );

      //! Get trigger block from header.
      /*!
       * Returns trigger block
      */
      uint * tiData ( );

      //! Get temperature values from header.
      /*!
       * Returns temperature value.
       * \param index temperature index, 0-11.
      */
      double temperature ( uint index, bool oldHybrid = false );

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
      DevboardSample *sample (uint index);

      //! Get sample at index
      /*!
       * Returns pointer to copy of sample object. A newly allocated sample object
       * is created and must be deleted after use.
       * \param index Sample index. 0 - count()-1.
      */
      DevboardSample *sampleCopy (uint index);

};

#endif
