//-----------------------------------------------------------------------------
// File          : TriggerSample.h
// Author        : Ryan Herbst  <rherbst@slac.stanford.edu>
// Created       : 08/26/2011
// Project       : Heavy Photon API
//-----------------------------------------------------------------------------
// Description :
// Sample Container
// Sample Data consists of the following: Z[xx:xx] = Zeros, O[xx:xx] = Ones
//    Sample[0] = O[0], Error[0], Hybrid[1:0], Drop[0], ApvChip[2:0], Z[0], Channel[6:0], FpgaAddress[15:0]
//    Sample[1] = Z[1:0], Sample1[13:0]], Z[1:0], Sample0[13:0]
//    Sample[2] = Z[1:0], Sample3[13:0]], Z[1:0], Sample2[13:0]
//    Sample[3] = Z[1:0], Sample5[13:0]], Z[1:0], Sample4[13:0]
//-----------------------------------------------------------------------------
// Copyright (c) 2011 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 08/26/2011: created
//-----------------------------------------------------------------------------
#ifndef __TRIGGER_SAMPLE_H__
#define __TRIGGER_SAMPLE_H__
#include "TrackerSample.h"
#include <sys/types.h>
using namespace std;

//! Trigger Event Container Class
class TriggerSample : public TrackerSample {

   public:

   static const uint kSampleSize = 4;
   static const uint kEventCode = 1;
   
      //! Constructor for static pointer
      TriggerSample ();

      //! Constructor with copy
      TriggerSample ( uint *data );

      /*!
       * Returns hybrid index.
      */
      uint hybrid ( );

      //! Get error flag
      /*!
       * Returns apv error flag
      */
      bool error ( );

      bool head();

      bool tail();

      //! Get apv index.
      /*!
       * Returns apv index.
      */
      uint apv ( );

      //! Get channel index.
      /*!
       * Returns channel index.
      */
      uint channel ( );

      //! Get Front End Board Address
      /*!
       * Returns febAddress
      */
      uint febAddress ( );

      uint rceAddress();

      uint noSync();

      //! Get adc value at index.
      /*!
       * Returns adc value value.
       * \param index sub-sample index, 0-5.
      */
      uint value ( uint index );

};

#endif
