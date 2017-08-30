//-----------------------------------------------------------------------------
// File          : Filter.h
// Author        : Ryan Herbst  <rherbst@slac.stanford.edu>
// Created       : 04/14/2012
// Project       : Heavy Photon API
//-----------------------------------------------------------------------------
// Description :
// Filter data container.
//-----------------------------------------------------------------------------
// Copyright (c) 2012 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 04/14/2012: created
//-----------------------------------------------------------------------------
#ifndef __FILTER_H__
#define __FILTER_H__
// Tracker Filter Container Class
class Filter {
   public:

      // Constants
      static const unsigned int IdLength    = 200;
      static const unsigned int FpgaCount   = 7;
      static const unsigned int HybridCount = 3;
      static const unsigned int ApvCount    = 5;
      static const unsigned int CoefCount   = 10;

      // Sizes
      static const unsigned int FpgaSize = HybridCount * ApvCount * CoefCount;
      static const unsigned int Size     = FpgaCount * FpgaSize;

      // ID Container
      char filterId[IdLength];

      // Filter container
      double filterData[FpgaCount][HybridCount][ApvCount][CoefCount];
};

#endif
