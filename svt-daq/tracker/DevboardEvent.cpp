//-----------------------------------------------------------------------------
// File          : DevboardEvent.cpp
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
#include "DevboardEvent.h"
using namespace std;

void DevboardEvent::update() { }

// Constructor
DevboardEvent::DevboardEvent () : Data() {
   double       temp;
   double       tk;
   double       res;
   double       volt;
   unsigned int idx;

   // Fill temperature lookup table
   temp = minTemp_;
   while ( temp < maxTemp_ ) {
      tk = k0_ + temp;
      //res = t25_ * exp(coeffA_+(coeffB_/tk)+(coeffC_/(tk*tk))+(coeffD_/(tk*tk*tk)));      
      res = constA_ * exp(beta_/tk);
      volt = (res*vmax_)/(rdiv_+res);
      idx = (uint)((volt / vref_) * (double)(adcCnt_-1));
      if ( idx < adcCnt_ ) tempTable_[idx] = temp; 
      temp += incTemp_;
   }

   // Fill new temperature lookup table
   temp = minTemp_;
   while ( temp < maxTemp_ ) {
      tk = k0_ + temp;
      //res = t25_ * exp(coeffA_+(coeffB_/tk)+(coeffC_/(tk*tk))+(coeffD_/(tk*tk*tk)));      
      res = constA_ * exp(beta_/tk);
      volt = (res*vmax_)/(rdiv_+res);
      idx = (uint)((volt / vrefNew_) * (double)adcCnt_);
      if ( idx < adcCnt_ ) tempTableNew_[idx] = temp; 
      temp += incTemp_;
   }
}

// Deconstructor
DevboardEvent::~DevboardEvent () {
}

// Get TI flag from header
bool DevboardEvent::isTiFrame ( ) {
   return((data_[0] & 0x80000000) != 0);
}

// Get FpgaAddress value from header.
uint DevboardEvent::fpgaAddress ( ) {
   return(data_[0] & 0xFFFF);
}

// Get sequence count from header.
uint DevboardEvent::sequence ( ) {
   return(data_[1]);
}

// Get trigger block from header.
uint * DevboardEvent::tiData ( ) {
   return(&(data_[2]));
}


// Get temperature values from header.
double DevboardEvent::temperature ( uint index, bool oldHybrid ) {
   uint adcValue;
   uint convValue;
   uint bitmask;

   if ( isTiFrame () ) return(0.0);

   if (oldHybrid) bitmask = 0xFFF;
   else bitmask = 0xFFFF;

   switch (index) {
      case  0: adcValue = (data_[2]&bitmask);
	break;
      case  1: adcValue = ((data_[2]>>16)&bitmask);
	break;
      case  2: adcValue = (data_[3]&bitmask);
	break;
      case  3: adcValue = ((data_[3]>>16)&bitmask);
	break;
      case  4: adcValue = (data_[4]&bitmask);
	break;
      case  5: adcValue = ((data_[4]>>16)&bitmask);
	break;
      case  6: adcValue = (data_[5]&bitmask);
	break;
      case  7: adcValue = ((data_[5]>>16)&bitmask);
	break;
      case  8: adcValue = (data_[6]&bitmask);
	break;
      case  9: adcValue = ((data_[6]>>16)&bitmask);
	break;
      case 10: adcValue = (data_[7]&bitmask);
	break;
      case 11: adcValue = ((data_[7]>>16)&bitmask);
	break;
      default: return(0.0);
   }

   if (oldHybrid) {
      return (tempTable_[adcValue]);
   } else {
      if ( adcValue & 0x8000 ) convValue = ((adcValue >> 3) & 0xFFF);
      else convValue = adcValue & 0xFFF;

      return (tempTableNew_[convValue]);
   }
}

// Get sample count
uint DevboardEvent::count ( ) {
   if ( isTiFrame () ) return(0);
   else return((size_-(headSize_ + tailSize_))/sampleSize_);
}

// Get sample at index
DevboardSample *DevboardEvent::sample (uint index) {
   if ( isTiFrame () ) return(NULL);
   else if ( index >= count() ) return(NULL);
   else {
      sample_.setData(&(data_[headSize_+(index*sampleSize_)]));
      return(&sample_);
   }
}

// Get sample at index
DevboardSample *DevboardEvent::sampleCopy (uint index) {
   DevboardSample *tmp;

   if ( isTiFrame () ) return(NULL);
   else if ( index >= count() ) return(NULL);
   else {
      tmp = new DevboardSample (&(data_[headSize_+(index*sampleSize_)]));
      return(tmp);
   }
}

