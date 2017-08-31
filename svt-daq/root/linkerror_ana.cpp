//-----------------------------------------------------------------------------
// File          : linkerror_ana.cpp
// Author        : Per Hansson Adrian  <phansson@slac.stanford.edu>
// Created       : 09/15/2015
// Project       : Kpix Software Package
//-----------------------------------------------------------------------------
// Description :
// File to look at link errors in data stream
//-----------------------------------------------------------------------------
// Copyright (c) 2009 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 09/15/2015: created
//-----------------------------------------------------------------------------

#include <iostream>
#include <sstream>
#include <cstdio>
#include "DataReadEvio.h"
#include "TiTriggerEvent.h"
#include "TriggerSample.h"


using namespace std;

TiTriggerEvent triggerEvent;
TriggerSample* triggerSample;
int eventCount;
int sampleCount;
int numEvents = -1;
bool isHeadMultiSample;
int rce;
int hybrid;
int apv;
int channel;
int feb;
int debug;
int errorCountAll;
int errorCountHead;

int main(int argc, char**argv ) {

  printf("JUST GO\n");
  
  if( argc < 2 ) {
    cout  << "need an evio file as argument." << endl;
    return 1;
  }

  if( argc > 2 ) {    
    cout << "set nr of events to " << argv[2] << endl;
    istringstream iss;
    iss.str(argv[2]);
    iss >> numEvents;
    cout << "set nr of events to " << numEvents << endl;
  }

  if( argc > 3 ) {       
    cout << "set debug level to " << argv[3] << endl;
    istringstream iss;
    iss.str(argv[3]);
    iss >> debug;
    cout << "set debug level to " << debug << endl;
  }


  DataReadEvio* dataRead = new DataReadEvio();
  dataRead->set_engrun(true);
  dataRead->set_bank_num(51);

  if( ! dataRead->open( argv[1] ) ) {
    cout << "couldn't open file " << argv[1] << endl;
    return 2;
  }
  
  cout << "Opened evio file \"" << argv[1] << "\"" << endl; 
  

  eventCount = 0;
  errorCountAll = 0;
  errorCountHead = 0;


  cout << "read events" << endl;
  

  triggerSample = new TriggerSample();

  int tiEventNumber;
  int tiEventNumberPrev = -1;
  int tiEvents = 0;
  
  while( dataRead->next(&triggerEvent) == true ) {
    
    if( debug > 0 || (tiEvents % 10000 == 0) ) cout << "read event " << eventCount << " tiEvents " << tiEvents << endl;

    if( numEvents > 0 && eventCount > numEvents )
      break;
    
    tiEventNumber = triggerEvent.tiEventNumber();
    if(tiEventNumber > tiEventNumberPrev) {
      tiEventNumberPrev = tiEventNumber;
      tiEvents++;
    }
    
    sampleCount = triggerEvent.count();
    
    if( debug > 0) cout << "sampleCount: " << sampleCount << endl;
    
    for( int iSampleCount=0; iSampleCount!=sampleCount; ++iSampleCount) {
      
      if( debug > 0) cout << "iSampleCount: " << iSampleCount << endl;
      
      triggerEvent.sample(iSampleCount, triggerSample);

      rce = triggerSample->rceAddress();
      feb = triggerSample->febAddress();
      hybrid = triggerSample->hybrid();
      apv = triggerSample->apv();
      channel = triggerSample->channel();
      


      if( debug > 0) 
        cout << "rce: " << rce << " feb: " << feb << " hybrid: " << hybrid << " apv: " << apv << " channel: " << channel << endl;
      
      
      if( debug > 0) {
        cout << "HEAD " << (triggerSample->head() ? 1 : 0) << " TAIL " << (triggerSample->tail() ? 1 : 0) << " ERR " << (triggerSample->error() ? 1 : 0);
        cout << " ADC samples: ";
        for(int y = 0; y < 6; ++y) {
          uint val = triggerSample->value( y );
          cout << " " << val;
        } 
        cout << endl;
      }

      
      // Count the error bits
      
      if( triggerSample->error() ) {
        cout << "error found: " << rce << " feb: " << feb << " hybrid: " << hybrid << " apv: " << apv << " channel: " << channel << endl;
        errorCountAll++;
      }
      
      if( triggerSample->head() ) {
        if( triggerSample->error() ) {
          cout << "head error found: " << rce << " feb: " << feb << " hybrid: " << hybrid << " apv: " << apv << " channel: " << channel << endl;
          errorCountHead++;
        }
      }
      
      
      
      
      
    } // iSampleCount
    
    
    eventCount++;
  } // while data read OK
  
  
  printf("Clean data read\n");
  dataRead->close();
  delete dataRead;
  delete triggerSample;
  
  cout << "File " << argv[1] << " tiEvents " << tiEvents << " errorCountAll " << errorCountAll << " errorCountHead: " << errorCountHead << endl;

  return 0;


} // main
