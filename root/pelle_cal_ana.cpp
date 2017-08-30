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
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TFile.h>


using namespace std;

#define MAX_RCE 14
#define MAX_FEB 10
#define MAX_HYB 4
#define MAX_DELAY 9

TiTriggerEvent triggerEvent;
TriggerSample* triggerSample;
int eventCount;
int sampleCount;
int numEvents = -1;
bool isHeadMultiSample;
int rce;
int hybrid;
int apv;
int apvch;
int channel;
int feb;
int debug;
int errorCountAll;
int errorCountHead;
int delay;
ostringstream oss;
char name[200];
TH2F *baselineSamplesHist2D[MAX_FEB][MAX_HYB][6];
TH2F *samplesHist2D[MAX_FEB][MAX_HYB];
TH2F *samplesDelayHist2D[MAX_DELAY][MAX_FEB][MAX_HYB];
TH1F*sampleCountHist;
TCanvas *c1, *c2;
TString inname;
bool flip_channels = true;
bool mux_channels = false;
int chanMap[128];




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

  if (inname=="")
    {
      inname=argv[optind];
      
      inname.ReplaceAll(".bin","");
      if (inname.Contains('/')) {
        inname.Remove(0,inname.Last('/')+1);
      }
    }

  for (int idx = 0; idx < 128; idx++ ) {
     chanMap[(32*(idx%4)) + (8*(idx/4)) - (31*(idx/16))] = idx;
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


  cout << "create histograms" << endl;
    for (int fpga = 0;fpga<MAX_FEB;fpga++) {
       for (int hyb = 0;hyb<MAX_HYB;hyb++) {
          oss.str("");        
          oss << "samples-fpga-"<<fpga<<"-hyb-"<<hyb;
          samplesHist2D[fpga][hyb] = new TH2F(oss.str().c_str(),oss.str().c_str(),6,0.0,6.0,50,0.,10000.0);
          for (int sample = 0;sample<6;sample++) {
             oss.str("");        
             oss << "baselinesamples-fpga-"<<fpga<<"-hyb-"<<hyb<<"-sample-"<<sample;
             baselineSamplesHist2D[fpga][hyb][sample] = new TH2F(oss.str().c_str(),oss.str().c_str(),640,0,640,50,0., 10000.0);
          }        
          for(int idelay=0;idelay<MAX_DELAY;++idelay) {
            oss.str("");        
            oss << "samples-delay-"<< idelay <<"-fpga-"<<fpga<<"-hyb-"<<hyb;
            samplesDelayHist2D[idelay][fpga][hyb] = new TH2F(oss.str().c_str(),oss.str().c_str(),6,0.0,6.0,50,0.,10000.0);
            
          }
       }
    }

  oss.str("");        
  oss << "sampleCountHist";
  sampleCountHist = new TH1F(oss.str().c_str(),oss.str().c_str(),50,0,50);


  cout << "read events" << endl;
  

  triggerSample = new TriggerSample();

  int tiEventNumber;
  int tiEventNumberPrev = -1;
  int tiEvents = 0;
  int delay = 0;
  while( dataRead->next(&triggerEvent) == true ) {
    
    if( debug > 0 || (tiEvents % 10000 == 0) ) cout << "read event " << eventCount << " tiEvents " << tiEvents << endl;

    if( numEvents > 0 && eventCount > numEvents )
      break;
    
    tiEventNumber = triggerEvent.tiEventNumber();
    if(tiEventNumber > tiEventNumberPrev) {

       if( tiEvents % 100 == 0 && tiEvents > 0) {
          delay++;
          if( delay >= MAX_DELAY) delay = MAX_DELAY-1; 
       }
          cout << "delay " << delay << " " << tiEvents << endl;

      tiEventNumberPrev = tiEventNumber;
      tiEvents++;
    }
    
    sampleCount = triggerEvent.count();
    
    if( debug > 0) cout << "sampleCount: " << sampleCount << endl;
    
    sampleCountHist->Fill(sampleCount);




    for( int iSampleCount=0; iSampleCount!=sampleCount; ++iSampleCount) {
      
      if( debug > 0) cout << "iSampleCount: " << iSampleCount << endl;
      
      triggerEvent.sample(iSampleCount, triggerSample);

      rce = triggerSample->rceAddress();
      feb = triggerSample->febAddress();
      hybrid = triggerSample->hybrid();
      apv = triggerSample->apv();
      apvch = triggerSample->channel();
      
      if (mux_channels) channel = chanMap[apvch];
      else channel = apvch;
      
      if (flip_channels)
        channel += (4-apv)*128;
      else
        channel += apv*128;

      
      
      if( debug > 0) 
        cout << "rce: " << rce << " feb: " << feb << " hybrid: " << hybrid << " apv: " << apv << " channel: " << channel << endl;
            
      if( debug > 0)
        cout << "HEAD " << (triggerSample->head() ? 1 : 0) << " TAIL " << (triggerSample->tail() ? 1 : 0) << " ERR " << (triggerSample->error() ? 1 : 0);

      if( debug > 0) cout << " ADC samples: ";
      for(int y = 0; y < 6; ++y) {
        uint val = triggerSample->value( y );
        baselineSamplesHist2D[feb][hybrid][y]->Fill(channel,val);
        samplesHist2D[feb][hybrid]->Fill(y,val);
        samplesDelayHist2D[delay][feb][hybrid]->Fill(y,val);
        
        if( debug > 0) cout << " " << val;
      } 
      if( debug > 0) cout << endl;
      

      
    } // iSampleCount
    
    
    eventCount++;
  } // while data read OK
  
  
  printf("Clean data read\n");
  dataRead->close();
  delete dataRead;
  delete triggerSample;
  
  cout << "File " << argv[1] << " tiEvents " << tiEvents << " errorCountAll " << errorCountAll << " errorCountHead: " << errorCountHead << endl;
  




  gROOT->SetStyle("Plain");
  gStyle->SetOptStat("emrou");
  gStyle->SetPalette(1,0);
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.1);                
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetMarkerStyle(6);
  c1 = new TCanvas("c1","c1",1200,900);
  c2 = new TCanvas("c2","c2",1200,900);

  sprintf(name,"%s_cal_ana.root",inname.Data());
  TFile* outputFile = new TFile(name,"RECREATE");
  
  
  sprintf(name,"%s_baselinesamples.ps[",inname.Data());
  c2->Print(name);

  c2->Clear();
  c2->cd();
  sampleCountHist->Draw();
  sprintf(name,"%s_baselinesamples.ps",inname.Data());
  c2->Print(name);
  
  outputFile->cd();
  sampleCountHist->Write();


  for (int fpga = 0;fpga<MAX_FEB;fpga++) {
     for (int hyb = 0;hyb<MAX_HYB;hyb++) {        
        c2->Clear();            
        c2->cd();
        sprintf(name,"%s_samples_F%d_H%d.png",inname.Data(),fpga,hyb);
        printf("Drawing %s\n",name);
        samplesHist2D[fpga][hyb]->Draw("colz");
        //c2->SaveAs(name);
        sprintf(name,"%s_baselinesamples.ps",inname.Data());
        c2->Print(name);
        outputFile->cd();
        samplesHist2D[fpga][hyb]->Write();

        for (int idelay=0; idelay<MAX_DELAY; ++idelay) {
          c2->Clear();            
          c2->cd();
          sprintf(name,"%s_samples_delay_D%d_F%d_H%d.png",inname.Data(),idelay,fpga,hyb);
          printf("Drawing %s\n",name);
          samplesDelayHist2D[idelay][fpga][hyb]->Draw("colz");
          sprintf(name,"%s_baselinesamples.ps",inname.Data());
          c2->Print(name);
          outputFile->cd();
          samplesDelayHist2D[idelay][fpga][hyb]->Write();       
        }
     }
     
  }


  for (int fpga = 0;fpga<MAX_FEB;fpga++) {
    for (int hyb = 0;hyb<MAX_HYB;hyb++) {        
      for(int sample=0; sample<6; ++sample) {
        c2->Clear();            
        c2->cd();
        sprintf(name,"%s_baselinesamples_F%d_H%d_S%d.png",inname.Data(),fpga,hyb,sample);
        printf("Drawing %s\n",name);
        baselineSamplesHist2D[fpga][hyb][sample]->Draw("colz");
        //c2->SaveAs(name);
        sprintf(name,"%s_baselinesamples.ps",inname.Data());
        c2->Print(name);
        outputFile->cd();
        baselineSamplesHist2D[fpga][hyb][sample]->Write();
      }
    }
  }
  
  sprintf(name,"%s_baselinesamples.ps]",inname.Data());
  c1->Print(name);
  

  outputFile->Close();
  


  return 0;
  

} // main
