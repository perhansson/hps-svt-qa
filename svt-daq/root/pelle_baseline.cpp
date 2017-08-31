//-----------------------------------------------------------------------------
// File          : cal_summary.cc
// Author        : Ryan Herbst  <rherbst@slac.stanford.edu>
// Created       : 03/03/2011
// Project       : Kpix Software Package
//-----------------------------------------------------------------------------
// Description :
// File to generate calibration summary plots.
//-----------------------------------------------------------------------------
// Copyright (c) 2009 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 03/03/2011: created
//-----------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <vector>
#include <TFile.h>
#include <TH1F.h>
#include <meeg_utils.hh>
#include <TProfile.h>
#include <TH2F.h>
#include <TF1.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
//#include <TApplication.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLegend.h>
#include <stdarg.h>
#include <DevboardEvent.h>
#include <DevboardSample.h>
#include <TiTriggerEvent.h>
#include <TriggerSample.h>
#include <Data.h>
#include <DataRead.h>
#include <DataReadEvio.h>
#include <unistd.h>
using namespace std;

#define MAX_RCE 14
#define MAX_FEB 10
#define MAX_HYB 4

//Functions
uint getTimeBin(double t);
int compareInts(const void* a, const void* b);
//const int nTimeBins = 18;
//double timeBins[nTimeBins] = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,35,100};
const int nTimeBins = 31;
//double timeBins[nTimeBins] = {0,5,10,15,20,25,30,100};
double timeBins[nTimeBins] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,40};




int main ( int argc, char **argv ) {

  printf("Just Go.\n");

  bool flip_channels = true;
  bool mux_channels = false;
  bool read_temp = true;
  int hybrid_type = 0;
  bool evio_format = false;
  bool triggerevent_format = false;
  int use_rce = -1;
  int use_fpga = -1;
  int use_hybrid = -1;
  int use_chan = -1;
  int num_events = -1;
  int c;
  bool doPulseShape = false;
  bool doTimeEvo = false;
  TCanvas         *c1, *c2, *c3;
  int chanMap[128];
  for (int idx = 0; idx < 128; idx++ ) {
    chanMap[(32*(idx%4)) + (8*(idx/4)) - (31*(idx/16))] = idx;
  }
  ostringstream oss;

  TH2F* pulseShapeHist[MAX_RCE][MAX_FEB][MAX_HYB][640];
  int hybridCount[MAX_RCE][MAX_FEB][MAX_HYB];
  int *channelCount[MAX_RCE][MAX_FEB][MAX_HYB];
  int *channelAllCount[MAX_RCE][MAX_FEB][MAX_HYB];
  TH2F *baselineHist2D[MAX_RCE][MAX_FEB][MAX_HYB];
  TH2F *baselineSamplesHist2D[MAX_RCE][MAX_FEB][MAX_HYB][6];
  TH2F *baselineSamples0TimeBinHist2D[MAX_RCE][MAX_FEB][MAX_HYB][nTimeBins];
  TH2F *maxSamplePhaseMod24Hist2D[MAX_RCE][MAX_FEB][MAX_HYB];
  TH1F* timeStampDiffHist;
  TH1F* timeStampDiffBinHist;
  TH1F* timeStampHistMod6;
  TH1F* timeStampHistMod12;
  TH1F* timeStampHistMod24;  

  DataRead        *dataRead;
  int svt_bank_num = 3;
  DevboardEvent    event;
  TiTriggerEvent    triggerevent;
  TriggerSample   *triggersample = new TriggerSample();
  int		samples[6];
  int            eventCount;
  int            tiEventCount;
  int runCount = -1;
  TString inname;
  TString outdir;
  char            name[200];
  unsigned long timeStamp;
  unsigned long tiEventNumber;
  unsigned long timeStampPrev;
  unsigned long tiEventNumberPrev;
  unsigned long timeStampDiff;
  double timeStampDiffD;
  uint tiTimeBin;
  int rce;
  int fpga;
  int samplecount;
  int hyb;
  int apv;
  int apvch;
  int channel;
  bool goodSample = true;
  int maxSample;
  double thresholdOcc = 20.0;

    
  while ((c = getopt(argc,argv,"ho:nmct:R:H:F:e:Es:b:VA:B:C:D:T")) !=-1)
    switch (c)
      {
      case 'h':
        printf("-h: print this help\n");
        printf("-o: use specified output filename\n");
        printf("-n: DAQ (Ryan's) channel numbering\n");
        printf("-m: number channels in raw mux order\n");
        printf("-c: don't compute correlations\n");
        printf("-t: hybrid type (1 for old test run hybrid, 2 for new 2014 hybrid)\n");
        printf("-R: use only specified RCE\n");
        printf("-F: use only specified FPGA\n");
        printf("-H: use only specified hybrid\n");
        printf("-C: use only specified channel were applicable\n");
        printf("-e: stop after specified number of events\n");
        printf("-E: use EVIO file format\n");
        printf("-b: EVIO bank number for SVT (default 3)\n");
        printf("-V: use TriggerEvent event format\n");
        printf("-P: Do pulse shape plots\n");    
        printf("-T: Do time b/w/ trigger evolution for chosen channels\n");
        return(0);
        break;
      case 'o':
        inname = optarg;
        outdir = optarg;
        if (outdir.Contains('/')) {
          outdir.Remove(outdir.Last('/')+1);
        }
        else outdir="";
        break;
      case 'n':
        flip_channels = false;
        break;
      case 'm':
        mux_channels = true;
        break;
      case 't':
        hybrid_type = atoi(optarg);
        break;
      case 'R':
        use_rce = atoi(optarg);
        break;
      case 'F':
        use_fpga = atoi(optarg);
        break;
      case 'H':
        use_hybrid = atoi(optarg);
        break;
      case 'C':
        use_chan = atoi(optarg);
        break;
      case 'e':
        num_events = atoi(optarg);
        break;
      case 'E':
        evio_format = true;
        break;
      case 'b':
        svt_bank_num = atoi(optarg);
        break;
      case 'V':
        triggerevent_format = true;
        break;
      case 'P':
        doPulseShape = true;
        break;
      case 'T':
        doTimeEvo = true;
        break;
      case '?':
        printf("Invalid option or missing option argument; -h to list options\n");
        return(1);
      default:
        abort();
      }
    
    
  if (hybrid_type==0) {
    printf("WARNING: no hybrid type set; use -t to specify old or new hybrid\n");
    printf("Configured for old (test run) hybrid\n");
    hybrid_type = 1;
  }
  
  if (evio_format) {
    DataReadEvio *tmpDataRead = new DataReadEvio();
    if (!triggerevent_format)
      tmpDataRead->set_bank_num(svt_bank_num);
    else
      {
        tmpDataRead->set_engrun(true);
        tmpDataRead->set_bank_num(51);
      }
    dataRead = tmpDataRead;
     cout << "use evio_format reader\n";
  } else 
    dataRead = new DataRead();

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
  c3 = new TCanvas("c3","c3",1200,900);

  // Start X11 view
  //TApplication theApp("App",NULL,NULL);

  // Root file is the first and only arg
  if ( argc-optind < 1 ) {
    cout << "Usage: pelle_baseline data_files\n";
    return(1);
  }
  
  if (inname=="")
    {
      inname=argv[optind];
      
      inname.ReplaceAll(".bin","");
      if (inname.Contains('/')) {
        inname.Remove(0,inname.Last('/')+1);
      }
    }
  

  printf("Init some histos\n");

  oss.str("");
  oss << "timeStampDiff";
  timeStampDiffHist = new TH1F(oss.str().c_str(),oss.str().c_str(),100,0,100.);
  oss.str("");
  oss << "timeStampDiffBin";
  timeStampDiffBinHist = new TH1F(oss.str().c_str(),oss.str().c_str(),nTimeBins,0,nTimeBins);
  oss.str("");
  oss << "timeStampMod6";
  timeStampHistMod6 = new TH1F(oss.str().c_str(),oss.str().c_str(),6,0,6);
  oss.str("");
  oss << "timeStampMod12";
  timeStampHistMod12 = new TH1F(oss.str().c_str(),oss.str().c_str(),12,0,12);
  oss.str("");
  oss << "timeStampMod24";
  timeStampHistMod24 = new TH1F(oss.str().c_str(),oss.str().c_str(),24,0,24);


  for (int rce = 0;rce<MAX_RCE;rce++) {
    for (int fpga = 0;fpga<MAX_FEB;fpga++) {
      for (int hyb = 0;hyb<MAX_HYB;hyb++) {

        hybridCount[rce][fpga][hyb] = 0;

        oss.str("");
        oss << "baseline-rce-"<<rce<<"-fpga-"<<fpga<<"-hyb-"<<hyb;
        baselineHist2D[rce][fpga][hyb] = new TH2F(oss.str().c_str(),oss.str().c_str(),640,0,640,10,0., 14000.0);

        
        oss.str("");
        oss << "maxSamplePhaseMod24-rce-"<<rce<<"-fpga-"<<fpga<<"-hyb-"<<hyb;
        maxSamplePhaseMod24Hist2D[rce][fpga][hyb] = new TH2F(oss.str().c_str(),oss.str().c_str(),6,0,6,24,0., 24.0);

        
        for (int sample = 0;sample<6;sample++) {
          oss.str("");        
          oss << "baselinesamples-rce-"<<rce<<"-fpga-"<<fpga<<"-hyb-"<<hyb<<"-sample-"<<sample;
          baselineSamplesHist2D[rce][fpga][hyb][sample] = new TH2F(oss.str().c_str(),oss.str().c_str(),640,0,640,10,0., 14000.0);
        }

        for (int tb = 0;tb<nTimeBins;tb++) {
          baselineSamples0TimeBinHist2D[rce][fpga][hyb][tb]=NULL;
          oss.str("");        
          oss << "baselinesample0timebin-rce-"<<rce<<"-fpga-"<<fpga<<"-hyb-"<<hyb<<"-sample0-timebin-"<<tb;

          // do hi-res for selected hybrid          
          if(use_rce!=-1 && use_rce==rce) {
            if(use_fpga!=-1 && use_fpga==fpga) {
              if(use_hybrid!=-1 && use_hybrid==hyb) {             
                baselineSamples0TimeBinHist2D[rce][fpga][hyb][tb] = new TH2F(oss.str().c_str(),oss.str().c_str(),640,0,640,100,1000., 14000.0);
              }
            }
          }
          if(baselineSamples0TimeBinHist2D[rce][fpga][hyb][tb]==NULL) {
            baselineSamples0TimeBinHist2D[rce][fpga][hyb][tb] = new TH2F(oss.str().c_str(),oss.str().c_str(),640,0,640,10,1000., 14000.0);
          }
          baselineSamples0TimeBinHist2D[rce][fpga][hyb][tb]->Sumw2();
          
        }
        
        
        for (int ch =0;ch<640;ch++) {          
          pulseShapeHist[rce][fpga][hyb][ch] = NULL;
        }        
      }
    }
  }

  int currentFile = optind;
  int nInputFilesToProcess = argc-optind+optind;
  
  while(currentFile<nInputFilesToProcess) {
      
      //cout << "Reading data file " <<argv[optind] << endl;
      cout << "Reading data file " <<argv[currentFile] << endl;
      // Attempt to open data file
      if ( ! dataRead->open(argv[currentFile]) ) return(2);
      
      
      bool readOK;
      cout << "Get first event " << endl;
      if (triggerevent_format) {
        dataRead->next(&triggerevent);
      } else {
        dataRead->next(&event);
      }
      
      eventCount = 0;
      tiEventCount = 0;
      timeStampPrev = 0;
      tiEventNumberPrev = 0;
      timeStampDiff = 0;
      timeStampDiffD = 0;
      tiTimeBin = 0;
      

  do {
    rce = 0;
    fpga = 0;

    if (triggerevent_format) {
      samplecount = triggerevent.count();
      //printf("datacode %d, sequence %d, samplecount %d\n",triggerevent.dataEventCode(),triggerevent.sequence(),samplecount);
    } else {
      fpga = event.fpgaAddress();
      samplecount = event.count();
      if (read_temp && !event.isTiFrame()) for (uint i=0;i<4;i++) {
          printf("Event %d, temperature #%d: %f\n",eventCount,i,event.temperature(i,hybrid_type==1));
          read_temp = false;
        }
    }

    if (!triggerevent_format && fpga==7) 
      {
        //printf("not a data event\n");
        continue;
      }
    //cout<<"  fpga #"<<event.fpgaAddress()<<"; number of samples = "<<event.count()<<endl;
    if (eventCount%1000==0) printf("Event %d\n",eventCount);    
    //if (num_events!=-1 && eventCount >= num_events) break;
    if (num_events!=-1 && tiEventCount >= num_events) break;

    for (int x=0; x < samplecount; x++) {
      
      goodSample = true;
      maxSample = -1;
      

      // Get sample
      if (triggerevent_format) {
        triggerevent.sample(x,triggersample);
        rce = triggersample->rceAddress();
        fpga = triggersample->febAddress();
        hyb = triggersample->hybrid();
        apv = triggersample->apv();
        apvch = triggersample->channel();
        goodSample = (!triggersample->head() && !triggersample->tail());
        for ( int y=0; y < 6; y++ ) {
          //printf("%x\n",sample->value(y));
          samples[y] = triggersample->value(y);

          //find max sample
          if(maxSample!=-1) {
            if(samples[y]>samples[maxSample]) {
              maxSample = y;
            }
          } else {
            maxSample = y;
          }
        }
      } else {
        DevboardSample *sample  = event.sample(x);
        hyb = sample->hybrid();
        apv = sample->apv();
        apvch = sample->channel();
        for ( int y=0; y < 6; y++ ) {
          //printf("%x\n",sample->value(y));
          samples[y] = sample->value(y) & 0x3FFF;
          if (samples[y]==0) goodSample = false;
          
        }
      }
      
      
      if (mux_channels) channel = chanMap[apvch];
      else channel = apvch;
      
      if (flip_channels)
        channel += (4-apv)*128;
      else
        channel += apv*128;
      
      //if (eventCount==0) printf("channel %d\n",channel);
      
      if ( channel >= (5 * 128) ) {
        cout << "Channel " << dec << channel << " out of range" << endl;
        cout << "Apv = " << dec << apv << endl;
        cout << "Chan = " << dec << apvch << endl;
        exit(1);
      }
      
      

      // calculate the time to the previous trigger
      if(triggerevent.tiEventNumber()>tiEventNumber) {
        tiEventCount++;
        // it's a new event, save the old ones
        tiEventNumberPrev = tiEventNumber;
        timeStampPrev = timeStamp;
        // update the current ones
          tiEventNumber = triggerevent.tiEventNumber();
          timeStamp = triggerevent.timeStamp();
          timeStampHistMod6->Fill(timeStamp%6);
          timeStampHistMod12->Fill(timeStamp%12);
          timeStampHistMod24->Fill(timeStamp%24);
          
          // difference in time
          timeStampDiff = timeStamp-timeStampPrev;
          timeStampDiffD = 4.0e-3*((double)timeStamp-timeStampPrev);
          
          timeStampDiffHist->Fill(timeStampDiffD);          
          //tiTimeBin = (uint)(timeStampDiffD/10.0);
          tiTimeBin = getTimeBin(timeStampDiffD);
          //if(tiTimeBin>=nTimeBins) tiTimeBin=9;
          timeStampDiffBinHist->Fill(tiTimeBin);
          
          
          if(tiEventCount%10000==0) {
            if(triggerevent_format) {
              printf("TI event count %d TI event number %ld (prev %ld)\n",tiEventCount, tiEventNumber, tiEventNumberPrev);
              printf("time stamp %ld (prev %ld diff %ld diffD %f [us] tiTimeBin %d)\n",timeStamp, timeStampPrev, timeStampDiff, timeStampDiffD, tiTimeBin);
            }
          }

          //printf("processing hybrid: rce = %d, feb = %d, hyb = %d\n",rce,fpga,hyb);        
          
          
      }
      //cout << timeStampDiff <<  " 4ns clocks ( "<< timeStampDiffD <<" us) -> timebin (25us bins) " <<  tiTimeBin << endl;
      
      
      if (use_rce!=-1 && rce!=use_rce) continue;
      if (use_fpga!=-1 && fpga!=use_fpga) continue;
      if (use_hybrid!=-1 && hyb!=use_hybrid) continue;
      if (!goodSample) continue;



      



      //printf("fpga %d, hybrid %d, apv %d, chan %d\n",fpga,hyb,apv,apvch);
      //printf("event %d\tx=%d\tR%d F%d H%d A%d channel %d, samples:\t%d\t%d\t%d\t%d\t%d\t%d\n",eventCount,x,rce,fpga,hyb,apv,apvch,samples[0],samples[1],samples[2],samples[3],samples[4],samples[5]);
      
      
      // Filter APVs
      //if ( eventCount >= 0 ) {
        if (hybridCount[rce][fpga][hyb]==0) {
          printf("found new hybrid: rce = %d, feb = %d, hyb = %d\n",rce,fpga,hyb);
          channelCount[rce][fpga][hyb] = new int[640];
          channelAllCount[rce][fpga][hyb] = new int[640];
          for (int i=0;i<640;i++) {
            channelCount[rce][fpga][hyb][i] = 0;
            channelAllCount[rce][fpga][hyb][i] = 0;
          }
        }
        hybridCount[rce][fpga][hyb]++;
        channelCount[rce][fpga][hyb][channel]++;
        

        if ( doPulseShape && pulseShapeHist[rce][fpga][hyb][channel] == NULL) {
          oss.str("");
          oss << "pulseShapeHist_" << rce << "_" << fpga << "_" << hyb << "_" << channel;
          cout << "Create histogram " << oss.str() << endl;
          pulseShapeHist[rce][fpga][hyb][channel] = new TH2F(oss.str().c_str(),oss.str().c_str(),6,0,6,50,0,10000);
        }

        for ( int y=0; y < 6; y++ ) {
          //printf("%x\n",sample->value(y));
          int value = samples[y];
          if (value<1000) {
            printf("out of range: event %d, rce = %d, feb = %d, hyb = %d, channel = %d, sample[%d] = %d\n",eventCount,rce,fpga,hyb,channel,y,samples[y]);
          }
                    
          baselineHist2D[rce][fpga][hyb]->Fill(channel,value);
          baselineSamplesHist2D[rce][fpga][hyb][y]->Fill(channel,value);

          channelAllCount[rce][fpga][hyb][channel]++;
          


          if ( doPulseShape) 
            pulseShapeHist[rce][fpga][hyb][channel]->Fill(y,value);
          
        } //y=samples


        if(maxSample>=3 && maxSample<=5) {
          baselineSamples0TimeBinHist2D[rce][fpga][hyb][tiTimeBin]->Fill(channel,samples[0]);
        }

        //printf("baselineSamples0TimeBinHist2D event %d, rce = %d, feb = %d, hyb = %d, channel = %d, sample[0] = %d timebin %d N so far: %f\n",eventCount,rce,fpga,hyb,channel,samples[0],tiTimeBin,baselineSamples0TimeBinHist2D[rce][fpga][hyb][tiTimeBin]->GetEntries());

        maxSamplePhaseMod24Hist2D[rce][fpga][hyb]->Fill(maxSample,timeStamp%24);
        
        //} // eventCount>20
    } // x=samplecount
    eventCount++;

    
    if (triggerevent_format) {
      readOK = dataRead->next(&triggerevent);
    } else {
      readOK = dataRead->next(&event);
    }
    
  } while (readOK);

  dataRead->close();
  
  if (!evio_format && eventCount != runCount)
    {
      printf("ERROR: events read = %d, runCount = %d\n",eventCount, runCount);
    }
  

  printf("Process summary [%d] :\nEventCount = %d\nTiEventCount=%d\n",currentFile,eventCount,tiEventCount);
  

  currentFile++;
    }
    



  c1->Clear();
  sprintf(name,"%s_timestampdiff.png",inname.Data());
  printf("Drawing %s\n",name);
  timeStampDiffHist->Draw();
  c1->SaveAs(name);

  c1->Clear();
  sprintf(name,"%s_timestampdiffbin.png",inname.Data());
  printf("Drawing %s\n",name);
  timeStampDiffBinHist->Draw();
  c1->SaveAs(name);

  c1->Clear();
  sprintf(name,"%s_timestampmod6.png",inname.Data());
  printf("Drawing %s\n",name);
  timeStampHistMod6->Draw("E");
  c1->SaveAs(name);
  c1->Clear();
  sprintf(name,"%s_timestampmod12.png",inname.Data());
  printf("Drawing %s\n",name);
  timeStampHistMod12->Draw("E");
  c1->SaveAs(name);
  c1->Clear();
  sprintf(name,"%s_timestampmod24.png",inname.Data());
  printf("Drawing %s\n",name);
  timeStampHistMod24->Draw("E");
  c1->SaveAs(name);


  
  TH2F* pulseShapeHistPerHybrid;
  for (int rce = 0;rce<MAX_RCE;rce++) {
    for (int fpga = 0;fpga<MAX_FEB;fpga++) {
      for (int hyb = 0;hyb<MAX_HYB;hyb++) {
        if (hybridCount[rce][fpga][hyb]) {

          bool doTimeEvo = false;
          
          if(use_rce!=-1 && use_rce==rce ) {
            if(use_fpga!=-1 && use_fpga==fpga ) {
              if(use_hybrid!=-1 && use_hybrid==hyb ) {
                doTimeEvo = true;
              }
            }
          }
          
          
          pulseShapeHistPerHybrid = NULL;
          
          if(doPulseShape) {
            for (int ch=0;ch<640;ch++) {
            
              if(  pulseShapeHist[rce][fpga][hyb][ch] != NULL) {             
                if(pulseShapeHistPerHybrid==NULL) {
                  sprintf(name,"%s_pulseshape_R%d_F%d_H%d.png",inname.Data(),rce,fpga,hyb);         
                  pulseShapeHistPerHybrid = (TH2F*)  pulseShapeHist[rce][fpga][hyb][ch]->Clone(name);
                } 
                else {
                  pulseShapeHistPerHybrid->Add(pulseShapeHist[rce][fpga][hyb][ch]);
                }
              
                // do individual channels if wanted
                if( use_chan !=-1 && use_chan == ch) {
                  sprintf(name,"%s_pulseshape_R%d_F%d_H%d_CH%d.png",inname.Data(),rce,fpga,hyb,ch);
                  printf("Drawing %s\n",name);
                  c1->Clear();
                  printf("Drawing %s\n",name);
                  pulseShapeHist[rce][fpga][hyb][ch]->Draw("colz");
                  printf("Drawing %s\n",name);
                  c1->SaveAs(name);
                  printf("Saved %s\n",name);
                }             
              }
            
            } //ch
          
          

            if( pulseShapeHistPerHybrid != NULL ) {
              sprintf(name,"%s_pulseshape_R%d_F%d_H%d.png",inname.Data(),rce,fpga,hyb);
              printf("Drawing %s\n",name);
              c1->Clear();
              printf("Drawing %s\n",name);
              pulseShapeHistPerHybrid->Draw("colz");
              printf("Drawing %s\n",name);
              c1->SaveAs(name);
              printf("Saved %s\n",name);
              delete pulseShapeHistPerHybrid;
            }
          }


          c1->Clear();          
          sprintf(name,"%s_baseline_R%d_F%d_H%d.png",inname.Data(),rce,fpga,hyb);
          printf("Drawing %s\n",name);
          baselineHist2D[rce][fpga][hyb]->Draw("colz");
          c1->SaveAs(name);


          c1->Clear();          
          sprintf(name,"%s_maxsamplephasemod24_R%d_F%d_H%d.png",inname.Data(),rce,fpga,hyb);
          printf("Drawing %s\n",name);
          maxSamplePhaseMod24Hist2D[rce][fpga][hyb]->SetStats(false);
          maxSamplePhaseMod24Hist2D[rce][fpga][hyb]->Draw("colz,text");
          c1->SaveAs(name);
          

          TH1F* occupancy;


          c1->Clear();          
          sprintf(name,"%s_occupancy_R%d_F%d_H%d.png",inname.Data(),rce,fpga,hyb);
          printf("Drawing %s\n",name);
          occupancy = (TH1F*) baselineHist2D[rce][fpga][hyb]->ProjectionX(name);
          occupancy->SetDirectory(0);
          if (tiEventCount>0) occupancy->Scale(1.0/tiEventCount);
          else  occupancy->Scale(0.0); // set to 0          
          occupancy->Draw();
          //c1->SetLogy();
          c1->SaveAs(name);
          // for(int bin=1; bin<prjY->GetNbinsX()+1; ++bin) {
          //   double n = prjY->GetBinContent(bin);
          //   double occ = 0.0;
          //   if (eventCount>0) {
          //     occ = n/eventCount;
          //   }
          //   occupancy[sample]->ProjectionY();
          // }
          // delete prjY;
          
          delete occupancy;




          for(int sample=0; sample<6; ++sample) {
            c1->Clear();            
            sprintf(name,"%s_baselinesamples_R%d_F%d_H%d_S%d.png",inname.Data(),rce,fpga,hyb,sample);
            printf("Drawing %s\n",name);
            baselineSamplesHist2D[rce][fpga][hyb][sample]->Draw("colz");
            c1->SaveAs(name);
          }

          
          c2->Clear();
          TLegend* leg = new TLegend(0.35,0.5,0.73,0.85);
          leg->SetFillColor(0);
          leg->SetFillStyle(0);
          c3->Clear();
          c3->cd();
          sprintf(name,"baselinesamples0timebins_profile_R%d_F%d_H%d",rce,fpga,hyb);
          TH2F* prf_temp = new TH2F(name,name,5,0,640,10,3000.0,7000.0);
          prf_temp->SetDirectory(0);
          prf_temp->SetStats(false);
          prf_temp->Draw();
          //TLegend* leg3 = new TLegend(0.35,0.5,0.73,0.85);
          //leg3->SetFillColor(0);
          Color_t color;
          TProfile* prfS0[nTimeBins];          

          sprintf(name,"timebin_vs_occupancy_R%d_F%d_H%d",rce,fpga,hyb);
          TH2F* histTimeEvo2D = new TH2F(name, name,640,0,640,nTimeBins,0,nTimeBins);
          histTimeEvo2D->SetDirectory(0);
       
          double largeTiDtOcc[640];
          double highestOcc[640];

          for(int i=0;i<640;++i) {
            largeTiDtOcc[i]=0.;
            highestOcc[i] = 0.;
          }
          
          
          TH1F* occupancytimebin[nTimeBins];

          TH1F* baselineSample0TimeBinPerCh[640][nTimeBins];
          for(int i=0;i<640;++i) {
            for(int j=0;j<nTimeBins;++j) {
              baselineSample0TimeBinPerCh[i][j] = NULL;
            }
          }
          
          


          for(int timebin=0; timebin<nTimeBins; ++timebin) {
            
            c1->cd();
            c1->Clear();            
            sprintf(name,"%s_baselinesamples0timebins_R%d_F%d_H%d_T%d.png",inname.Data(),rce,fpga,hyb,timebin);
            printf("Drawing %s\n",name);
            baselineSamples0TimeBinHist2D[rce][fpga][hyb][timebin]->Draw("colz");
            c1->SaveAs(name);

            
            // do per channel plot for selected hybrid only
            if(doTimeEvo) {
              for(int xbin=1; xbin<=baselineSamples0TimeBinHist2D[rce][fpga][hyb][timebin]->GetNbinsX();++xbin) {
                int ch = xbin-1;
                sprintf(name,"%s_baselinesample0_timebins_R%d_F%d_H%d_CH%d_T%d.png",inname.Data(),rce,fpga,hyb,ch,timebin);
                baselineSample0TimeBinPerCh[ch][timebin] = (TH1F*) baselineSamples0TimeBinHist2D[rce][fpga][hyb][timebin]->ProjectionY(name,xbin,xbin,"E");   
                baselineSample0TimeBinPerCh[ch][timebin]->SetDirectory(0);
                
                /*
                c1->cd();
                c1->Clear();          
                printf("Drawing %s\n",name); 
                baselineSample0TimeBinPerCh[ch][timebin]->Draw("hist");
                //myText(0.6,0.6,TString::Format("%d events in timebin",tiEventCountTimeBin).Data(),0.05,1);
                //c1->SetLogy();
                c1->SaveAs(name);


                if(grbaselinesample0timebin[ch]==NULL) {
                  grbaselinesample0timebin[ch] = new TGraphErrors();
                  grbaselinermssample0timebin[ch] = new TGraphErrors();
                }
                if(baselineSample0TimeBinPerCh[ch][timebin]->GetEntries()>2) {
                  double t = timeBins[timebin];
                  int n =  grbaselinesample0timebin[ch]->GetN();
                  grbaselinesample0timebin[ch]->SetPoint(n,t,baselineSample0TimeBinPerCh[ch][timebin]->GetMean());
                  grbaselinesample0timebin[ch]->SetPointError(n,0,baselineSample0TimeBinPerCh[ch][timebin]->GetMeanError());
                  grbaselinermssample0timebin[ch]->SetPoint(n,t,baselineSample0TimeBinPerCh[ch][timebin]->GetRMS());
                  grbaselinermssample0timebin[ch]->SetPointError(n,0,baselineSample0TimeBinPerCh[ch][timebin]->GetRMSError());
                }
                */

              }
            }
            
            

            sprintf(name,"%s_occupancy_timebins_R%d_F%d_H%d_T%d.png",inname.Data(),rce,fpga,hyb,timebin);
            occupancytimebin[timebin] = (TH1F*) baselineSamples0TimeBinHist2D[rce][fpga][hyb][timebin]->ProjectionX(name);   
            //occupancytimebin[timebin] = (TH1F*)baselinesample0timebin[timebin]->Clone(name);
            occupancytimebin[timebin]->SetDirectory(0);
            // find the number of ti events in this time bin to normalize correctly
            uint tiEventCountTimeBin = timeStampDiffBinHist->GetBinContent(timebin+1);
            if (tiEventCountTimeBin>0) occupancytimebin[timebin]->Scale(1.0/tiEventCountTimeBin);
            else  occupancytimebin[timebin]->Scale(0.0); // set to 0          
            
            c1->cd();
            c1->Clear();          
            occupancytimebin[timebin]->Draw("hist");
            myText(0.6,0.6,TString::Format("%d events in timebin",tiEventCountTimeBin).Data(),0.05,1);
            //c1->SetLogy();
            c1->SaveAs(name);
            

            for(int xbin=1; xbin<=occupancytimebin[timebin]->GetNbinsX();++xbin) {
              double occ = occupancytimebin[timebin]->GetBinContent(xbin);
              double err = occupancytimebin[timebin]->GetBinError(xbin);
              histTimeEvo2D->SetBinContent(xbin,timebin+1,occ);
              histTimeEvo2D->SetBinError(xbin,timebin+1,err);
              
              // find the largest occupancy among large trig time diff
              if(timeBins[timebin]>30.) {                
                if(occ > largeTiDtOcc[xbin-1]) {
                  largeTiDtOcc[xbin-1] = occ;
                }
              }
              // find the largest fluctutation in occupancy
              if(occ > highestOcc[xbin-1]) {
                highestOcc[xbin-1] = occ;
              }
              
            }
            
            prfS0[timebin] = baselineSamples0TimeBinHist2D[rce][fpga][hyb][timebin]->ProfileX();

            if(timebin<9) color = (Color_t)timebin+1;
            else color = (Color_t)timebin+2;
            occupancytimebin[timebin]->SetLineColor(color);
            prfS0[timebin]->SetLineColor(color);
            prfS0[timebin]->SetMarkerColor(color);


            if(timebin==0) {
              c2->cd();
              occupancytimebin[timebin]->Draw("hist");
              c3->cd();
              prfS0[timebin]->Draw("same");
              
            }
            else {
              c2->cd();
              occupancytimebin[timebin]->Draw("hist,same");
              c3->cd();
              prfS0[timebin]->Draw("same");
            }
            leg->AddEntry( occupancytimebin[timebin], TString::Format("#DeltaT_{bin} %d",timebin).Data(),"L");
            

          } // timebin

          
          if(doTimeEvo) {

            
            
            TH1F* baselineSample0TimeBinPerChAll[nTimeBins];
            TGraphErrors* grbaselinesample0timebin[640];
            TGraphErrors* grbaselinermssample0timebin[640];
            for(int i=0;i<640;++i) {
              grbaselinesample0timebin[i] = NULL;
              grbaselinermssample0timebin[i] = NULL;
            }
            
            for(int timebin=0;timebin<nTimeBins;++timebin) baselineSample0TimeBinPerChAll[timebin] = NULL;
            
            for(int ch=0;ch<640;++ch) {
              if(largeTiDtOcc>0) {
                if(highestOcc[ch]/largeTiDtOcc[ch]>thresholdOcc) {
                  c1->cd();
                  c1->Clear();


                  
                  
                  if(grbaselinesample0timebin[ch]==NULL) {
                    grbaselinesample0timebin[ch] = new TGraphErrors();
                    grbaselinermssample0timebin[ch] = new TGraphErrors();
                  }

                  for(int timebin=0;timebin<nTimeBins;++timebin) {
                    
                    if(baselineSample0TimeBinPerCh[ch][timebin]->GetEntries()>2) {
                      double t = timeBins[timebin];
                      int n =  grbaselinesample0timebin[ch]->GetN();
                      grbaselinesample0timebin[ch]->SetPoint(n,t,baselineSample0TimeBinPerCh[ch][timebin]->GetMean());
                      grbaselinesample0timebin[ch]->SetPointError(n,0,baselineSample0TimeBinPerCh[ch][timebin]->GetMeanError());
                      grbaselinermssample0timebin[ch]->SetPoint(n,t,baselineSample0TimeBinPerCh[ch][timebin]->GetRMS());
                      grbaselinermssample0timebin[ch]->SetPointError(n,0,baselineSample0TimeBinPerCh[ch][timebin]->GetRMSError());
                    }

                    // merge them
                    if( baselineSample0TimeBinPerChAll[timebin] == NULL) {
                      baselineSample0TimeBinPerChAll[timebin] = (TH1F*) baselineSample0TimeBinPerCh[ch][timebin]->Clone(TString::Format("%sAll",baselineSample0TimeBinPerCh[ch][timebin]->GetName()).Data());
                      baselineSample0TimeBinPerChAll[timebin]->SetDirectory(0);
                      baselineSample0TimeBinPerChAll[timebin]->Reset();
                    }
                    baselineSample0TimeBinPerChAll[timebin]->Add(baselineSample0TimeBinPerCh[ch][timebin]);
                    
                  }
                  
                  sprintf(name,"%s_baselinesample0mean_timebins_R%d_F%d_H%d_CH%d.png",inname.Data(),rce,fpga,hyb,ch);
                  grbaselinesample0timebin[ch]->SetMarkerSize(1.0);
                  grbaselinesample0timebin[ch]->SetMarkerStyle(20);
                  grbaselinesample0timebin[ch]->Draw("ALP");                  
                  c1->SaveAs(name);
                  
                  c1->cd();
                  c1->Clear();
                  sprintf(name,"%s_baselinesample0rms_timebins_R%d_F%d_H%d_CH%d.png",inname.Data(),rce,fpga,hyb,ch);
                  grbaselinermssample0timebin[ch]->SetMarkerSize(1.0);
                  grbaselinermssample0timebin[ch]->SetMarkerStyle(20);
                  grbaselinermssample0timebin[ch]->Draw("ALP");
                  c1->SaveAs(name);
                  

                }
              } //largeOccDiff
            } //ch
            
            for(int i=0;i<640;++i) {
              if(grbaselinesample0timebin[i]!=NULL) {
                delete grbaselinesample0timebin[i];
                delete grbaselinermssample0timebin[i];
              }
            }

            TGraphErrors* grbaselinesample0timebinAll = new TGraphErrors();
            TGraphErrors* grbaselinermssample0timebinAll = new TGraphErrors();

            
            for(int timebin=0;timebin<nTimeBins;++timebin) {
              if( baselineSample0TimeBinPerChAll[timebin] != NULL) {
                c1->Clear();    
                sprintf(name,"%s_baselinesamples0_timebins_R%d_F%d_H%d_CHALL_T%d.png",inname.Data(),rce,fpga,hyb,timebin);
                printf("Drawing %s\n",name);
                c1->cd();
                baselineSample0TimeBinPerChAll[timebin]->Draw();
                c1->SaveAs(name);
                
                if(baselineSample0TimeBinPerChAll[timebin]->GetEntries()>2) {
                  
                  double t = timeBins[timebin];
                  int n =  grbaselinesample0timebinAll->GetN();
                  grbaselinesample0timebinAll->SetPoint(n,t,baselineSample0TimeBinPerChAll[timebin]->GetMean());
                  grbaselinesample0timebinAll->SetPointError(n,0,baselineSample0TimeBinPerChAll[timebin]->GetMeanError());
                  grbaselinermssample0timebinAll->SetPoint(n,t,baselineSample0TimeBinPerChAll[timebin]->GetRMS());
                  grbaselinermssample0timebinAll->SetPointError(n,0,baselineSample0TimeBinPerChAll[timebin]->GetRMSError());
                }
                delete baselineSample0TimeBinPerChAll[timebin];
              }
            }

            c1->Clear();    
            sprintf(name,"%s_baselinesamples0mean_timebins_R%d_F%d_H%d_CHALL.png",inname.Data(),rce,fpga,hyb);
            printf("Drawing %s\n",name);
            c1->cd();
            grbaselinesample0timebinAll->SetMarkerSize(1.0);
            grbaselinesample0timebinAll->SetMarkerStyle(20);
            grbaselinesample0timebinAll->Draw("ALP");
            c1->SaveAs(name);
            delete grbaselinesample0timebinAll;

            c1->Clear();    
            sprintf(name,"%s_baselinesamples0rms_timebins_R%d_F%d_H%d_CHALL.png",inname.Data(),rce,fpga,hyb);
            printf("Drawing %s\n",name);
            c1->cd();
            grbaselinermssample0timebinAll->SetMarkerSize(1.0);
            grbaselinermssample0timebinAll->SetMarkerStyle(20);
            grbaselinermssample0timebinAll->Draw("ALP");
            c1->SaveAs(name);
            delete grbaselinermssample0timebinAll;
                
            
            
          } //doTimeEvo
          

          for(int i=0;i<640;++i) {
            for(int j=0;j<nTimeBins;++j) {
              if(baselineSample0TimeBinPerCh[i][j]!=NULL) {
                delete baselineSample0TimeBinPerCh[i][j];
              }
            }
          }


          sprintf(name,"%s_occupancy_timebinsAll_R%d_F%d_H%d.png",inname.Data(),rce,fpga,hyb);
          printf("Drawing %s\n",name);
          c2->cd();
          // set y-scale
          double m=-1;
          for(int i=0;i<nTimeBins;++i) {
            if(occupancytimebin[i]->GetEntries()>0 && occupancytimebin[i]->GetMaximum()>m) {
              m = occupancytimebin[i]->GetMaximum();
            }
          }
          occupancytimebin[0]->SetMaximum(m);
          leg->Draw();
          c2->SaveAs(name);
          
          
          
          
          sprintf(name,"%s_baselinesamples0_profile_timebinsAll_R%d_F%d_H%d.png",inname.Data(),rce,fpga,hyb);
          printf("Drawing %s\n",name);
          c3->cd();
          leg->Draw();
          c3->SaveAs(name);
          

          c1->Clear();       
          sprintf(name,"%s_timebin_vs_occupancy_R%d_F%d_H%d.png",inname.Data(),rce,fpga,hyb);
          printf("Drawing %s\n",name);
          histTimeEvo2D->Draw("colz");
          c1->SaveAs(name);


          for(int i=0;i<640;++i) {
            printf("ch %d: largeDt occ %f \n",i, largeTiDtOcc[i]);
          }
          
          
          for(int bin=1; bin<=histTimeEvo2D->GetNbinsX();++bin) {
            int ch = bin-1;
            
            if(doTimeEvo==false) {
              // plot time evolution for channels with large occupancy compared to at large trigger delta's              
              if(largeTiDtOcc[ch] > 0.) {
                if( highestOcc[ch]/largeTiDtOcc[ch] >thresholdOcc) {
                  doTimeEvo = true;
                  printf("largeDt occ %f highest_occ %f\n",largeTiDtOcc[ch], highestOcc[ch]);                  
                } 
              }
            }
            
            //if(ch%127==0 || ch%50==0 || diTimeEvo) {
            if( doTimeEvo) {
              
              if(largeTiDtOcc[ch] > 0.) {
                if( highestOcc[ch]/largeTiDtOcc[ch] >thresholdOcc) {
                  printf("largeDt occ %f highest_occ %f\n",largeTiDtOcc[ch], highestOcc[ch]);                  

                  TGraphErrors* grTimeEvo = new TGraphErrors();
                  for(int ybin=1; ybin<=histTimeEvo2D->GetNbinsY();++ybin) {
                    double t = timeBins[ybin-1];
                    grTimeEvo->SetPoint(ybin-1,t,histTimeEvo2D->GetBinContent(bin,ybin));
                    grTimeEvo->SetPointError(ybin-1,0., histTimeEvo2D->GetBinError(bin,ybin));
                  }
                  c1->Clear();    
                  sprintf(name,"occupancy_vs_timebin_ch%d_R%d_F%d_H%d.png",ch,rce,fpga,hyb);
                  printf("Drawing %s\n",name);
                  grTimeEvo->SetTitle(name);
                  grTimeEvo->SetMarkerSize(1.0);
                  grTimeEvo->SetMarkerStyle(20);
                  grTimeEvo->Draw("ALP");
                  sprintf(name,"%s_occupancy_vs_timebin_ch%d_R%d_F%d_H%d.png",inname.Data(),ch,rce,fpga,hyb);
                  c1->SaveAs(name);
                  delete grTimeEvo;
                  
                } 
              }

            }

            
          }



          delete histTimeEvo2D;
          

          // for(int sample=0; sample<6; ++sample) {
          //   c1->Clear();            
          //   sprintf(name,"%s_baselinesamples_R%d_F%d_H%d_S%d.png",inname.Data(),rce,fpga,hyb,sample);
          //   printf("Drawing %s\n",name);
          //   c1->SaveAs(name);
          // }


          delete leg;
          
          for(int timebin=0;timebin<nTimeBins;++timebin) {
            delete prfS0[timebin];
            delete occupancytimebin[timebin];
          }


          
        }
      }
    }
  }
  // Start X-Windows
  //theApp.Run();
  
  // Close file
  delete dataRead;
  return(0);
}

int compareInts(const void* a, const void* b) {
  if ( *(uint*)a <  *(uint*)b ) return -1;
  else if ( *(uint*)a == *(uint*)b ) return 0;
  else if ( *(uint*)a >  *(uint*)b ) return 1;
  else {
    printf("Should neven get here!!\n");
    exit(1);
  }
}

uint getTimeBin(double t) {
  // want to use bsearch but need to return index not just if it's in the array?

  if(t<0.0) {
    printf("I should never be here\n");
    exit(1);
  }
  
  uint bin = -1;
  for(uint i =0; i<nTimeBins-1; ++i) {
    if( t < timeBins[i+1] && t >= timeBins[i] ) {
      bin = i;
      break;
    }     
  }
  if(bin==-1) {
    //printf("this time diff is large? %f \n", t);
    bin = nTimeBins-1;
    //exit(1);
  }
  return bin;

  // uint *pTime;
  // uint* key = &t;
  // pTime = (uint*) bsearch(key, timeBins, nTimeBins, sizeof(uint), compareTimes);
  // if(pTime!=NULL) {
  //   printf("pTime for t=%d is %d!\n",t,*pTime);
    
  // } else {
  //   printf("pTime for t=%d is null!\n",t);
  //   exit(1);
  // }
  // return *pTime;
}

