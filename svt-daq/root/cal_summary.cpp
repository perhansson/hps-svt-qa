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
#include <iomanip>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TStyle.h>
#include <stdarg.h>
#include <DevboardEvent.h>
#include <DevboardSample.h>
#include <Data.h>
#include <DataRead.h>
using namespace std;

int chanMap[128];

void initChan ( ) {
   int idx;
   int chan;

   for ( idx = 0; idx < 128; idx++ ) {
      chan = (32*(idx%4)) + (8*(idx/4)) - (31*(idx/16));
      chanMap[chan] = idx;
   }
}


int convChan ( int chan ) {
   //return(chanMap[chan]);
   return(chan);
}


// Process the data
// Pass root file to open as first and only arg.
int main ( int argc, char **argv ) {
   TCanvas         *c1, *c2, *c3, *c4, *c5;
   TH2F            *histAll;
   TH1F            *histSng[128];
   double          histMin[128];
   double          histMax[128];
   TGraph          *gr;
   TGraph          *mean;
   TGraph          *sigma;
   double          grChan[128];
   double          grMean[128];
   double          grRange[128];
   DataRead        dataRead;
   DevboardEvent    event;
   DevboardSample   *sample;
   double          grY[6000];
   double          grX[6000];
   uint            grCnt; 
   uint            x;
   uint            y;
   uint            value;
   uint            channel;
   uint            tar;
   uint            eventCount;
   double          avg;
   char            name[100];

   initChan();

   gStyle->SetOptStat(kFALSE);

   // Start X11 view
   TApplication theApp("App",NULL,NULL);

   // Root file is the first and only arg
   if ( argc != 3 ) {
      cout << "Usage: hist_summary channel data_file\n";
      return(1);
   }
   tar = atoi(argv[1]);

   // 2d histogram
   histAll = new TH2F("Value_Hist_All","Value_Hist_All",16384,0,16384,128,0,128);

   for (channel=0; channel < 128; channel++) {
      sprintf(name,"%i",channel);
      histSng[channel] = new TH1F(name,name,16384,0,16384);
      histMin[channel] = 16384;
      histMax[channel] = 0;
   }

   // Attempt to open data file
   if ( ! dataRead.open(argv[2]) ) return(2);

   // Process each event
   eventCount = 0;
   grCnt = 0;
   while ( dataRead.next(&event) ) {

      for (x=0; x < event.count(); x++) {

         // Get sample
         sample  = event.sample(x);
         channel = sample->channel();

         // Filter APVs
         if ( eventCount > 20 ) {

            avg = 0;
            for ( y=0; y < 6; y++ ) {
               value = sample->value(y);
               
               histAll->Fill(value,channel);
               histSng[channel]->Fill(value);

               if ( value < histMin[channel] ) histMin[channel] = value;
               if ( value > histMax[channel] ) histMax[channel] = value;

               if ( channel == tar && eventCount < 1000 ) {
                  grY[grCnt] = value;
                  grX[grCnt] = grCnt;
                  grCnt++;
               }
            }
         }
      }
      eventCount++;

   }

   for(channel = 0; channel < 128; channel++) {
      grMean[channel]  = (histMax[channel] + histMin[channel])/2.0;
      grRange[channel] = (histMax[channel] - histMin[channel])/2.0;
      grChan[channel]  = channel;
   }

   c1 = new TCanvas("c11","c11");
   c1->cd();
   histAll->Draw("colz");

   c2 = new TCanvas("c12","c12");
   c2->cd();
   histSng[tar]->GetXaxis()->SetRangeUser(histMin[tar],histMax[tar]);
   histSng[tar]->Draw();

   c3 = new TCanvas("c13","c13");
   c3->cd();
   gr = new TGraph(grCnt,grX,grY);
   gr->Draw("a*");

   c4 = new TCanvas("c14","c14");
   c4->cd();
   mean = new TGraph(128,grChan,grMean);
   mean->Draw("a*");

   c5 = new TCanvas("c15","c15");
   c5->cd();
   sigma = new TGraph(128,grChan,grRange);
   sigma->Draw("a*");

   // Start X-Windows
   theApp.Run();

   // Close file
   dataRead.close();
   return(0);
}

