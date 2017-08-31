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

// Process the data
// Pass root file to open as first and only arg.
int main ( int argc, char **argv ) {
   TCanvas         *c1, *c2;
   double          grY[128*6];
   double          grX[128*6];
   TGraph          *plot;
   DataRead        dataRead;
   DevboardEvent    event;
   DevboardSample   *sample;
   uint            idx;
   uint            x;
   uint            y;
   uint            value;
   uint            channel;
   uint            tarFpga;
   uint            tarHybrid;
   uint            tarApv;
   TH1F            *hist[60];
   double          histMin[60];
   double          histMax[60];
   stringstream    title;

   gStyle->SetOptStat(kFALSE);

   // Start X11 view
   TApplication theApp("App",NULL,NULL);

   // Root file is the first and only arg
   if ( argc != 5 ) {
      cout << "Usage: plot_readout fpga hybrid apv data_file\n";
      return(1);
   }
   tarFpga   = atoi(argv[1]);
   tarHybrid = atoi(argv[2]);
   tarApv    = atoi(argv[3]);

   for (x=0; x<60; x++) {
      title.str("");
      title << x;
      hist[x] = new TH1F(title.str().c_str(),title.str().c_str(),16384,0,16384);
      histMin[x] = 16384;
      histMax[x] = 0;
    }

   // Attempt to open data file
   if ( ! dataRead.open(argv[4]) ) return(2);

   c1 = new TCanvas("1","1");

   while ( dataRead.next(&event) ) {

      // Check for matching FPGA
      if ( event.fpgaAddress() == tarFpga ) {

         for (x=0; x<128*6; x++) {
            grX[x] = x;
            grY[x] = 0;
         }

         for (x=0; x < 128; x++) {

            // Get sample
            sample  = event.sample(x);

            // Check for matching hybrid
            if ( sample->hybrid() == tarHybrid && sample->apv() == tarApv ) {
               channel = sample->channel();

               //for (y=0;y<6;y++) {
                  //idx   = y * 20 + x;
                  idx   = x;
                  //value = sample->value(y);
                  value = sample->value(0);

                  grX[idx] = idx;
                  grY[idx] = value;
//
                  //hist[idx]->Fill(value);
                  ////if ( value < histMin[idx] ) histMin[idx] = value;
                  //if ( value > histMax[idx] ) histMax[idx] = value;
               //}
            }
         }

         c1->cd();
         plot = new TGraph(128,grX,grY);
         plot->GetYaxis()->SetRangeUser(5000,9000);
         plot->Draw("a*");
         c1->Update();
         cout << "Plot" << endl;
         //break;
         sleep(1);

      }
   }

   c2 = new TCanvas("2","2");
   c2->Divide(4,6,.01,.01);

   for (x=0; x< 24; x++) {
      c2->cd(x+1);
      hist[x]->Draw();
      hist[x]->GetXaxis()->SetRangeUser(histMin[x],histMax[x]);
   }

   // Start X-Windows
   theApp.Run();

   // Close file
   dataRead.close();
   return(0);
}

