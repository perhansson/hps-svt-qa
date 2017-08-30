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
#include "cosmic_utils.hh"
using namespace std;

// Process the data
// Pass root file to open as first and only arg.
int main ( int argc, char **argv ) {
  TCanvas         *c1, *c2, *c3, *c4, *c5, *c6, *c7, *c8, *c9,*c10, *c11;
   TH2F            *histAll;
   TH2F            *histCos;
   TH2F            *histCosHgh;
   TH2F            *histCosPeak;
   TH2F            *histCosPos;
   TH1F            *histSng[640];
   double          histMin[640];
   double          histMax[640];
   double          histCosMin;
   double          histCosMax;
   double          histCosPosMin = 8192*3;
   double          histCosPosMax = -8192*3;
   double          allMin;
   double          allMax;
   TGraph          *mean;
   TGraph          *sigma;
   double          grChan[640];
   double          grMean[640];
   double          grSigma[640];
   DataRead        dataRead;
   DevboardEvent    event;
   DevboardSample   *sample;
   uint            x;
   uint            y;
   double          value;
   double          svalue[6];
   uint            channel;
   uint            tarChan;
   uint            tarFpga;
   uint            tarHybrid;
   uint            eventCount;
   char            name[100];
   bool            valid[640];
   uint            grCount;
   uint            scnt;
   uint            scnt_consec;
   uint            ecnt;
   bool            posSig;
   uint            verbose_level = 0;
   double          maxSample = -99.;
   double          histCosHghMin = 8192*3;
   double          histCosHghMax = -8192*3;
   double          histCosPeakMin = 8192*3;
   double          histCosPeakMax = -8192*3;
   double          min_sng = 99999.;
   bool            save = false;
   uint            Nevents = 0;
   uint            samplesAboveThres = 3;
   uint            consecutiveSamplesAboveThres = 3;
   uint            thresh = 3;
   double          peakSample = -1;
   bool            hasZeroSample = false;

   gStyle->SetOptStat(kFALSE);

   // Start X11 view
   TApplication theApp("App",NULL,NULL);

   // Root file is the first and only arg
   if ( argc < 6 ) {
      cout << "Usage: cosmic_summary fpga hybrid channel saveFile nevents base_file data_files \n";
      return(1);
   }
   tarFpga   = atoi(argv[1]);
   tarHybrid = atoi(argv[2]);
   tarChan   = atoi(argv[3]);
   save = atoi(argv[4]) != 0 ? true : false;
   Nevents = atoi(argv[5]);

   // 2d histogram
   histAll = new TH2F("Value_Hist_All","Value_Hist_All",16384,0,16384,640,0,640);
   histCos = new TH2F("Value_Hist_Cos","Value_Hist_Cos",16384,-8192,8192,640,0,640);
   histCosHgh = new TH2F("Value_Hist_Cos_Hgh","Value_Hist_Cos_Hgh",16384,-8192,8192,640,0,640);
   histCosPeak = new TH2F("Value_Hist_Cos_Peak","Value_Hist_Cos_Peak",16384,-8192,8192,640,0,640);
   histCosPos = new TH2F("Value_Hist_Cos_Pos","Value_Hist_Cos_Pos",16384,-8192,8192,640,0,640);

   for (channel=0; channel < 640; channel++) {
      sprintf(name,"%i",channel);
      histSng[channel] = new TH1F(name,name,16384,0,16384);
      histMin[channel] = 16384;
      histMax[channel] = -16384;
      valid[channel]   = false;
   }
   allMin = 16384;
   allMax = 0;
   histCosMin = 8192*3;
   histCosMax = -8192*3;

   // Attempt to open data file
   if ( ! dataRead.open(argv[6]) ) return(2);

   // Process each event
   eventCount = 0;

   while ( dataRead.next(&event) ) {
      if ( eventCount == 0 ) {
         if ( dataRead.getConfig("cntrlFpga:hybrid:apv25:PreampPolarity") == "Inverting") {
            posSig = true;
            cout << "Detected positive signals" << endl;
         }
         else {
            posSig = false;
            cout << "Detected negative signals" << endl;
         }
      }

      //if ( eventCount >= Nevents && Nevents > 0 ) break;
      if ( eventCount >= 10000 ) break;

      if ( eventCount % 1000 ==0 ) cout << "Processing event " << eventCount << endl;

      for (x=0; x < event.count(); x++) {

         // Check for matching FPGA
         if ( event.fpgaAddress() == tarFpga ) {

            // Get sample
            sample  = event.sample(x);

            // Check for matching hybrid
            if ( sample->hybrid() == tarHybrid ) {
               channel = (sample->apv() * 128) + (127-sample->channel());

               // Check for out of range channels
               if ( channel >= (5 * 128) ) {
                  cout << "Channel " << dec << channel << " out of range" << endl;
                  cout << "Apv = " << dec << sample->apv() << endl;
                  cout << "Chan = " << dec << sample->channel() << endl;
               } else {
                  valid[channel] = true;

                  // Filter APVs
                  if ( eventCount > 100 ) {
                     for ( y=0; y < 6; y++ ) {
                        value = sample->value(y);
                        histAll->Fill(value,channel);
                        histSng[channel]->Fill(value);
			if ( value < 1 ) {
			  cout << "Warning: Sample " << y << " channel " << channel << " has " << value << endl; 
			}
			if ( value == 8192 ) {
			  cout << "Warning: Sample " << y << " channel " << channel << " has " << value << endl; 
			}
                        if ( value < histMin[channel] ) histMin[channel] = value;
                        if ( value > histMax[channel] ) histMax[channel] = value;
                        if ( value < allMin           ) allMin = value;
                        if ( value > allMax           ) allMax = value;
                     }
                  }
               }
            }
         }
      }
      eventCount++;
   }

   // Fit histograms
   ofstream ofs("thresholds.dat",ios::out);
   if (ofs.is_open() ) {
     ofs << tarFpga << "," << tarHybrid << "," << 0 << endl;
   }
   grCount = 0;
   for(channel = 0; channel < 640; channel++) {
      if ( valid[channel] ) {
         histSng[channel]->Fit("gaus");
         histSng[channel]->GetXaxis()->SetRangeUser(histMin[channel],histMax[channel]);
         histSng[channel]->Draw();
         grMean[channel]  = histSng[channel]->GetFunction("gaus")->GetParameter(1);
         grSigma[channel] = histSng[channel]->GetFunction("gaus")->GetParameter(2);
         grChan[channel]  = channel;
      } else {
         grMean[channel]  = 0;
         grSigma[channel] = 0;
         grChan[channel]  = channel;
      }
      
      ofs << channel << "," << grMean[channel] << endl;

   }

   ofs.close();


   c1 = new TCanvas("c1","c1");
   c1->cd();
   histAll->Draw("colz");
   histAll->GetXaxis()->SetRangeUser(allMin,allMax);
   
   c2 = new TCanvas("c2","c2");
   c2->cd();
   histSng[tarChan]->Draw();

   c3 = new TCanvas("c3","c3");
   c3->cd();
   mean = new TGraph(640,grChan,grMean);
   mean->Draw("a*");

   c4 = new TCanvas("c4","c4");
   c4->cd();
   sigma = new TGraph(640,grChan,grSigma);
   sigma->Draw("a*");

   // Attempt to open data file
   dataRead.close();

   // loop over files
   int nfread = 0;
   while( (nfread+7) < argc) {
     const char* fileName = argv[nfread+7];
     cout << fileName << endl;
     


   if ( ! dataRead.open(fileName) ) return(2);

   // Process each event
   eventCount = 0;
   ecnt       = 0;

   while ( dataRead.next(&event) ) {
      if ( eventCount == 0 ) {
         if ( dataRead.getConfig("cntrlFpga:hybrid:apv25:PreampPolarity") == "Inverting") {
            posSig = true;
            cout << "Detected positive signals" << endl;
         }
         else {
            posSig = false;
            cout << "Detected negative signals" << endl;
         }
      }

      if ( eventCount >= Nevents && Nevents > 0 ) break;


      if ( eventCount % 1000 == 0 ) cout << "Processing event " << eventCount << endl;

      for (x=0; x < event.count(); x++) {

         // Check for matching FPGA
         if ( event.fpgaAddress() == tarFpga ) {

            // Get sample
            sample  = event.sample(x);

            // Check for matching hybrid
            if ( sample->hybrid() == tarHybrid ) {
               channel = (sample->apv() * 128) + (127-sample->channel());

               // Check for out of range channels
               if ( channel >= (5 * 128) ) {
                  cout << "Channel " << dec << channel << " out of range" << endl;
                  cout << "Apv = " << dec << sample->apv() << endl;
                  cout << "Chan = " << dec << sample->channel() << endl;
               } else {
                  valid[channel] = true;

                  // Filter APVs
                  scnt = 0;
		  maxSample = -99.;  
		  hasZeroSample = false;
                  for ( y=0; y < 6; y++ ) {
		    
		    //if ( sample->value(y) == 0 ) hasZeroSample = true;
		    
		    svalue[y] = (double)sample->value(y) - grMean[channel];
		    if ( sample->apv() != 0 && svalue[y] > (thresh*grSigma[channel]) ) {
		      scnt++;		     
		    }
		    if ( svalue[y] > maxSample ) maxSample = svalue[y];
		    
		    if ( svalue[y] < -5000. ) cout << "Sample " << y << " channel " << channel << " has " << svalue[y] << " subtr signal from " << sample->value(y) << " mean " << grMean[channel] << endl; 
		    
		    
                  }
		  
		  if ( !hasZeroSample ) { // remove events which has a zero sample...?
		    
		    if ( scnt >= samplesAboveThres ) {
		      
		      scnt_consec = cosmic_utils::getNrOfConsecutiveSamplesAboveThresh(svalue,thresh*grSigma[channel]);
		      
		      if ( scnt_consec >= consecutiveSamplesAboveThres ) {
			
			bool shapeOk = true;
			if(shapeOk) {
			  if ( verbose_level > 0 ) {
			    cout << "Found hit."
				 << " V0=" << svalue[0]
				 << ", V1=" << svalue[1]
				 << ", V2=" << svalue[2]
				 << ", V3=" << svalue[3]
				 << ", V4=" << svalue[4]
				 << ", V5=" << svalue[5]
				 << ", Event=" << dec << eventCount
				 << ", Count=" << dec << ecnt
				 << ", Channel=" << dec << channel
				 << ", Apv=" << dec << sample->apv() << endl;
			  }
		      
			  ecnt++;
		      
			  histCosHgh->Fill(maxSample,channel);
			  if ( maxSample < histCosHghMin ) histCosHghMin = maxSample;
			  if ( maxSample > histCosHghMax ) histCosHghMax = maxSample;
		      
			  min_sng = 99999.;
			  for ( y=0; y < 6; y++ ) {
			    histCos->Fill(svalue[y],channel);
			    if ( svalue[y] < histCosMin ) histCosMin = svalue[y];
			    if ( svalue[y] > histCosMax ) histCosMax = svalue[y];
			    if( svalue[y] < min_sng) min_sng = svalue[y];
			  }
		      
			  if ( min_sng > 0 ) {
			    for ( y=0; y < 6; y++ ) {
			      histCosPos->Fill(sample->value(y),channel);
			      if ( sample->value(y) < histCosPosMin ) histCosPosMin = sample->value(y);
			      if ( sample->value(y) > histCosPosMax ) histCosPosMax = sample->value(y);
			    }
			  }

			  // look for peak in 3 samples above threshold
			  if ( scnt_consec >= 3 ) {
			    peakSample = cosmic_utils::getMaxSampleOfThree(svalue);
			    if(peakSample>0) {
			      histCosPeak->Fill(peakSample,channel);
			      if ( maxSample < histCosPeakMin ) histCosPeakMin = maxSample;
			      if ( maxSample > histCosPeakMax ) histCosPeakMax = maxSample;			    
			    } 
			  }
			} //shapeOk
		      } //samplesabove
		    } // samplesabove
		  } //zero sample
               }
            }
         }
      }
      eventCount++;
   } //event loop
   ++nfread;
   } //files

   c5 = new TCanvas("c5","c5");
   c5->cd();
   histCos->GetXaxis()->SetRangeUser(histCosMin,histCosMax);
   histCos->Draw("colz");

   c6 = new TCanvas("c6","c6");
   c6->cd();
   histCosHgh->GetXaxis()->SetRangeUser(histCosHghMin,histCosHghMax);
   histCosHgh->Draw("colz");

   c7 = new TCanvas("c7","c7");
   c7->cd();
   histCosPos->GetXaxis()->SetRangeUser(histCosPosMin,histCosPosMax);
   histCosPos->Draw("colz");

   char buf[50];

   c8 = new TCanvas("c8","c8");
   c8->cd();
   sprintf(buf,"%s_prjy",histCos->GetName());
   TH1F* histCosPrj = (TH1F*)histCos->ProjectionX(buf,1+128*1,1+128*3);
   //histCosPrj->GetXaxis()->SetRangeUser(histCosMin,histCosMax);
   histCosPrj->GetXaxis()->SetRangeUser(-100,1000);
   histCosPrj->Draw("");

   c9 = new TCanvas("c9","c9");
   c9->cd();
   sprintf(buf,"%s_prjy",histCosHgh->GetName());
   TH1F* histCosHghPrj = (TH1F*)histCosHgh->ProjectionX(buf,1+128*1,1+128*3);
   //histCosHghPrj->GetXaxis()->SetRangeUser(histCosHghMin,histCosHghMax);
   histCosHghPrj->GetXaxis()->SetRangeUser(-100,1000);
   histCosHghPrj->Draw("");

   c10 = new TCanvas("c10","c10");
   c10->cd();
   histCosPeak->GetXaxis()->SetRangeUser(histCosPeakMin,histCosPeakMax);
   histCosPeak->Draw("colz");

   c11 = new TCanvas("c11","c11");
   c11->cd();
   sprintf(buf,"%s_prjy",histCosPeak->GetName());
   TH1F* histCosPeakPrj = (TH1F*)histCosPeak->ProjectionX(buf,1+128*1,1+128*3);
   //histCosPeakPrj->GetXaxis()->SetRangeUser(histCosPeakMin,histCosPeakMax);
   histCosPeakPrj->GetXaxis()->SetRangeUser(-100,1000);
   histCosPeakPrj->Draw("");



   if ( save ) {
     sprintf(buf,"%s.pdf",c1->GetName());
     c1->SaveAs(buf);
     sprintf(buf,"%s.pdf",c2->GetName());
     c2->SaveAs(buf);
     sprintf(buf,"%s.pdf",c3->GetName());
     c3->SaveAs(buf);
     sprintf(buf,"%s.pdf",c4->GetName());
     c4->SaveAs(buf);
     sprintf(buf,"%s.pdf",c5->GetName());
     c5->SaveAs(buf);
     sprintf(buf,"%s.pdf",c6->GetName());
     c6->SaveAs(buf);
     sprintf(buf,"%s.pdf",c7->GetName());
     c7->SaveAs(buf);
     sprintf(buf,"%s.pdf",c8->GetName());
     c8->SaveAs(buf);
     sprintf(buf,"%s.pdf",c9->GetName());
     c9->SaveAs(buf);
     sprintf(buf,"%s.pdf",c10->GetName());
     c10->SaveAs(buf);
     sprintf(buf,"%s.pdf",c11->GetName());
     c11->SaveAs(buf);
   }
   // Start X-Windows
   theApp.Run();

   // Close file
   dataRead.close();
   return(0);
}

