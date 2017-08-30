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
#include <TFile.h>
#include <TH1F.h>
#include <meeg_utils.hh>
#include <TH2I.h>
#include <TF1.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
//#include <TApplication.h>
#include <TGraph.h>
#include <TStyle.h>
#include <stdarg.h>
#include <DevboardEvent.h>
#include <DevboardSample.h>
#include <TriggerEvent.h>
#include <TriggerSample.h>
#include <Data.h>
#include <DataRead.h>
#include <DataReadEvio.h>
#include <unistd.h>
using namespace std;

//#define corr1 599
//#define corr1 604
//#define corr2 500

// Process the data
// Pass root file to open as first and only arg.
int main ( int argc, char **argv ) {
    bool debug = false;
    bool flip_channels = true;
    bool mux_channels = false;
    bool skip_corr = false;
    bool read_temp = true;
    int hybrid_type = 0;
    bool evio_format = false;
    bool triggerevent_format = false;
    bool subtract_reference = false;
    int use_fpga = -1;
    int use_hybrid = -1;
    int num_events = -1;
    int ignore_count = 20;
    int c;
    TCanvas         *c1;
    TH2I            *histAll[7];
    short *allSamples[640][7];
    for (int i=0;i<640;i++) for (int j=0;j<7;j++)
    {
        allSamples[i][j] = new short[16384];
        for (int k=0;k<16384;k++) allSamples[i][j][k]=0;
    }
    int chanMap[128];
    for (int idx = 0; idx < 128; idx++ ) {
        chanMap[(32*(idx%4)) + (8*(idx/4)) - (31*(idx/16))] = idx;
    }

    bool channelActive[640];
    double channelValue[640];
    int eventSamples[640][6];
    int channelCount[640];
    double channelMean[640];
    double channelVariance[640];
    double channelCovar[640][640];
    for (int i=0;i<640;i++) {
        channelCount[i] = 0;
        channelMean[i] = 0.0;
        channelVariance[i] = 0.0;
        for (int j=0;j<640;j++) {
            channelCovar[i][j] = 0.0;
        }
    }

    int reference_offset_mean = 0;
    int reference_offset[6];
    int reference_delta[6];

    int apvEventCount = 0;
    double apvEventMean[5];
    double hybridEventMean;
    double apvMean[5];
    double hybridMean = 0.0;
    double apvVariance[5];
    double hybridVariance = 0.0;

    //TH2S *corrHist;
    //corrHist = new TH2S("corrHist","Channel correlation",16384,-0.5,16383.5,16384,-0.5,16383.5);

    double          histMin[640];
    double          histMax[640];
    int hybridMin, hybridMax;
    //TGraph          *sigma;
    int ni = 0;
    double          grChan[640];
    double          grMean[7][640];
    double          grSigma[7][640];
    for (int i=0;i<640;i++) {
        grChan[i]=i;
    }
    DataRead        *dataRead;
    DevboardEvent    event;
    TriggerEvent    triggerevent;
    TriggerSample   *triggersample = new TriggerSample();
    int            eventCount;
    int runCount;
    TString inname;
    TString outdir;
    char            name[200];
    char title[200];
    TGraph          *graph[7];
    TMultiGraph *mg;

    while ((c = getopt(argc,argv,"ho:nmct:H:F:e:EdVS")) !=-1)
        switch (c)
        {
            case 'h':
                printf("-h: print this help\n");
                printf("-d: turn on debug\n");
                printf("-o: use specified output filename\n");
                printf("-n: DAQ (Ryan's) channel numbering\n");
                printf("-m: number channels in raw mux order\n");
                printf("-c: don't compute correlations\n");
                printf("-t: hybrid type (1 for old test run hybrid, 2 for new 2014 hybrid)\n");
                printf("-F: use only specified FPGA\n");
                printf("-H: use only specified hybrid\n");
                printf("-e: stop after specified number of events\n");
                printf("-E: use EVIO file format\n");
                printf("-V: use TriggerEvent event format\n");
                printf("-S: subtract channel 639\n");
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
            case 'c':
                skip_corr = true;
                break;
            case 'm':
                mux_channels = true;
                break;
            case 't':
                hybrid_type = atoi(optarg);
                break;
            case 'F':
                use_fpga = atoi(optarg);
                break;
            case 'H':
                use_hybrid = atoi(optarg);
                break;
            case 'e':
                num_events = atoi(optarg);
                break;
            case 'E':
                evio_format = true;
                break;
            case 'V':
                triggerevent_format = true;
                break;
            case 'S':
                subtract_reference = true;
                break;
            case 'd':
                debug = true;
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

    if (evio_format)
        dataRead = new DataReadEvio();
    else 
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

    // Start X11 view
    // TApplication theApp("App",NULL,NULL);

    // Root file is the first and only arg
    if ( argc-optind != 1 ) {
        cout << "Usage: meeg_baseline data_file\n";
        return(1);
    }


    hybridMin = 16384;
    hybridMax = 0;
    for (int channel=0; channel < 640; channel++) {
        histMin[channel] = 16384;
        histMax[channel] = 0;
    }

    if (inname=="")
    {
        inname=argv[optind];

        inname.ReplaceAll(".bin","");
        if (inname.Contains('/')) {
            inname.Remove(0,inname.Last('/')+1);
        }
    }

    sprintf(name,"%s_baseline.root",inname.Data());
    TFile *myFile = new TFile(name,"RECREATE");

    // 2d histogram
    for (int i=0;i<6;i++)
    {
        sprintf(name,"Value_Hist_s%d",i);
        sprintf(title,"Baseline values, sample %d;ADC counts;Channel",i);
        histAll[i] = new TH2I(name,title,16384,-0.5,16383.5,640,-0.5,639.5);
    }
    sprintf(title,"Baseline values, all samples;ADC counts;Channel");
    histAll[6] = new TH2I("Value_Hist_All",title,16384,-0.5,16383.5,640,-0.5,639.5);


    ofstream outfile;
    cout << "Writing calibration to " << inname+".base" << endl;
    outfile.open(inname+".base");

    cout << "Reading data file " <<argv[optind] << endl;
    // Attempt to open data file
    bool result = dataRead->open(argv[optind]);    
    cout << "DataRead open result: " << result << endl;

    if ( !result ) {
      cout << "Failed to open file" << endl;
      return(2);
    }
    
    cout << "Successfully opened file" << endl;

    cout << "pos " << dataRead->pos() << endl;

    //if ( ! dataRead->open(argv[optind]) ) return(2);

    TString confname=argv[optind];
    confname.ReplaceAll(".bin","");
    confname.Append(".conf");
    if (confname.Contains('/')) {
        confname.Remove(0,confname.Last('/')+1);
    }

    ofstream outconfig;
    cout << "Writing configuration to " <<outdir<<confname << endl;
    outconfig.open(outdir+confname);
    cout << "1 " << endl;
    bool readOK;

    cout << "file position " << dataRead->pos() << endl;

    if (triggerevent_format) {
        dataRead->next(&triggerevent);
    } else {
        dataRead->next(&event);
    }

    cout << "Read config xml " << endl;
    outconfig << dataRead->getConfigXml();
    outconfig << endl;
    outconfig << dataRead->getStatusXml();
    outconfig.close();

    cout << "2 " << endl;

    runCount = atoi(dataRead->getConfig("RunCount").c_str());
    int max_count = runCount==0 ? 10000 : runCount;
    double *apv_means[5];
    for (int i=0;i<5;i++) apv_means[i] = new double[max_count];
    double *moving_ti = new double[max_count];
    /*
       double *moving_yi = new double[runCount];
       double *moving_yi2 = new double[runCount];
       */

    if(debug) printf("%d events expected from config\n",runCount);


    // Process each event
    eventCount = 0;

    do {

        int rce = 0;
        int fpga = 0;
        int samplecount;

        if (triggerevent_format) {
            samplecount = triggerevent.count();
        } else {
            fpga = event.fpgaAddress();
            samplecount = event.count();
            if (read_temp && !event.isTiFrame()) for (uint i=0;i<4;i++) {
                printf("Event %d, temperature #%d: %f\n",eventCount,i,event.temperature(i,hybrid_type==1));
                read_temp = false;
            }
        }

        if(debug) printf("Event %d\n",eventCount);

        //if(debug) printf("fpga %d\n",event.fpgaAddress());
        if (!triggerevent_format && fpga==7) 
        {
            //printf("not a data event\n");
            continue;
        }
        //if(debug) cout<<"  fpga #"<<fpga<<"; number of samples = "<<event.count() << " is" <<endl;
        //if(debug) cout<<"  fpga #"<<fpga<<"; number of samples = "<<event.count()<<endl;
        if (eventCount%1000==0) printf("Event %d\n",eventCount);
        if (num_events!=-1 && eventCount >= num_events) break;
        if (eventCount<max_count) {
            moving_ti[eventCount] = eventCount;
            for (int i=0;i<5;i++) apv_means[i][eventCount] = 0.0;
        }
        for (int i=0;i<640;i++)
        {
            channelActive[i] = false;
        }

        for (int x=0; x < samplecount; x++) {
            int hyb;
            int apv;
            int apvch;
            int channel;
            int samples[6];

            bool goodSample = true;

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
                printf("event %d\tx=%d\tF%d H%d A%d channel %d, samples:\t%d\t%d\t%d\t%d\t%d\t%d\n",eventCount,x,event.fpgaAddress(),sample->hybrid(),sample->apv(),sample->channel(),sample->value(0),sample->value(1),sample->value(2),sample->value(3),sample->value(4),sample->value(5));
            }
            if (use_fpga!=-1 && fpga!=use_fpga) continue;
            if (use_hybrid!=-1 && hyb!=use_hybrid) continue;
            if (!goodSample) continue;
            //printf("hybrid %d\n",sample->hybrid());

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
            }

            if (subtract_reference && apv==0 && apvch==127) {
                if (eventCount==ignore_count-1) {
                    //for (int y=0;y<6;y++) reference_offset[y]=samples[y];
                    for (int y=0;y<6;y++) reference_offset_mean += samples[y];
                    reference_offset_mean /= 6;
                } else if (eventCount >= ignore_count)
                    //for (int y=0;y<6;y++) reference_delta[y]=samples[y]-reference_offset[y];
                    for (int y=0;y<6;y++) reference_delta[y]=samples[y]-reference_offset_mean;
            }

            for ( int y=0; y < 6; y++ ) {
                eventSamples[channel][y] = samples[y];
            }
            if ( eventCount >= ignore_count ) {
                channelCount[channel]++;
                channelActive[channel] = true;
            }
        }

        if (subtract_reference && eventCount>=ignore_count) {
            //int mean_delta = 0;
            //for (int y=0;y<6;y++) mean_delta+=reference_delta[y];
            //mean_delta /= 6;
            for (int i=0;i<640;i++) for ( int y=0; y<6; y++ ) {
                eventSamples[i][y]-=reference_delta[y];
                //eventSamples[i][y]-=mean_delta;
            }
        }

        for (int i=0;i<5;i++) apvEventMean[i] = 0;
        hybridEventMean = 0;
        for (int i=0;i<640;i++) {
            // Filter APVs
            if ( eventCount >= ignore_count ) {
                channelValue[i] = 0;
                for ( int y=0; y < 6; y++ ) {
                    int value = eventSamples[i][y];
                    channelValue[i]+=value;

                    //vhigh = (value << 1) & 0x2AAA;
                    //vlow  = (value >> 1) & 0x1555;
                    //value = vlow | vhigh;

                    histAll[y]->Fill(value,i);
                    histAll[6]->Fill(value,i);
                    allSamples[i][y][value]++;
                    allSamples[i][6][value]++;

                    if ( value < histMin[i] ) histMin[i] = value;
                    if ( value > histMax[i] ) histMax[i] = value;
                    if ( value < hybridMin ) hybridMin = value;
                    if ( value > hybridMax ) hybridMax = value;
                }
                channelValue[i]/=6.0;
                channelValue[i] = eventSamples[i][0];
                //if (channel==17) printf("%d %d %d %d %d %d\n",eventCount,x,rce,fpga,hyb,samples[0]);

            }
            double mean = 0;
            for (int y=0;y<6;y++) mean+=eventSamples[i][y];
            mean/=6.0;
            int apv = i/128;
            if (flip_channels) apv = 4-apv;
            hybridEventMean += eventSamples[i][5];
            apvEventMean[apv] += eventSamples[i][5];
            if (eventCount<max_count) {
                apv_means[apv][eventCount] += mean;
            }
        }
        /*
           if (channel==corr1)
           {
           moving_yi[eventCount] = mean;
           }
           if (channel==corr2)
           {
           moving_yi2[eventCount] = mean;
           }
           */
        if (hybridEventMean!=0) {
            apvEventCount++;
            for (int i=0;i<5;i++) {
                apvEventMean[i] /= 128;
                double delta = apvEventMean[i] - apvMean[i];
                apvMean[i] += delta/apvEventCount;
                apvVariance[i] += delta*(apvEventMean[i]-apvMean[i]);
            }
            hybridEventMean /= 640;
            double delta = hybridEventMean - hybridMean;
            hybridMean += delta/apvEventCount;
            hybridVariance += delta*(hybridEventMean-hybridMean);
        }
        if (eventCount<max_count) {
            for (int i=0;i<5;i++) apv_means[i][eventCount] /= 128;
        }
        if (!skip_corr) {
            //double grDelta[640];
            //ni=0;
            for (int i=0;i<640;i++) if (channelActive[i])
            {
                double delta = channelValue[i]-channelMean[i];
                if (channelCount[i]==1)
                {
                    channelMean[i] = channelValue[i];
                }
                else
                {
                    channelMean[i] += delta/channelCount[i];
                }
                channelVariance[i] += delta*(channelValue[i]-channelMean[i]);
                for (int j=0;j<i;j++) if (channelActive[j])
                {
                    channelCovar[i][j] += (channelValue[i]-channelMean[i])*(channelValue[j]-channelMean[j]);
                    /*
                       if (i==corr1&&j==corr2)
                       {
                    //if (channelValue[i]<7620)
                    corrHist->Fill(channelValue[i],channelValue[j]);
                    //else
                    //printf("event %d\n",eventCount);
                    }
                    */
                }
                /*
                   if (i>255 && i<384) {
                   grChan[ni]=i;
                   grDelta[ni]=delta;
                //printf("%d: %f, %f\n",ni,grChan[ni],grDelta[ni]);
                ni++;
                }
                */
            }
        }
        /*
           int startscope=1000;
           if (eventCount==startscope-1) {
           mg = new TMultiGraph();
           }
           if (ni>0 && eventCount>=startscope && eventCount<startscope+20) {
           printf("%f\n",grDelta[17]);
           int i = eventCount-startscope;
           TGraph* tempgraph;
           tempgraph = new TGraph(ni,grChan,grDelta);
           tempgraph->SetMarkerColor(i+1);
           tempgraph->SetLineColor(i+1);
           mg->Add(tempgraph);
           }
           if (eventCount==startscope+20) {
           c1->Clear();
           sprintf(name,"%s_base_delta.png",inname.Data());
           printf("%s\n",name);
           mg->Draw("ALP");
           mg->GetYaxis()->SetRangeUser(-500,500);
           c1->SaveAs(name);
           for (int i=0;i<7;i++)
           delete graph[i];
           delete mg;
           }
           */

        eventCount++;

        if (triggerevent_format) {
            readOK = dataRead->next(&triggerevent);
        } else {
            readOK = dataRead->next(&event);
        }
    } while (readOK);
    dataRead->close();

    if (eventCount != runCount)
    {
        printf("ERROR: events read = %d, runCount = %d\n",eventCount, runCount);
    }

    /*
       TGraph *movingGraph = new TGraph(eventCount,moving_ti,moving_yi);
       TGraph *movingGraph2 = new TGraph(eventCount,moving_ti,moving_yi2);
       movingGraph->SetMarkerColor(2);
       mg = new TMultiGraph();
       mg->Add(movingGraph);
       mg->Add(movingGraph2);

       mg->Draw("a*");
       c1->SaveAs("meeg.png");
       c1->Clear();
       */

    mg = new TMultiGraph();
    for (int i=0;i<5;i++)
    {
        graph[i] = new TGraph(min(eventCount,max_count),moving_ti,apv_means[i]);
        graph[i]->SetMarkerColor(i+1);
        if (i==4) graph[i]->SetMarkerColor(6);
        mg->Add(graph[i]);
    }
    mg->Draw("ap");
    sprintf(name,"%s_base_apvmeans.png",inname.Data());
    printf("%s\n",name);
    c1->SaveAs(name);
    c1->Clear();
    //for (int i=0;i<5;i++) delete graph[i];
    delete mg;
    for (int i=0;i<5;i++) delete[] apv_means[i];
    delete[] moving_ti;

    for (int i=0;i<5;i++) {
        printf("APV %d: common-mode noise %f\n",i,sqrt(apvVariance[i]/(apvEventCount-1)));
    }
    printf("Hybrid: common-mode noise %f\n",sqrt(hybridVariance/(apvEventCount-1)));

    /*
       c1->Clear();
       corrHist->GetXaxis()->SetRangeUser(histMin[corr1],histMax[corr1]);
       corrHist->GetYaxis()->SetRangeUser(histMin[corr2],histMax[corr2]);
       corrHist->Draw("colz");
       c1->SaveAs("meeg2.png");
       */

    int deadAPV = -1;
    for (int i=0;i<640;i++)
    {
        if (channelCount[i])
        {
            if (i/128==deadAPV)
            {
                printf("Counted %d events on channel %d, even though we thought APV %d was dead\n",channelCount[i],i,deadAPV);
            }
            if (channelCount[i]!=(int)eventCount-20)
            {
                printf("Counted %d events for channel %d; expected %d\n",channelCount[i],i,eventCount-20);
            }
            if (!skip_corr)
            {
                channelVariance[i]/=channelCount[i];
                for (int j=0;j<i;j++) if (channelCount[j])
                {
                    channelCovar[i][j]/=eventCount-20;
                    channelCovar[i][j]/=sqrt(channelVariance[i]);
                    channelCovar[i][j]/=sqrt(channelVariance[j]);
                    //printf("%d, %d, %f\n",i,j,channelCovar[i][j]);
                }
            }
        }
        else
        {
            if (i/128!=deadAPV)
            {
                deadAPV = i/128;
                printf("No events on channel %d; assuming APV %d is dead\n",i,deadAPV);
            }
        }
    }

    ni=0;
    for (int channel = 0; channel < 640; channel++) {
        int count;
        if (histMin[channel]<=histMax[channel])
        {
            grChan[ni]  = channel;
            outfile <<channel<<"\t";
            for (int i=0;i<7;i++)
            {
                doStats_mean(16384,(int)histMin[channel],(int)histMax[channel],allSamples[channel][i],count,grMean[i][ni],grSigma[i][ni]);
                outfile<<grMean[i][ni]<<"\t"<<grSigma[i][ni]<<"\t";
                //printf("%d:\t%f\t%f\n",channel,grSigma[i][ni],sqrt(channelVariance[channel]));
            }
            outfile<<endl;
            ni++;
        }
    }


    if (!skip_corr)
    {
        //c1->SetLogz();
        TH2F *channelCorr = new TH2F("channelCorr","Channel correlation - positive correlations;channel 1;channel 2",640,-0.5,639.5,640,-0.5,639.5);
        channelCorr->SetStats(kFALSE);
        channelCorr->SetNdivisions(0,"XY");
        for (int i=0;i<640;i++) for (int j=0;j<i;j++) if (channelCovar[i][j]>0)
        {
            channelCorr->Fill(i,j,channelCovar[i][j]);
            channelCorr->Fill(j,i,channelCovar[i][j]);
        }
        channelCorr->Draw("colz");
        sprintf(name,"%s_pos_corr.pdf",inname.Data());
        c1->SaveAs(name);
        sprintf(name,"%s_pos_corr.png",inname.Data());
        c1->SaveAs(name);

        channelCorr->Reset();
        channelCorr->SetTitle("Channel correlation - negative correlations;channel 1;channel 2");
        for (int i=0;i<640;i++) for (int j=0;j<i;j++) if (channelCovar[i][j]<0)
        {
            channelCorr->Fill(i,j,-1*channelCovar[i][j]);
            channelCorr->Fill(j,i,-1*channelCovar[i][j]);
        }
        channelCorr->Draw("colz");
        sprintf(name,"%s_neg_corr.pdf",inname.Data());
        c1->SaveAs(name);
        sprintf(name,"%s_neg_corr.png",inname.Data());
        c1->SaveAs(name);
        //c1->SetLogz(0);

        c1->Clear();
        mg = new TMultiGraph();
        double yi[6][640], xi[6][640];
        for (int i=0;i<5;i++) 
        {
            for (int j=0;j<127;j++)
            {
                yi[i][j] = channelCovar[i*128+j+1][i*128+j];
                xi[i][j] = i*128+j+0.5;
            }
            graph[i] = new TGraph(127,xi[i],yi[i]);
            graph[i]->SetMarkerColor(i+2);
            mg->Add(graph[i]);
        }
        for (int i=1;i<5;i++)
        {
            yi[5][i-1] = channelCovar[i*128][i*128-1];
            xi[5][i-1] = i*128-0.5;
        }
        graph[5] = new TGraph(4,xi[5],yi[5]);
        graph[5]->SetMarkerColor(1);
        mg->Add(graph[5]);

        mg->SetTitle("Adjacent channel correlations");
        mg->Draw("a*");
        //mg->GetXaxis()->SetRangeUser(0,640);
        sprintf(name,"%s_corr_adj.png",inname.Data());
        c1->SaveAs(name);
        for (int i=0;i<6;i++) delete graph[i];
        delete mg;

        //printf("correlation %f %f\n",corrHist->GetCorrelationFactor(),channelCovar[corr1][corr2]);
        //printf("correlation %f\n",channelCovar[corr1][corr2]);
    }

    histAll[6]->GetXaxis()->SetRangeUser(hybridMin,hybridMax);
    histAll[6]->Draw("colz");
    sprintf(name,"%s_base.png",inname.Data());
    c1->SaveAs(name);

    for (int i=0;i<6;i++)
    {
        histAll[i]->GetXaxis()->SetRangeUser(hybridMin,hybridMax);
        histAll[i]->Draw("colz");
        sprintf(name,"%s_base_%d.png",inname.Data(),i);
        c1->SaveAs(name);
    }

    c1->Clear();
    mg = new TMultiGraph();
    for (int i=0;i<7;i++)
    {
        graph[i] = new TGraph(ni,grChan,grMean[i]);
        graph[i]->SetMarkerColor(i+2);
        mg->Add(graph[i]);
    }
    graph[6]->SetMarkerColor(1);
    mg->SetTitle("Pedestal;Channel;ADC counts");
    //mg->GetXaxis()->SetRangeUser(0,640);
    mg->Draw("a*");
    sprintf(name,"%s_base_pedestal.png",inname.Data());
    c1->SaveAs(name);
    for (int i=0;i<7;i++)
        delete graph[i];
    delete mg;

    c1->Clear();
    mg = new TMultiGraph();
    double max = 0;
    for (int i=0;i<7;i++)
        for (int j=0;j<ni;j++)
            if (grSigma[i][j]>max) max = grSigma[i][j];

    for (int i=0;i<7;i++)
    {
        graph[i] = new TGraph(ni,grChan,grSigma[i]);
        graph[i]->SetMarkerColor(i+2);
        mg->Add(graph[i]);

        if(i==0) {
            for (int ch=0;ch<graph[i]->GetN();ch++) {
                double x,y;
                graph[i]->GetPoint(ch,x,y);
                if (y > 50) printf("ch %f: noise %f above 50\n",x,y);
                if (y < 20) printf("ch %f: noise %f below 20\n",x,y);
            }
        }

    }
    TGraph *CMnoise[5];
    for (int i=0;i<5;i++)
    {
        int apv = flip_channels?4-i:i;
        double chan[2] = {128.0*apv,128.0*(apv+1)-1};
        double noise = sqrt(apvVariance[apv]/(apvEventCount-1));
        double yi[2] = {noise,noise};
        CMnoise[i] = new TGraph(2,chan,yi);
        CMnoise[i]->SetLineColor(2);
        mg->Add(CMnoise[i],"L");
    }
    TGraph *hybridNoise;
    {
        double chan[2] = {0,639};
        double noise = sqrt(hybridVariance/(apvEventCount-1));
        double yi[2] = {noise,noise};
        hybridNoise = new TGraph(2,chan,yi);
        mg->Add(hybridNoise,"L");
    }
    graph[6]->SetMarkerColor(1);
    mg->SetTitle("Noise;Channel;ADC counts");
    mg->Draw("a*");
    mg->GetYaxis()->SetRangeUser(0,1.2*max);
    sprintf(name,"%s_base_noise.png",inname.Data());
    c1->SaveAs(name);
    for (int i=0;i<7;i++)
        delete graph[i];
    for (int i=0;i<5;i++)
        delete CMnoise[i];
    delete hybridNoise;
    delete mg;

    c1->Clear();
    mg = new TMultiGraph();
    double grShift[640];
    for (int i=0;i<6;i++)
    {
        for (int n=0;n<ni;n++) grShift[n]=grMean[i][n]-grMean[6][n];
        graph[i] = new TGraph(ni,grChan,grShift);
        graph[i]->SetMarkerColor(i+1);
        mg->Add(graph[i]);
    }
    mg->SetTitle("Sample-to-sample shift;Channel;ADC counts");
    //mg->GetXaxis()->SetRangeUser(0,640);
    mg->Draw("a*");
    sprintf(name,"%s_base_shift.png",inname.Data());
    c1->SaveAs(name);
    for (int i=0;i<6;i++)
        delete graph[i];
    delete mg;

    // Start X-Windows
    //theApp.Run();

    // Close file
    outfile.close();
    delete dataRead;
    for (int i=0;i<640;i++) for (int j=0;j<7;j++)
    {
        delete[] allSamples[i][j];
        //for (int k=0;k<16384;k++) allSamples[i][j][k]=0;
    }
    myFile->Write();
    myFile->Close();
    return(0);
}

