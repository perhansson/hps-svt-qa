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
#include <TH2S.h>
#include <TH2D.h>
#include <TF1.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TApplication.h>
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
#include <TMath.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <unistd.h>
#include "meeg_utils.hh"

#define N_TIME_CONSTS 2

using namespace std;

// Process the data
// Pass root file to open as first and only arg.
int main ( int argc, char **argv ) {
    int c;
    bool plot_tp_fits = false;
    bool plot_fit_results = false;
    bool force_cal_grp = false;
    bool use_baseline_cal = false;
    bool flip_channels = true;
    bool move_fitstart = false;
    bool read_temp = true;
    int hybrid_type = 0;
    bool evio_format = false;
    bool triggerevent_format = false;
    int use_fpga = -1;
    int use_hybrid = -1;
    int num_events = -1;
    double fit_shift;
    ifstream calfile;
    TString inname = "";
    TString outdir = "";
    int cal_grp = -1;
    int cal_delay = 0;
    double delay_step = SAMPLE_INTERVAL/8;
    TCanvas         *c1;
    //TH2I            *histAll;
    short *allSamples[2][640][48] = {{{NULL}}};
    //bool hasSamples[2][640][48] = {{{false}}};
    //TH1D            *histSamples1D;
    int          histMin[640];
    for (int i=0;i<640;i++) histMin[i]=16384;
    int          histMax[640];
    for (int i=0;i<640;i++) histMax[i]=0;
    //TGraph          *mean;
    //TGraph          *sigma;
    DataRead        *dataRead;
    DevboardEvent    event;
    TriggerEvent    triggerevent;
    TriggerSample   *triggersample = new TriggerSample();
    int		samples[6];
    int            eventCount;
    int runCount;
    double          sum;
    char            name[100];
    char            name2[100];
    char title[200];
    int nChan[2] = {0};
    double grChan[2][640];
    double grTp[N_TIME_CONSTS][2][640];
    double grA[2][640];
    double grT0[2][640];
    double grChisq[2][640];
    double          calMean[640][7] = {{0.0}};
    double          calSigma[640][7] = {{1.0}};
    for (int i=0;i<640;i++) {
        for (int j=0;j<7;j++) {
            calMean[i][j] = 0.0;
            calSigma[i][j] = 1.0;
        }
    }

    while ((c = getopt(argc,argv,"hfrg:o:b:d:s:nt:H:F:e:EV")) !=-1)
        switch (c)
        {
            case 'h':
                printf("-h: print this help\n");
                printf("-f: plot Tp fits for each channel\n");
                printf("-g: force use of specified cal group\n");
                printf("-r: plot fit results\n");
                printf("-o: use specified output filename\n");
                printf("-b: use specified baseline cal file\n");
                printf("-d: use specified dtrig baseline cal file\n");
                printf("-n: DAQ (Ryan's) channel numbering\n");
                printf("-s: start fit at given delay after a first guess at T0\n");
                printf("-t: hybrid type (1 for old test run hybrid, 2 for new 2014 hybrid)\n");
                printf("-F: use only specified FPGA\n");
                printf("-H: use only specified hybrid\n");
                printf("-e: stop after specified number of events\n");
                printf("-E: use EVIO file format\n");
                printf("-V: use TriggerEvent event format\n");
                return(0);
                break;
            case 'f':
                plot_tp_fits = true;
                break;
            case 'r':
                plot_fit_results = true;
                break;
            case 'n':
                flip_channels = false;
                break;
            case 'g':
                force_cal_grp = true;
                cal_grp = atoi(optarg);
                break;
            case 'o':
                inname = optarg;
                outdir = optarg;
                if (outdir.Contains('/')) {
                    outdir.Remove(outdir.Last('/')+1);
                }
                else outdir="";
                break;
            case 'b':
                use_baseline_cal = true;
                cout << "Reading baseline calibration from " << optarg << endl;
                calfile.open(optarg);
                while (!calfile.eof()) {
                    int channel;
                    calfile >> channel;
                    for (int i=0;i<7;i++)
                    {
                        double temp;
                        calfile >> temp;
                        calMean[channel][i]+=temp;
                        calfile >> calSigma[channel][i];
                    }
                }
                calfile.close();
                break;
            case 'd':
                cout << "Reading dtrig baseline calibration from " << optarg << endl;
                calfile.open(optarg);
                while (!calfile.eof()) {
                    int channel;
                    calfile >> channel;
                    double temp;
                    for (int i=0;i<6;i++)
                    {
                        calfile >> temp;
                        calMean[channel][i]-=temp;
                        calfile >> calSigma[channel][i];
                    }
                    calfile >> temp;
                    for (int i=0;i<6;i++)
                    {
                        calMean[channel][i]+=temp;
                    }
                    calfile >> calSigma[channel][6];
                }
                calfile.close();
                break;
            case 's':
                move_fitstart = true;
                fit_shift = atof(optarg);
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
    gStyle->SetPalette(1,0);
    gStyle->SetOptStat("emrou");
    gStyle->SetStatW(0.2);                
    gStyle->SetStatH(0.1);                
    gStyle->SetTitleOffset(1.4,"y");
    gStyle->SetPadLeftMargin(0.15);
    c1 = new TCanvas("c1","c1",1200,900);

    // Start X11 view
    //TApplication theApp("App",NULL,NULL);

    // Root file is the first and only arg
    if ( argc-optind==0) {
        cout << "Usage: meeg_tp data_file\n";
        return(1);
    }


    if (inname == "")
    {
        inname=argv[optind];

        inname.ReplaceAll(".bin","");
        if (inname.Contains('/')) {
            inname.Remove(0,inname.Last('/')+1);
        }
    }

    sprintf(name,"%s_tp.root",inname.Data());
    TFile *myFile = new TFile(name,"RECREATE");

    ofstream tpfile[2];
    cout << "Writing Tp calibration to " << inname+".tp_pos" << endl;
    tpfile[0].open(inname+".tp_pos");
    cout << "Writing Tp calibration to " << inname+".tp_neg" << endl;
    tpfile[1].open(inname+".tp_neg");

    ofstream shapefile[2];
    cout << "Writing pulse shape to " << inname+".shape_pos" << endl;
    shapefile[0].open(inname+".shape_pos");
    cout << "Writing pulse shape to " << inname+".shape_neg" << endl;
    shapefile[1].open(inname+".shape_neg");

    ofstream noisefile[2];
    cout << "Writing pulse noise to " << inname+".noise_pos" << endl;
    noisefile[0].open(inname+".noise_pos");
    cout << "Writing pulse noise to " << inname+".noise_neg" << endl;
    noisefile[1].open(inname+".noise_neg");


    while (optind<argc)
    {
        cout << "Reading data file " <<argv[optind] << endl;
        // Attempt to open data file
        if ( ! dataRead->open(argv[optind]) ) return(2);

        bool readOK;

        if (triggerevent_format) {
            dataRead->next(&triggerevent);
        } else {
            dataRead->next(&event);
        }

        if (!evio_format) {
            TString confname=argv[optind];
            confname.ReplaceAll(".bin","");
            confname.Append(".conf");
            if (confname.Contains('/')) {
                confname.Remove(0,confname.Last('/')+1);
            }

            ofstream outconfig;
            cout << "Writing configuration to " <<outdir<<confname << endl;
            outconfig.open(outdir+confname);

            outconfig << dataRead->getConfigXml();
            outconfig << endl;
            outconfig << dataRead->getStatusXml();
            outconfig.close();

            runCount = atoi(dataRead->getConfig("RunCount").c_str());

            if (!force_cal_grp)
            {
                string cgrp = dataRead->getConfig("cntrlFpga:hybrid:apv25:CalGroup");
                if (cgrp.length()==0) cgrp = dataRead->getConfig("FrontEndTestFpga:FebCore:Hybrid:apv25:CalGroup");
                cal_grp = atoi(cgrp.c_str());
                cout<<"Read calibration group "<<cal_grp<<" from data file"<<endl;
            }

            string csel = dataRead->getConfig("cntrlFpga:hybrid:apv25:Csel");
            if (csel.length()==0) csel = dataRead->getConfig("FrontEndTestFpga:FebCore:Hybrid:apv25:Csel");
            cal_delay = atoi(csel.substr(4,1).c_str());
            cout<<"Read calibration delay "<<cal_delay<<" from data file"<<endl;
            if (cal_delay==0)
            {
                cal_delay=8;
                cout<<"Force cal_delay=8 to keep sample time in range"<<endl;
            }
        } else
        {
            runCount = 0;
            if (!force_cal_grp) cal_grp = 0;
            cal_delay = 1;
        }


        // Process each event
        eventCount = 0;
        //bool goodEvent;
        int pulsePolarity = -1;

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

            //if(debug) printf("Event %d\n",eventCount);

            if (!triggerevent_format && fpga==7) 
            {
                //printf("not a data event\n");
                continue;
            }

            //goodEvent = true;
            //for (x=0; x < event.count(); x++) {
            //	sample  = event.sample(x);
            //	channel = (sample->apv() * 128) + sample->channel();
            //	if (channel==32 && sample->value(1)>8482 && sample->value(1)>7800) goodEvent = false; 
            //}
            //if (eventCount%2==0) printf("Event %d is %s\n",eventCount,goodEvent?"good":"bad");
            //goodEvent = true;
            if (eventCount%1000==0) printf("Event %d\n",eventCount);
            if (num_events!=-1 && eventCount >= num_events) break;
            //if (goodEvent) 
            for (int x=0; x < samplecount; x++) {
                int hyb;
                int apv;
                int apvch;
                int channel;

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
                }
                //printf("event %d\tx=%d\tF%d H%d A%d channel %d, samples:\t%d\t%d\t%d\t%d\t%d\t%d\n",eventCount,x,event.fpgaAddress(),sample->hybrid(),sample->apv(),sample->channel(),sample->value(0),sample->value(1),sample->value(2),sample->value(3),sample->value(4),sample->value(5));
                if (use_fpga!=-1 && fpga!=use_fpga) continue;
                if (use_hybrid!=-1 && hyb!=use_hybrid) continue;
                if (!goodSample) continue;

                channel = apvch;

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

                if ((apvch-cal_grp)%8!=0) continue;

                // Filter APVs
                if ( eventCount >= 20 ) {
                    bool bad_event = false;
                    for ( int y=0; y < 6; y++ ) if (samples[y]==0) {
                        printf("sample is zero: event %d, channel %d, sample %d\n",eventCount,channel,y);
                        bad_event = true;
                    }
                    if (bad_event) continue;

                    sum = 0;
                    for ( int y=0; y < 6; y++ ) {
                        //vhigh = (value << 1) & 0x2AAA;
                        //vlow  = (value >> 1) & 0x1555;
                        //value = vlow | vhigh;

                        //histAll->Fill(value,channel);
                        //histSng[channel]->Fill(value);

                        if ( samples[y] < histMin[channel] ) histMin[channel] = samples[y];
                        if ( samples[y] > histMax[channel] ) histMax[channel] = samples[y];
                        sum+=samples[y];
                    }
                    sum-=6*samples[0];
                    int sgn = sum>0?0:1;
                    if (pulsePolarity==-1)
                    {
                        pulsePolarity=((sgn-eventCount)%2 + 2)%2;
                        printf("Saw a %s pulse, event %d, channel %d\n",sgn?"positive":"negative",eventCount,channel);
                    }
                    //int sgn = eventCount%2;
                    for ( int y=0; y < 6; y++ ) {
                        if (allSamples[sgn][channel][8*y+8-cal_delay]==NULL)
                        {
                            allSamples[sgn][channel][8*y+8-cal_delay] = new short[16384];
                            for (int i=0;i<16384;i++) allSamples[sgn][channel][8*y+8-cal_delay][i]=0;
                        }
                        allSamples[sgn][channel][8*y+8-cal_delay][samples[y]]++;
                    }
                    //tpfile<<"T0 " << fit_par[0] <<", A " << fit_par[1] << "Fit chisq " << chisq << ", DOF " << dof << ", prob " << TMath::Prob(chisq,dof) << endl;
                }
            }
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
        optind++;
    }

    double yi[48], ey[48], ti[48];
    int ni;
    TGraphErrors *fitcurve;
    /*
    TF1 *shapingFunction = new TF1("Shaping Function",
            "[0]+[1]*(x>[2])*((x-[2])/[3])*exp(1-((x-[2])/[3]))",-1.0*SAMPLE_INTERVAL,5*SAMPLE_INTERVAL);
            */

    //TF1 *shapingFunction = new TF1("Shaping Function",
            //"[0]+\
            //[1]*(x>[2])*(\
                //([3]/(([3]-[4])*([3]-[5])))*exp(([2]-x)/[3])+\
                //([4]/(([4]-[5])*([4]-[3])))*exp(([2]-x)/[4])+\
                //([5]/(([5]-[3])*([5]-[4])))*exp(([2]-x)/[5]))",
            //-1.0*SAMPLE_INTERVAL,5*SAMPLE_INTERVAL);

    //TF1 *shapingFunction = new TF1("Shaping Function",
            //"[0]+\
            //[1]*(x>[2])*\
            //([3]/(([3]-[4])*([3]-[4])))*(\
                //exp(([2]-x)/[3])-\
                //(1+(([3]-[4])/([3]*[4]))*(x-[2]))*exp(([2]-x)/[4]))",
            //-1.0*SAMPLE_INTERVAL,5*SAMPLE_INTERVAL);

    //TF1 *shapingFunction = new TF1("Shaping Function",
            //"[0]+\
            //[1]*(x>[2])*(\
                //([3]*[3]/(([3]-[4])*([3]-[5])*([3]-[6])))*exp(([2]-x)/[3])+\
                //([4]*[4]/(([4]-[5])*([4]-[6])*([4]-[3])))*exp(([2]-x)/[4])+\
                //([5]*[5]/(([5]-[6])*([5]-[3])*([5]-[4])))*exp(([2]-x)/[5])+\
                //([6]*[6]/(([6]-[3])*([6]-[4])*([6]-[5])))*exp(([2]-x)/[6]))",
            //-1.0*SAMPLE_INTERVAL,5*SAMPLE_INTERVAL);

    //TF1 *shapingFunction = new TF1("Shaping Function",
            //"[0]+\
            //[1]*(x>[2])*(\
                //([3]*[3]/(([3]-[4])*([3]-[5])*([3]-[5])))*exp(([2]-x)/[3])+\
                //([4]*[4]/(([4]-[5])*([4]-[5])*([4]-[3])))*exp(([2]-x)/[4])+\
                //(([3]*[3]/(([4]-[3])*([3]-[5])*([3]-[5]))) - ([4]*[4]/(([4]-[3])*([4]-[5])*([4]-[5]))))*exp(([2]-x)/[5])+\
                //((x-[2])/(([5]-[3])*([5]-[4])))*exp(([2]-x)/[5]))",
            //-1.0*SAMPLE_INTERVAL,5*SAMPLE_INTERVAL);

    TF1 *shapingFunction = new TF1("Shaping Function",
            "[0]+\
            [1]*(x>[2])*\
            ([3]*[3]/(([3]-[4])*([3]-[4])*([3]-[4])))*(\
                exp(([2]-x)/[3])-\
                (1+\
                 (([3]-[4])/([3]*[4]))*(x-[2])+\
                 (([3]-[4])*([3]-[4])/(2*[3]*[4]*[3]*[4]))*(x-[2])*(x-[2]))*exp(([2]-x)/[4]))",
            -1.0*SAMPLE_INTERVAL,5*SAMPLE_INTERVAL);
            
            
    double chanNoise[2][640]={{0.0}};
    double chanChan[640];
    //TH1S *histSamples1D = new TH1S("h1","h1",16384,-0.5,16383.5);
    TH2S *histSamples;
    double A, T0, A0, fit_start;

    for (int i=0;i<640;i++) chanChan[i] = i;

    for (int channel=0;channel<640;channel++) for (int sgn=0;sgn<2;sgn++) {
        ni=0;
        for (int i=0;i<48;i++) if (allSamples[sgn][channel][i]!=NULL) 
        {
            int nsamples = 0;
            double rms;
            doStats(16384, histMin[channel], histMax[channel], allSamples[sgn][channel][i], nsamples, yi[ni], rms);
            if (use_baseline_cal) yi[ni] -= calMean[channel][i/8];
            //if (use_baseline_cal) yi[ni] += calMean[channel][i/8]-2*calMean[channel][6];

            ey[ni] = rms/sqrt((double)nsamples);
            chanNoise[sgn][channel]+=rms;
            //histSamples1D->Reset();
            //printf("TH!: %f, %f\n", histSamples1D->GetEntries(),(histSamples1D->GetRMS()*histSamples1D->GetRMS()-yi[ni]*yi[ni]));
            //for (int j=histMin[channel];j<=histMax[channel];j++)
            //{
            //histSamples1D->SetBinContent(j+1,allSamples[sgn][channel][i][j]);
            //}
            //if (histSamples1D->Fit("gaus","Q0")==0) {
            //yi[ni]  = histSamples1D->GetFunction("gaus")->GetParameter(1);
            //ey[ni]  = histSamples1D->GetFunction("gaus")->GetParError(1);
            //}
            //else {
            //	printf("Could not fit channel %d, polarity %d, sample %d\n",channel,sgn,i);
            //}

            ti[ni] = (i-8)*delay_step;
            ni++;
        }
        if (ni==0) continue;
        chanNoise[sgn][channel]/=ni;

        if (plot_tp_fits)
        {
            gStyle->SetOptStat(0);
            sprintf(name,"samples_%s_%i",sgn?"neg":"pos",channel);
            sprintf(title,"APV25 pulse shape, channel %d, %s pulses;Time [ns];Amplitude [ADC counts]",channel,sgn?"negative":"positive");
            if (use_baseline_cal)
                histSamples = new TH2S(name,title,48,-8.5*delay_step,39.5*delay_step,16384,-0.5-calMean[channel][6],16383.5-calMean[channel][6]);
            else 
                histSamples = new TH2S(name,title,48,-8.5*delay_step,39.5*delay_step,16384,-0.5,16383.5);
            c1->Clear();
            for (int i=0;i<48;i++) if (allSamples[sgn][channel][i]!=NULL) 
            {
                for (int j=histMin[channel];j<=histMax[channel];j++)
                {
                    if (use_baseline_cal)
                        histSamples->Fill((i-8)*delay_step,j-calMean[channel][i/8],allSamples[sgn][channel][i][j]);
                    //histSamples->Fill((i-8)*delay_step,j+calMean[channel][i/8]-2*calMean[channel][6],allSamples[sgn][channel][i][j]);
                    else
                        histSamples->Fill((i-8)*delay_step,j,allSamples[sgn][channel][i][j]);
                }
            }
            if (use_baseline_cal)
                histSamples->GetYaxis()->SetRangeUser(histMin[channel]-calMean[channel][6],histMax[channel]-calMean[channel][6]);
            else
                histSamples->GetYaxis()->SetRangeUser(histMin[channel],histMax[channel]);
            histSamples->Draw("colz");
        }

        fitcurve = new TGraphErrors(ni,ti,yi,NULL,ey);
        if (sgn==0) shapingFunction->SetParameter(1,TMath::MaxElement(ni,yi)-yi[0]);
        else shapingFunction->SetParameter(1,TMath::MinElement(ni,yi)-yi[0]);
        shapingFunction->SetParameter(2,-10.0);
        //shapingFunction->SetParameter(3,50.0);
        shapingFunction->SetParameter(3,80.0);
        shapingFunction->SetParameter(4,12.0);
        //shapingFunction->SetParameter(5,10.0);
        //shapingFunction->SetParameter(6,70.0);
        shapingFunction->FixParameter(0,yi[0]);
        A0 = yi[0];
        if (ni>0)
        {
            if (fitcurve->Fit(shapingFunction,"Q0","",-1*SAMPLE_INTERVAL,5*SAMPLE_INTERVAL)==0)
            {
                A = shapingFunction->GetParameter(1);
                T0 = shapingFunction->GetParameter(2);
                for (int i=0;i<N_TIME_CONSTS;i++) {
                    grTp[i][sgn][nChan[sgn]]=shapingFunction->GetParameter(3+i);
                }
                //printf("%f, %f, %f, %f\n",shapingFunction->GetParameter(0),shapingFunction->GetParameter(1),shapingFunction->GetParameter(2),shapingFunction->GetParameter(3));
                //printf("%f, %f, %f, %f, %f, %f\n",shapingFunction->GetParameter(0),shapingFunction->GetParameter(1),shapingFunction->GetParameter(2),shapingFunction->GetParameter(3),shapingFunction->GetParameter(4),shapingFunction->GetParameter(5));
                //printf("%f, %f, %f, %f, %f\n",shapingFunction->GetParameter(0),shapingFunction->GetParameter(1),shapingFunction->GetParameter(2),shapingFunction->GetParameter(3),shapingFunction->GetParameter(4));
                //printf("%f, %f, %f, %f, %f, %f, %f\n",shapingFunction->GetParameter(0),shapingFunction->GetParameter(1),shapingFunction->GetParameter(2),shapingFunction->GetParameter(3),shapingFunction->GetParameter(4),shapingFunction->GetParameter(5),shapingFunction->GetParameter(6));
                double residuals[8];
                //for (int i =0;i<8;i++) {
                //    residuals[i]=0;
                //}
                for (int i =0;i<8;i++) {
                    double dataX,dataY;
                    fitcurve->GetPoint(i,dataX,dataY);
                    double res = dataY - shapingFunction->Eval(dataX);
                    residuals[i] = res;
                    //residuals[i] += res/6.0;
                    //printf("%f, %f, %f\n",dataX,dataY,res);
                }
                for (int i =0;i<48;i++) {
                    double dataX,dataY;
                    fitcurve->GetPoint(i,dataX,dataY);
                    fitcurve->SetPoint(i,dataX,dataY-residuals[i%8]);
                }
                fitcurve->Fit(shapingFunction,"Q0","",-1*SAMPLE_INTERVAL,5*SAMPLE_INTERVAL);
                if (move_fitstart)
                {
                    fit_start = T0+fit_shift;
                    fitcurve->Fit(shapingFunction,"Q0","",fit_start,5*SAMPLE_INTERVAL);
                    A = shapingFunction->GetParameter(1);
                    T0 = shapingFunction->GetParameter(2);
                    for (int i=0;i<N_TIME_CONSTS;i++) {
                        grTp[i][sgn][nChan[sgn]]=shapingFunction->GetParameter(3+i);
                    }
                }
                printf("%f, %f, %f",shapingFunction->GetParameter(0),shapingFunction->GetParameter(1),shapingFunction->GetParameter(2));
                for (int i=0;i<N_TIME_CONSTS;i++) {
                    printf(", %f",shapingFunction->GetParameter(3+i));
                }
                printf("\n");
                grChan[sgn][nChan[sgn]]=channel;
                grA[sgn][nChan[sgn]]=A;
                if (sgn==1) grA[sgn][nChan[sgn]]*=-1;
                grT0[sgn][nChan[sgn]]=T0;
                grChisq[sgn][nChan[sgn]]=shapingFunction->GetChisquare();
                nChan[sgn]++;
            }
            else
            {
                printf("Could not fit pulse shape for channel %d, polarity %d\n",channel,sgn);
            }
        }
        if (plot_tp_fits)
        {
            fitcurve->SetLineWidth(3);
            fitcurve->Draw();
            if (move_fitstart)
            {
                shapingFunction->SetLineStyle(1);
                shapingFunction->SetLineWidth(3);
                shapingFunction->SetLineColor(2);
                shapingFunction->SetRange(fit_start,5*SAMPLE_INTERVAL);
                shapingFunction->DrawCopy("LSAME");
                shapingFunction->SetRange(-1*SAMPLE_INTERVAL,fit_start);
                shapingFunction->SetLineStyle(2);
                shapingFunction->Draw("LSAME");
            }
            else
            {
                shapingFunction->SetLineStyle(1);
                shapingFunction->SetLineWidth(3);
                shapingFunction->SetLineColor(2);
                shapingFunction->SetRange(-1*SAMPLE_INTERVAL,5*SAMPLE_INTERVAL);
                shapingFunction->Draw("LSAME");
            }
            sprintf(name,"%s_tp_fit_%s_%i.png",inname.Data(),sgn?"neg":"pos",channel);
            c1->SaveAs(name);
            delete histSamples;
            gStyle->SetOptStat("emrou");
        }
        delete fitcurve;

        shapefile[sgn]<<channel<<"\t"<<ni;
        for (int i=0;i<ni;i++)
        {
            shapefile[sgn]<<"\t"<<ti[i]-T0<<"\t"<<(yi[i]-A0)/A<<"\t"<<ey[i]/A;
        }
        shapefile[sgn]<<endl;
    }
    for (int sgn=0;sgn<2;sgn++)
    {
        for (int i=0;i<nChan[sgn];i++)
        {
            tpfile[sgn] <<grChan[sgn][i]<<"\t"<<grA[sgn][i]<<"\t"<<grT0[sgn][i]<<"\t";
            for (int j=0;j<N_TIME_CONSTS;j++) {
                tpfile[sgn] <<grTp[j][sgn][i]<<"\t";
            }
            tpfile[sgn] <<grChisq[sgn][i]<<endl;
        }
        for (int i=0;i<640;i++) noisefile[sgn]<<i<<"\t"<<chanNoise[sgn][i]<<endl;
    }

    if (plot_fit_results) for (int sgn=0;sgn<2;sgn++)
    {
        c1->SetLogy(0);
        sprintf(name,"A_%s",sgn?"neg":"pos");
        sprintf(name2,"%s_tp_A_%s.png",inname.Data(),sgn?"neg":"pos");
        sprintf(title,"Fitted amplitude, %s pulses;Channel;Amplitude [ADC counts]",sgn?"negative":"positive");
        plotResults(title, name, name2, nChan[sgn], grChan[sgn], grA[sgn], c1);

        c1->SetLogy(0);
        sprintf(name,"T0_%s",sgn?"neg":"pos");
        sprintf(name2,"%s_tp_T0_%s.png",inname.Data(),sgn?"neg":"pos");
        sprintf(title,"Fitted T0, %s pulses;Channel;T0 [ns]",sgn?"negative":"positive");
        plotResults(title, name, name2, nChan[sgn], grChan[sgn], grT0[sgn], c1);

        for (int j=0;j<N_TIME_CONSTS;j++) {
            c1->SetLogy(0);
            sprintf(name,"Tp%d_%s",j+1,sgn?"neg":"pos");
            sprintf(name2,"%s_tp_Tp%d_%s.png",inname.Data(),j+1,sgn?"neg":"pos");
            sprintf(title,"Fitted Tp%d, %s pulses;Channel;Tp [ns]",j+1,sgn?"negative":"positive");
            plotResults(title, name, name2, nChan[sgn], grChan[sgn], grTp[j][sgn], c1);
        }

        c1->SetLogy(0);
        sprintf(name,"Chisq_%s",sgn?"neg":"pos");
        sprintf(name2,"%s_tp_Chisq_%s.png",inname.Data(),sgn?"neg":"pos");
        sprintf(title,"Fit chisq, %s pulses;Channel;Chisq",sgn?"negative":"positive");
        plotResults(title, name, name2, nChan[sgn], grChan[sgn], grChisq[sgn], c1);

        c1->SetLogy(0);
        sprintf(name,"Noise_%s",sgn?"neg":"pos");
        sprintf(name2,"%s_tp_Noise_%s.png",inname.Data(),sgn?"neg":"pos");
        sprintf(title,"Mean RMS noise per sample, %s pulses;Channel;Noise [ADC counts]",sgn?"negative":"positive");
        plotResults(title, name, name2, 640, chanChan, chanNoise[sgn], c1);

    }

    // Start X-Windows
    //theApp.Run();

    // Close file
    tpfile[0].close();
    tpfile[1].close();
    shapefile[0].close();
    shapefile[1].close();
    noisefile[0].close();
    noisefile[1].close();
    myFile->Write();
    myFile->Close();
    return(0);
}

