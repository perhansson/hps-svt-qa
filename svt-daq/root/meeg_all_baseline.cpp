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
#include <TH2F.h>
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

#define MAX_RCE 14
#define MAX_FEB 10
#define MAX_HYB 4

//#define corr1 599
//#define corr1 604
//#define corr2 500

// Process the data
// Pass root file to open as first and only arg.
int main ( int argc, char **argv ) {
    bool flip_channels = true;
    bool mux_channels = false;
    bool read_temp = true;
    int hybrid_type = 0;
    bool evio_format = false;
    bool triggerevent_format = false;
    int use_fpga = -1;
    int use_hybrid = -1;
    int num_events = -1;
    double threshold_sigma = 2.0;
    int c;
    TCanvas         *c1;
    int chanMap[128];
    for (int idx = 0; idx < 128; idx++ ) {
        chanMap[(32*(idx%4)) + (8*(idx/4)) - (31*(idx/16))] = idx;
    }

    int hybridCount[MAX_RCE][MAX_FEB][MAX_HYB];
    int *channelCount[MAX_RCE][MAX_FEB][MAX_HYB];
    int *channelAllCount[MAX_RCE][MAX_FEB][MAX_HYB];
    double *channelMean[7][MAX_RCE][MAX_FEB][MAX_HYB];
    double *channelVariance[7][MAX_RCE][MAX_FEB][MAX_HYB];
    for (int rce = 0;rce<MAX_RCE;rce++)
        for (int fpga = 0;fpga<MAX_FEB;fpga++)
            for (int hyb = 0;hyb<MAX_HYB;hyb++) {
                hybridCount[rce][fpga][hyb] = 0;
            }

    double          grChan[640];
    for (int i=0;i<640;i++)
    {
        grChan[i] = i;
    }
    DataRead        *dataRead;
    int svt_bank_num = -1;
    DevboardEvent    event;
    TriggerEvent    triggerevent;
    TriggerSample   *triggersample = new TriggerSample();
    int		samples[6];
    int            eventCount;
    int runCount;
    TString inname;
    TString outdir;
    char            name[200];
    TGraph          *graph[7];
    TMultiGraph *mg;

    while ((c = getopt(argc,argv,"ho:nmct:H:F:e:Es:b:V")) !=-1)
        switch (c)
        {
            case 'h':
                printf("-h: print this help\n");
                printf("-o: use specified output filename\n");
                printf("-n: DAQ (Ryan's) channel numbering\n");
                printf("-m: number channels in raw mux order\n");
                printf("-c: don't compute correlations\n");
                printf("-t: hybrid type (1 for old test run hybrid, 2 for new 2014 hybrid)\n");
                printf("-F: use only specified FPGA\n");
                printf("-H: use only specified hybrid\n");
                printf("-e: stop after specified number of events\n");
                printf("-E: use EVIO file format\n");
                printf("-s: number of sigmas for threshold file (default 2)\n");
                printf("-b: EVIO bank number for SVT (default 3)\n");
                printf("-V: use TriggerEvent event format\n");
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
            case 's':
                threshold_sigma = atof(optarg);
                break;
            case 'b':
                svt_bank_num = atoi(optarg);
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

    if (evio_format) {
        DataReadEvio *tmpDataRead = new DataReadEvio();
        if (triggerevent_format)
            tmpDataRead->set_engrun(true);
        if (svt_bank_num>0)
            tmpDataRead->set_bank_num(svt_bank_num);
        dataRead = tmpDataRead;
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

    // Start X11 view
    //   TApplication theApp("App",NULL,NULL);

    // Root file is the first and only arg
    if ( argc-optind != 1 ) {
        cout << "Usage: meeg_baseline data_file\n";
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
    ofstream basefile;
    cout << "Writing calibration to " << inname+".base" << endl;
    basefile.open(inname+".base");
    basefile << "#" << inname << endl;

    ofstream threshfile;
    cout << "Writing thresholds to " << inname+".thresholds" << endl;
    threshfile.open(inname+".thresholds");
    threshfile << "%" << inname << endl;

    ofstream outfile;
    cout << "Writing calibration to " << inname+".basecal" << endl;
    outfile.open(inname+".basecal");
    outfile << "#" << inname << endl;

    cout << "Reading data file " <<argv[optind] << endl;
    // Attempt to open data file
    if ( ! dataRead->open(argv[optind]) ) return(2);

    TString confname=argv[optind];
    confname.ReplaceAll(".bin","");
    confname.Append(".conf");
    if (confname.Contains('/')) {
        confname.Remove(0,confname.Last('/')+1);
    }

    ofstream outconfig;
    cout << "Writing configuration to " <<outdir<<confname << endl;
    outconfig.open(outdir+confname);

    bool readOK;

    if (triggerevent_format) {
        dataRead->next(&triggerevent);
    } else {
        dataRead->next(&event);
    }
    outconfig << dataRead->getConfigXml();
    outconfig << endl;
    outconfig << dataRead->getStatusXml();
    outconfig.close();

    runCount = atoi(dataRead->getConfig("RunCount").c_str());

    // Process each event
    eventCount = 0;

    do {
        int rce = 0;
        int fpga = 0;
        int samplecount;

        //printf("fpga %d\n",event.fpgaAddress());
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
        if (num_events!=-1 && eventCount >= num_events) break;
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
            if (use_fpga!=-1 && fpga!=use_fpga) continue;
            if (use_hybrid!=-1 && hyb!=use_hybrid) continue;
            if (!goodSample) continue;
            //printf("fpga %d, hybrid %d, apv %d, chan %d\n",fpga,hyb,apv,apvch);
            //printf("event %d\tx=%d\tR%d F%d H%d A%d channel %d, samples:\t%d\t%d\t%d\t%d\t%d\t%d\n",eventCount,x,rce,fpga,hyb,apv,apvch,samples[0],samples[1],samples[2],samples[3],samples[4],samples[5]);

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


            // Filter APVs
            if ( eventCount >= 20 ) {
                if (hybridCount[rce][fpga][hyb]==0) {
                    printf("found new hybrid: rce = %d, feb = %d, hyb = %d\n",rce,fpga,hyb);
                    channelCount[rce][fpga][hyb] = new int[640];
                    channelAllCount[rce][fpga][hyb] = new int[640];
                    for (int j=0;j<7;j++) {
                        channelMean[j][rce][fpga][hyb] = new double[640];
                        channelVariance[j][rce][fpga][hyb] = new double[640];
                    }
                    for (int i=0;i<640;i++) {
                        channelCount[rce][fpga][hyb][i] = 0;
                        channelAllCount[rce][fpga][hyb][i] = 0;
                        for (int j=0;j<7;j++) {
                            channelMean[j][rce][fpga][hyb][i] = 0.0;
                            channelVariance[j][rce][fpga][hyb][i] = 0.0;
                        }
                    }
                }
                hybridCount[rce][fpga][hyb]++;
                channelCount[rce][fpga][hyb][channel]++;
                for ( int y=0; y < 6; y++ ) {
                    //printf("%x\n",sample->value(y));
                    int value = samples[y];
                    if (value<1000)
                        printf("out of range: event %d, rce = %d, feb = %d, hyb = %d, channel = %d, sample[%d] = %d\n",eventCount,rce,fpga,hyb,channel,y,samples[y]);
                    double delta = value-channelMean[y][rce][fpga][hyb][channel];
                    if (channelCount[rce][fpga][hyb][channel]==1)
                    {
                        channelMean[y][rce][fpga][hyb][channel] = value;
                    }
                    else
                    {
                        channelMean[y][rce][fpga][hyb][channel] += delta/channelCount[rce][fpga][hyb][channel];
                    }
                    //printf("%d\n",value);
                    //printf("%f\n",channelMean[y][fpga][hyb][channel]);
                    channelVariance[y][rce][fpga][hyb][channel] += delta*(value-channelMean[y][rce][fpga][hyb][channel]);

                    channelAllCount[rce][fpga][hyb][channel]++;
                    delta = value-channelMean[6][rce][fpga][hyb][channel];
                    if (channelAllCount[rce][fpga][hyb][channel]==1)
                    {
                        channelMean[6][rce][fpga][hyb][channel] = value;
                    }
                    else
                    {
                        channelMean[6][rce][fpga][hyb][channel] += delta/channelAllCount[rce][fpga][hyb][channel];
                    }
                    channelVariance[6][rce][fpga][hyb][channel] += delta*(value-channelMean[6][rce][fpga][hyb][channel]);
                }
            }
        }
        eventCount++;
        /*
           for (int i=0;i<640;i++)
           {
           if (channelMean[6][i]<200) printf("event %d, channel %d, %f\n",eventCount,i,channelMean[6][i]);
           }*/

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

    for (int rce = 0;rce<MAX_RCE;rce++)
        for (int fpga = 0;fpga<MAX_FEB;fpga++)
            for (int hyb = 0;hyb<MAX_HYB;hyb++) if (hybridCount[rce][fpga][hyb])
            {
                printf("processing hybrid: rce = %d, feb = %d, hyb = %d\n",rce,fpga,hyb);
                int deadAPV = -1;
                for (int i=0;i<640;i++)
                {
                    if (channelCount[rce][fpga][hyb][i])
                    {
                        if (i/128==deadAPV)
                        {
                            printf("Counted %d events on channel %d, even though we thought APV %d was dead\n",channelCount[rce][fpga][hyb][i],i,deadAPV);
                        }
                        if (!evio_format && channelCount[rce][fpga][hyb][i]!=(int)eventCount-20)
                        {
                            printf("Counted %d events for channel %d; expected %d\n",channelCount[rce][fpga][hyb][i],i,eventCount-20);
                        }
                        for (int j=0;j<6;j++)
                        {
                            channelVariance[j][rce][fpga][hyb][i]/=channelCount[rce][fpga][hyb][i];
                            channelVariance[j][rce][fpga][hyb][i]=sqrt(channelVariance[j][rce][fpga][hyb][i]);
                        }
                        channelVariance[6][rce][fpga][hyb][i]/=channelAllCount[rce][fpga][hyb][i];
                        channelVariance[6][rce][fpga][hyb][i]=sqrt(channelVariance[6][rce][fpga][hyb][i]);
                    }
                    else
                    {
                        if (i/128!=deadAPV)
                        {
                            deadAPV = i/128;
                            printf("No events on FPGA %d, hybrid %d, channel %d; assuming APV %d is dead\n",fpga,hyb,i,deadAPV);
                        }
                    }
                }
            }

    //printf("correlation %f %f\n",corrHist->GetCorrelationFactor(),channelCovar[corr1][corr2]);
    //printf("correlation %f\n",channelCovar[corr1][corr2]);

    for (int rce = 0;rce<MAX_RCE;rce++)
        for (int fpga = 0;fpga<MAX_FEB;fpga++)
            for (int hyb = 0;hyb<MAX_HYB;hyb++) if (hybridCount[rce][fpga][hyb]) 
            {
                c1->Clear();
                mg = new TMultiGraph();
                for (int i=0;i<7;i++)
                {
                    graph[i] = new TGraph(640,grChan,channelMean[i][rce][fpga][hyb]);
                    graph[i]->SetMarkerColor(i+2);
                    mg->Add(graph[i]);
                }
                graph[6]->SetMarkerColor(1);
                mg->SetTitle("Pedestal;Channel;ADC counts");
                //mg->GetXaxis()->SetRangeUser(0,640);
                mg->Draw("a*");
                sprintf(name,"%s_pedestal_R%d_F%d_H%d.png",inname.Data(),rce,fpga,hyb);
                c1->SaveAs(name);
                for (int i=0;i<7;i++)
                    delete graph[i];
                delete mg;

                c1->Clear();
                mg = new TMultiGraph();
                for (int i=0;i<7;i++)
                {
                    graph[i] = new TGraph(640,grChan,channelVariance[i][rce][fpga][hyb]);
                    graph[i]->SetMarkerColor(i+2);
                    mg->Add(graph[i]);
                }
                graph[6]->SetMarkerColor(1);
                mg->SetTitle("Noise;Channel;ADC counts");
                //mg->GetXaxis()->SetRangeUser(0,640);
                mg->Draw("a*");
                double max = 0;
                for (int i=0;i<7;i++)
                    for (int j=0;j<640;j++)
                        if (channelVariance[i][rce][fpga][hyb][j]>max) max = channelVariance[i][rce][fpga][hyb][j];
                mg->GetYaxis()->SetRangeUser(0,1.2*max);
                sprintf(name,"%s_noise_R%d_F%d_H%d.png",inname.Data(),rce,fpga,hyb);
                c1->SaveAs(name);
                for (int i=0;i<7;i++)
                    delete graph[i];
                delete mg;

                c1->Clear();
                mg = new TMultiGraph();
                double grShift[640];
                for (int i=0;i<6;i++)
                {
                    for (int n=0;n<640;n++) grShift[n]=channelMean[i][rce][fpga][hyb][n]-channelMean[6][rce][fpga][hyb][n];
                    graph[i] = new TGraph(640,grChan,grShift);
                    graph[i]->SetMarkerColor(i+1);
                    mg->Add(graph[i]);
                }
                mg->SetTitle("Sample-to-sample shift;Channel;ADC counts");
                //mg->GetXaxis()->SetRangeUser(0,640);
                mg->Draw("a*");
                sprintf(name,"%s_shift_R%d_F%d_H%d.png",inname.Data(),rce,fpga,hyb);
                c1->SaveAs(name);
                for (int i=0;i<6;i++)
                    delete graph[i];
                delete mg;

                for (int i=0;i<640;i++)
                {
                    int apv = i/128;
                    int channel = i%128;
                    if (flip_channels) apv = 4-apv;
                    if (channel == 0)
                        threshfile << fpga << "," << hyb << "," << apv << endl;
                    threshfile << apv*128+channel << "," << channelMean[6][rce][fpga][hyb][i] + threshold_sigma*channelVariance[6][rce][fpga][hyb][i] << endl;
                }
                for (int i=0;i<640;i++) if (channelCount[rce][fpga][hyb][i]>0)
                {
                    int apv = i/128;
                    int channel = i%128;
                    if (flip_channels) apv = 4-apv;
                    int sensorChannel =  apv*128+channel;
                    outfile << fpga << "\t" << hyb << "\t" << sensorChannel << "\t" << channelMean[6][rce][fpga][hyb][i] << "\t" << channelVariance[6][rce][fpga][hyb][i] << endl;

                    basefile << rce << "\t" << fpga << "\t" << hyb << "\t" << sensorChannel << "\t";
                    for (int j=0;j<7;j++)
                    {
                        basefile<<channelMean[j][rce][fpga][hyb][i]<<"\t"<<channelVariance[j][rce][fpga][hyb][i]<<"\t";
                    }
                    basefile<<endl;
                }
            }
    // Start X-Windows
    //theApp.Run();

    // Close file
    threshfile.close();
    outfile.close();
    basefile.close();
    delete dataRead;
    return(0);
}

