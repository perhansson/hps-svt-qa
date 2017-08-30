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
#include <TF1.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <stdarg.h>
#include <DevboardEvent.h>
#include <DevboardSample.h>
#include <Data.h>
#include <DataRead.h>
#include <DataReadEvio.h>
#include "ShapingCurve.hh"
#include "SmoothShapingCurve.hh"
#include "Samples.hh"
#include "Fitter.hh"
#include "AnalyticFitter.hh"
#include "LinFitter.hh"
#include <TMath.h>
#include "meeg_utils.hh"
#include <unistd.h>
#include <TMVA/TSpline1.h>
using namespace std;

// Process the data
// Pass root file to open as first and only arg.
int main ( int argc, char **argv ) {
	int c;
	bool subtract_T0 = false;
	bool print_fit_status = false;
	bool force_cal_grp = false;
	bool ignore_cal_grp = false;
	bool flip_channels = true;
	bool use_shape = false;
	bool shift_t0 = false;
	bool use_dist = false;
	bool evio_format = false;
	int fpga = -1;
	int hybrid = -1;
	int num_events = -1;
	double dist_window;
	bool make_shape = false;
	int single_channel = -1;
	TString inname = "";
	TString outdir = "";
	int cal_grp = -1;
	int cal_delay = 0;
	double delay_step = SAMPLE_INTERVAL/8;
	TCanvas         *c1;
	//TH2F            *histAll;
	double          histMin[640];
	double          histMax[640];
	//TGraph          *mean;
	//TGraph          *sigma;
	int nChan[2] = {0, 0};
	double          grChan[2][640];
	double grT0[2][640], grT0_sigma[2][640], grT0_err[2][640];
	double grA[2][640], grA_sigma[2][640], grA_err[2][640];
	double          calMean[640][7] = {{0.0}};
	double          calSigma[640][7] = {{1.0}};
	double calTp[2][640] = {{0.0}};
	double calA[2][640] = {{0.0}};
	double calT0[2][640] = {{0.0}};
	double calChisq[2][640] = {{0.0}};
	DataRead        *dataRead;
	DevboardEvent    event;
	DevboardSample   *sample;
	uint            x;
	uint            y;
	uint            value;
	int            channel;
	int            eventCount;
	int runCount;
	double          sum;
	char            filename[100];
	char            name[100];
	char            name2[100];
	char title[200];
	ShapingCurve *myShape[2][640];
	Fitter *myFitter[2][640];

	TH1F *histT0[2][640];
	TH1F *histA[2][640];
	TH1F *histT0_err[2][640];
	TH1F *histA_err[2][640];
	TH2F *histChiProb[2];
	TH2F *histT0_2d[2];
	TH2F *histA_2d[2];
	TH2I *pulse2D[2];
	TH2I *T0_A[2];
	double maxA[2] = {0, 0};
	double minA[2] = {16384, 16384};
	double maxT0[2] = {0.0, 0.0};
	double minT0[2] = {200.0, 200.0};
	TGraph *T0_dist[2];
	double T0_dist_b[2],T0_dist_m[2];


	while ((c = getopt(argc,argv,"hfsg:o:auc:d:tbnH:F:e:E")) !=-1)
		switch (c)
		{
			case 'h':
				printf("-h: print this help\n");
				printf("-f: print fit status to .fits\n");
				printf("-g: force use of specified cal group\n");
				printf("-c: use only specified channel\n");
				printf("-a: use all cal groups\n");
				printf("-o: use specified output filename\n");
				printf("-s: shift T0 by value from Tp cal file\n");
				printf("-u: use .shape cal instead of .tp\n");
				printf("-b: subtract value of T0 in .tp cal\n");
				printf("-t: make new .shape cal based on T0-A distribution\n");
				printf("-d: read and use .dist file, starting T0 window at specified value\n");
				printf("-n: DAQ (Ryan's) channel numbering\n");
				printf("-F: use only specified FPGA\n");
				printf("-H: use only specified hybrid\n");
				printf("-e: stop after specified number of events\n");
				printf("-E: use EVIO file format\n");
				return(0);
				break;
			case 'f':
				print_fit_status = true;
				break;
			case 'u':
				use_shape = true;
				break;
			case 'n':
				flip_channels = false;
				break;
			case 'a':
				ignore_cal_grp = true;
				break;
			case 's':
				shift_t0 = true;
				break;
			case 'b':
				subtract_T0 = true;
				break;
			case 'c':
				single_channel = atoi(optarg);
				break;
			case 'g':
				cal_grp = atoi(optarg);
				force_cal_grp=true;
				break;
			case 'o':
				inname = optarg;
				outdir = optarg;
				if (outdir.Contains('/')) {
					outdir.Remove(outdir.Last('/')+1);
				}
				else outdir="";
				break;
			case 'd':
				dist_window = atof(optarg);
				use_dist=true;
				break;
			case 't':
				make_shape = true;
				break;
			case 'F':
				fpga = atoi(optarg);
				break;
			case 'H':
				hybrid = atoi(optarg);
				break;
			case 'e':
				num_events = atoi(optarg);
				break;
			case 'E':
				evio_format = true;
				break;
			case '?':
				printf("Invalid option or missing option argument; -h to list options\n");
				return(1);
			default:
				abort();
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
	c1 = new TCanvas("c1","c1",1200,900);
	//c1->SetLogz();
	// Start X11 view
	//TApplication theApp("App",NULL,NULL);

	// Root file is the first and only arg
	if ( argc-optind<3 ) {
		cout << "Usage: meeg_t0res baseline_cal tp_cal data_file\n";
		return(1);
	}

	ifstream calfile;
	cout << "Reading baseline calibration from " << argv[optind] << endl;
	calfile.open(argv[optind]);
	while (!calfile.eof()) {
		calfile >> channel;
		if (calfile.eof()) break;
		for (int i=0;i<7;i++)
		{
			calfile >> calMean[channel][i];
			calfile >> calSigma[channel][i];
		}
	}
	calfile.close();

	optind++;

	for (int sgn=0;sgn<2;sgn++)
	{
		sprintf(name,"%s.tp_%s",argv[optind],sgn?"neg":"pos");
		cout << "Reading Tp calibration from "<<name<<endl;
		calfile.open(name);
		while (!calfile.eof()) {
			calfile >> channel;
			if (calfile.eof()) break;
			calfile >> calA[sgn][channel];
			calfile >> calT0[sgn][channel];
			calfile >> calTp[sgn][channel];
			calfile >> calChisq[sgn][channel];

		}
		calfile.close();
	}

	for (int sgn=0;sgn<2;sgn++)
	{
		if (use_shape)
		{
			int ni;
			double ti[300], yi[300], ey[300];
			double shape_t0, shape_y0, shape_e0;
			sprintf(name,"%s.shape_%s",argv[optind],sgn?"neg":"pos");
			cout << "Reading shape calibration from "<<name<<endl;
			calfile.open(name);
			while (!calfile.eof()) {
				calfile >> channel;
				calfile >> ni;
				calfile >> shape_t0;
				calfile >> shape_y0;
				calfile >> shape_e0;
				if (calfile.eof()) break;
				int j=0;
				do
				{
					ti[j] = -100.0+5.0*j;
					yi[j] = 0.0;
					ey[j] = 1.0;
					j++;
				} while (ti[j]<shape_t0);
				ti[j] = shape_t0;
				yi[j] = shape_y0;
				ey[j] = shape_e0;
				for (int i=1;i<ni;i++)
				{
					calfile >> ti[i+j];
					calfile >> yi[i+j];
					calfile >> ey[i+j];
					//			if (ti[i]<0) yi[i]=0;
				}
				ni+=j;
				do
				{
					ti[ni] = ti[ni-1]+5.0;
					yi[ni] = yi[ni-1];
					ey[ni] = 1.0;
					ni++;
				} while (ti[ni-1]<300.0);

				myShape[sgn][channel] = new SmoothShapingCurve(ni,ti,yi);
				myFitter[sgn][channel] = new LinFitter(myShape[sgn][channel],6,1,TMath::Mean(6,calSigma[channel]));


				if (single_channel!=-1 && channel==single_channel)
					//if (channel==639)
				{
					c1->Clear();
					TSpline3 *tempspline = new TSpline3("tempspline",ti,yi,ni,"",0.0,yi[ni-1]);
					tempspline->Draw();
					tempspline->Print();
					sprintf(name,"spline_%s.png",sgn?"neg":"pos");
					c1->SaveAs(name);
				}
			}
			calfile.close();
		}
		else
		{
			for (channel = 0; channel<640; channel++)
			{
				//myShape[sgn][channel] = new SmoothShapingCurve(calTp[sgn][channel]);
				//myFitter[sgn][channel] = new LinFitter(myShape[sgn][channel],6,1,calSigma[channel]);
				myShape[sgn][channel] = new ShapingCurve(calTp[sgn][channel]);
				myFitter[sgn][channel] = new AnalyticFitter(myShape[sgn][channel],6,1,TMath::Mean(6,calSigma[channel]));
			}
		}

		if (use_dist)
		{
			int ni = 0;
			double ti[500], integral[500], amplitude[500];
			sprintf(name,"%s.dist_%s",argv[optind],sgn?"neg":"pos");
			cout << "Reading T0-A calibration from "<<name<<endl;
			/*
			   ti[0] = -100.0;
			   integral[0] = -100.0;
			   amplitude[0] = 0.0;
			   ni++;
			   */
			calfile.open(name);
			while (!calfile.eof()) {
				calfile >> ti[ni];
				calfile >> amplitude[ni];
				calfile >> integral[ni];
				if (calfile.eof()) break;
				ni++;
			}
			calfile.close();
			/*
			   ti[ni] = ti[ni-1]+100.0;
			   integral[ni] = integral[ni-1]*2;
			   amplitude[ni] = 0.0;
			   ni++;
			   */
			//sprintf(name,"T0_dist_%s",sgn?"neg":"pos");
			//T0_dist[sgn] = new TSpline3(name,ti,integral,ni);
			T0_dist[sgn] = new TGraph(ni,ti,integral);
			T0_dist_m[sgn] = SAMPLE_INTERVAL/(T0_dist[sgn]->Eval(dist_window+SAMPLE_INTERVAL) - T0_dist[sgn]->Eval(dist_window));
			T0_dist_b[sgn] = dist_window - T0_dist[sgn]->Eval(dist_window)*T0_dist_m[sgn];

			c1->Clear();
			T0_dist[sgn]->Draw("al");
			sprintf(name,"T0_dist_%s.png",sgn?"neg":"pos");
			c1->SaveAs(name);
		}
	}
	optind++;


	// 2d histogram
	for (int sgn=0;sgn<2;sgn++)
	{
		sprintf(name,"Pulse_shape_%s",sgn?"Neg":"Pos");
		sprintf(title,"Normalized pulse shape, %s pulses;Time [ns];Amplitude [normalized]",sgn?"negative":"positive");
		pulse2D[sgn] = new TH2I(name,title,520,-1.5*SAMPLE_INTERVAL,5*SAMPLE_INTERVAL,500,-0.1,1.2);
		sprintf(name,"A_vs_T0_%s",sgn?"Neg":"Pos");
		sprintf(title,"Amplitude vs. T0, %s pulses;Time [ns];Amplitude [ADC counts]",sgn?"negative":"positive");
		T0_A[sgn] = new TH2I(name,title,500,-4*SAMPLE_INTERVAL,5*SAMPLE_INTERVAL,500,0.0,2000.0);
		//T0_A[sgn] = new TH2I(name,name,500,0,2*SAMPLE_INTERVAL,500,800.0,1400.0);
		if (force_cal_grp)
		{
			sprintf(name,"Chisq_Prob_%s",sgn?"Neg":"Pos");
			sprintf(title,"Chisq probability of fit, %s pulses;Channel;Probability",sgn?"negative":"positive");
			histChiProb[sgn] = new TH2F(name,title,80,0,640,100,0,1.0);
			sprintf(name,"T0_%s",sgn?"Neg":"Pos");
			sprintf(title,"Fitted values of T0 relative to channel average, %s pulses;Channel;T0 [ns]",sgn?"negative":"positive");
			histT0_2d[sgn] = new TH2F(name,title,80,0,640,1000,-1*SAMPLE_INTERVAL,6*SAMPLE_INTERVAL);
			sprintf(name,"A_%s",sgn?"Neg":"Pos");
			sprintf(title,"Fitted values of amplitude, %s pulses;Channel;Amplitude [ADC counts]",sgn?"negative":"positive");
			histA_2d[sgn] = new TH2F(name,title,80,0,640,1000,0,2000.0);
		}
		else
		{
			sprintf(name,"Chisq_Prob_%s",sgn?"Neg":"Pos");
			sprintf(title,"Chisq probability of fit, %s pulses;Channel;Probability",sgn?"negative":"positive");
			histChiProb[sgn] = new TH2F(name,title,640,0,640,100,0,1.0);
			sprintf(name,"T0_%s",sgn?"Neg":"Pos");
			sprintf(title,"Fitted values of T0 relative to channel average, %s pulses;Channel;T0 [ns]",sgn?"negative":"positive");
			histT0_2d[sgn] = new TH2F(name,title,640,0,640,1000,-1*SAMPLE_INTERVAL,6*SAMPLE_INTERVAL);
			sprintf(name,"A_%s",sgn?"Neg":"Pos");
			sprintf(title,"Fitted values of amplitude, %s pulses;Channel;Amplitude [ADC counts]",sgn?"negative":"positive");
			histA_2d[sgn] = new TH2F(name,title,640,0,640,1000,0,2000.0);
		}
	}

	for (int n=0;n<640;n++) {
		histMin[n] = 16384;
		histMax[n] = 0;
		for (int sgn=0;sgn<2;sgn++)
		{
			sprintf(name,"T0_%s_%i",sgn?"neg":"pos",n);
			histT0[sgn][n] = new TH1F(name,name,1000,0,3*SAMPLE_INTERVAL);
			sprintf(name,"T0err_%s_%i",sgn?"neg":"pos",n);
			histT0_err[sgn][n] = new TH1F(name,name,1000,0,10.0);
			sprintf(name,"A_%s_%i",sgn?"neg":"pos",n);
			histA[sgn][n] = new TH1F(name,name,1000,0,2000.0);
			sprintf(name,"Aerr_%s_%i",sgn?"neg":"pos",n);
			histA_err[sgn][n] = new TH1F(name,name,1000,0,100.0);
		}
	}


	if (inname=="")
	{
		inname=argv[optind];

		inname.ReplaceAll(".bin","");
		if (inname.Contains('/')) {
			inname.Remove(0,inname.Last('/')+1);
		}
	}



	ofstream fitfile;
	if (print_fit_status) {
		cout << "Writing fit status to " << inname+".fits" << endl;
		fitfile.open(inname+".fits");
	}

	ofstream outfile[2];
	cout << "Writing T0 calibration to " << inname+".t0_pos" << endl;
	outfile[0].open(inname+".t0_pos");
	cout << "Writing T0 calibration to " << inname+".t0_neg" << endl;
	outfile[1].open(inname+".t0_neg");

	ofstream outshapefile[2];
	if (make_shape)
	{
		cout << "Writing pulse shape to " << inname+".shape_pos" << endl;
		outshapefile[0].open(inname+".shape_pos");
		cout << "Writing pulse shape to " << inname+".shape_neg" << endl;
		outshapefile[1].open(inname+".shape_neg");
	}

	ofstream outdistfile[2];
	cout << "Writing T0 and A distribution to " << inname+".dist_pos" << endl;
	outdistfile[0].open(inname+".dist_pos");
	cout << "Writing T0 and A distribution to " << inname+".dist_neg" << endl;
	outdistfile[1].open(inname+".dist_neg");

	double samples[6];
	Samples *mySamples = new Samples(6,SAMPLE_INTERVAL);
	double fit_par[2], fit_err[2], chisq, chiprob;
	int dof;


	while (optind<argc)
	{
		cout << "Reading data file " <<argv[optind] << endl;
		// Attempt to open data file
		if ( ! dataRead->open(argv[optind]) ) return(2);

		TString confname=argv[optind];
		confname.ReplaceAll(".bin",".conf");
		if (confname.Contains('/')) {
			confname.Remove(0,confname.Last('/')+1);
		}

		ofstream outconfig;
		cout << "Writing configuration to " <<outdir<<confname << endl;
		outconfig.open(outdir+confname);

		dataRead->next(&event);
		dataRead->dumpConfig(outconfig);
		outconfig.close();
		//dataRead.dumpStatus();
		
		runCount = atoi(dataRead->getConfig("RunCount").c_str());
		if (!force_cal_grp)
		{
			cal_grp = atoi(dataRead->getConfig("cntrlFpga:hybrid:apv25:CalGroup").c_str());
			cout<<"Read calibration group "<<cal_grp<<" from data file"<<endl;
		}

		cal_delay = atoi(dataRead->getConfig("cntrlFpga:hybrid:apv25:Csel").substr(4,1).c_str());
		cout<<"Read calibration delay "<<cal_delay<<" from data file"<<endl;
		if (cal_delay==0)
		{
			cal_delay=8;
			cout<<"Force cal_delay=8 to keep sample time in range"<<endl;
		}


		// Process each event
		eventCount = 0;

		do {
			if (fpga!=-1 && event.fpgaAddress()!=fpga) continue;
			if (eventCount%1000==0) printf("Event %d\n",eventCount);
			for (x=0; x < event.count(); x++) {
				// Get sample
				sample  = event.sample(x);
				if (hybrid!=-1 && sample->hybrid()!=hybrid) continue;

				channel = sample->channel();
				if (flip_channels)
					channel += (4-sample->apv())*128;
				else
					channel += sample->apv()*128;
				//				if (sample->apv()==0 || sample->apv()==4) continue;

				if (single_channel!=-1 && channel !=single_channel) continue;
				if ( channel >= (5 * 128) ) {
					cout << "Channel " << dec << channel << " out of range" << endl;
					cout << "Apv = " << dec << sample->apv() << endl;
					cout << "Chan = " << dec << sample->channel() << endl;
				}

				if (!ignore_cal_grp && cal_grp!=-1 && ((int)sample->channel()-cal_grp)%8!=0) continue;
				int n = channel;
				// Filter APVs
				if ( eventCount >= 20 ) {

					int samplesAbove = 0;
					int samplesBelow = 0;
					sum = 0;
					for ( y=0; y < 6; y++ ) {
						value = sample->value(y);

						//vhigh = (value << 1) & 0x2AAA;
						//vlow  = (value >> 1) & 0x1555;
						//value = vlow | vhigh;


						if ( value < histMin[n] ) histMin[n] = value;
						if ( value > histMax[n] ) histMax[n] = value;
						samples[y] = value;
						samples[y] -= calMean[channel][y];
						//samples[y] -= sample->value(0);
						sum+=samples[y];
						if (samples[y]>5*calSigma[channel][y]) samplesAbove++;
						if (samples[y]<-5*calSigma[channel][y]) samplesBelow++;
					}
					if (sum<0)
						for ( y=0; y < 6; y++ ) {
							samples[y]*=-1;
						}
					if (samplesAbove>1 || samplesBelow>1)
					{
						int sgn = sum>0?0:1;
						mySamples->readEvent(samples,0.0);
						myFitter[sgn][n]->readSamples(mySamples);
						myFitter[sgn][n]->doFit();
						myFitter[sgn][n]->getFitPar(fit_par);
						myFitter[sgn][n]->getFitErr(fit_err);
						chisq = myFitter[sgn][n]->getChisq(fit_par);
						dof = myFitter[sgn][n]->getDOF();
						//if (fit_par[1]!=fit_par[1]) 
						//{
						//	printf("%d\t%d\t%f\t%f meeg\n",channel,sgn,fit_par[0],fit_par[1]);
						//	mySamples->print();
						//	myFitter[sgn][n]->print_fit();
						//	myFitter[sgn][n]->printResiduals();
						//}

						if (use_dist)
						{
							if (fit_par[0]>dist_window && fit_par[0]<dist_window+SAMPLE_INTERVAL)
							{
								fit_par[0] = T0_dist_m[sgn]*T0_dist[sgn]->Eval(fit_par[0]) + T0_dist_b[sgn];
							}
							else continue;
						}

						if (subtract_T0) fit_par[0]-=calT0[sgn][channel];


						histT0[sgn][n]->Fill(fit_par[0]);
						histT0_err[sgn][n]->Fill(fit_err[0]);
						histA[sgn][n]->Fill(fit_par[1]);
						histA_err[sgn][n]->Fill(fit_err[1]);
						T0_A[sgn]->Fill(fit_par[0],fit_par[1]);
						chiprob = TMath::Prob(chisq,dof);
						if (print_fit_status) fitfile<<"Channel "<<channel << ", T0 " << fit_par[0] <<", A " << fit_par[1] << ", Fit chisq " << chisq << ", DOF " << dof << ", prob " << chiprob << endl;
						if (fit_par[0]>maxT0[sgn]) maxT0[sgn] = fit_par[0];
						if (fit_par[0]<minT0[sgn]) minT0[sgn] = fit_par[0];
						if (fit_par[1]>maxA[sgn]) maxA[sgn] = fit_par[1];
						if (fit_par[1]<minA[sgn]) minA[sgn] = fit_par[1];
						histT0_2d[sgn]->Fill(channel,fit_par[0]);
						histA_2d[sgn]->Fill(channel,fit_par[1]);
						histChiProb[sgn]->Fill(channel,chiprob);
						for ( y=0; y < 6; y++ ) {
							pulse2D[sgn]->Fill(y*SAMPLE_INTERVAL-fit_par[0],samples[y]/fit_par[1]);
						}
					}
					else
					{
						//					fitfile << "Pulse below threshold on channel "<<channel<<endl;
					}
				}
			}
			eventCount++;

		} while ( dataRead->next(&event) );
		dataRead->close();
		if (eventCount != runCount)
		{
			printf("ERROR: events read = %d, runCount = %d\n",eventCount, runCount);
		}
		optind++;
	}

	for (int n=0;n<640;n++) for (int sgn=0;sgn<2;sgn++) if (histT0[sgn][n]->GetEntries()>0) {
		grChan[sgn][nChan[sgn]]=n;

		/*
		   histT0[sgn][n]->Fit("gaus","Q0");
		   grT0[sgn][nChan[sgn]] = histT0[sgn][n]->GetFunction("gaus")->GetParameter(1);
		   grT0_sigma[sgn][nChan[sgn]] = histT0[sgn][n]->GetFunction("gaus")->GetParameter(2);
		   histT0_err[sgn][n]->Fit("gaus","Q0");
		   grT0_err[sgn][nChan[sgn]] = histT0_err[sgn][n]->GetFunction("gaus")->GetParameter(1);
		   histA[sgn][n]->Fit("gaus","Q0");
		   grA[sgn][nChan[sgn]] = histA[sgn][n]->GetFunction("gaus")->GetParameter(1);
		   grA_sigma[sgn][nChan[sgn]] = histA[sgn][n]->GetFunction("gaus")->GetParameter(2);
		   histA_err[sgn][n]->Fit("gaus","Q0");
		   grA_err[sgn][nChan[sgn]] = histA_err[sgn][n]->GetFunction("gaus")->GetParameter(1);
		   */
		grT0[sgn][nChan[sgn]] = histT0[sgn][n]->GetMean();
		if (shift_t0) grT0[sgn][nChan[sgn]] -= cal_delay*delay_step + calT0[sgn][n];
		grT0_sigma[sgn][nChan[sgn]] = histT0[sgn][n]->GetRMS();
		grT0_err[sgn][nChan[sgn]] = histT0_err[sgn][n]->GetMean();
		grA[sgn][nChan[sgn]] = histA[sgn][n]->GetMean();
		grA_sigma[sgn][nChan[sgn]] = histA[sgn][n]->GetRMS();
		grA_err[sgn][nChan[sgn]] = histA_err[sgn][n]->GetMean();

		outfile[sgn] <<n<<"\t"<<grT0[sgn][nChan[sgn]]<<"\t\t"<<grT0_sigma[sgn][nChan[sgn]]<<"\t\t"<<grT0_err[sgn][nChan[sgn]]<<"\t\t";
		outfile[sgn]<<grA[sgn][nChan[sgn]]<<"\t\t"<<grA_sigma[sgn][nChan[sgn]]<<"\t\t"<<grA_err[sgn][nChan[sgn]]<<endl;     
		nChan[sgn]++;
	}

	for (int sgn=0;sgn<2;sgn++)
	{
		pulse2D[sgn]->Draw("colz");
		//sprintf(name,"Pulse_profile_%s_1",sgn?"Neg":"Pos");
		TH2I *tempPulse;
		int *tempArray;
		double *pulse_yi, *pulse_ti, *pulse_ey;
		TGraphErrors *graphErrors;
		if (make_shape)
		{
			tempPulse = (TH2I*) pulse2D[sgn]->Clone("temp_pulse");
			tempArray = new int[tempPulse->GetNbinsY()];
			pulse_yi = new double[tempPulse->GetNbinsX()];
			pulse_ti = new double[tempPulse->GetNbinsX()];
			pulse_ey = new double[tempPulse->GetNbinsX()];
			int pulse_ni = 0;

			double center, spread;
			int count;
			tempPulse->Rebin2D(10,5);
			for (int i=0;i<tempPulse->GetNbinsX();i++)
			{
				for (int j=0;j<tempPulse->GetNbinsY();j++) tempArray[j] = (int) tempPulse->GetBinContent(i+1,j+1);

				doStats(tempPulse->GetNbinsY(),0,tempPulse->GetNbinsY()-1,tempArray,count,center,spread);
				if (count)
				{
					pulse_yi[pulse_ni] = tempPulse->GetYaxis()->GetBinCenter(1)+center*tempPulse->GetYaxis()->GetBinWidth(1);
					pulse_ey[pulse_ni] = spread*tempPulse->GetYaxis()->GetBinWidth(1)/sqrt(count);
					pulse_ti[pulse_ni] = tempPulse->GetXaxis()->GetBinCenter(1)+i*tempPulse->GetXaxis()->GetBinWidth(1);
					pulse_ni++;
				}
			}
			for (channel=0;channel<640;channel++)
			{
				outshapefile[sgn] << channel << "\t";
				outshapefile[sgn] << pulse_ni << "\t";
				for (int i=0; i<pulse_ni; i++)
				{
					outshapefile[sgn] << pulse_ti[i] << "\t";
					outshapefile[sgn] << pulse_yi[i] << "\t";
					outshapefile[sgn] << pulse_ey[i] << "\t";
				}
				outshapefile[sgn] << endl;
			}
			outshapefile[sgn].close();
			graphErrors = new TGraphErrors(pulse_ni,pulse_ti,pulse_yi,0,pulse_ey);
			graphErrors->Draw();
		}
		sprintf(name,"%s_t0_pulse_%s.png",inname.Data(),sgn?"neg":"pos");
		c1->SaveAs(name);
		if (make_shape)
		{
			delete pulse_yi;
			delete pulse_ti;
			delete pulse_ey;
			delete tempPulse;
			delete tempArray;
			delete graphErrors;
		}
		/*
		   pulseProfile[sgn] = pulse2D[sgn]->ProfileX(name);
		   pulseProfile[sgn]->Rebin(10);
		   pulseProfile[sgn]->SetLineWidth(2.0);
		   pulseProfile[sgn]->Draw("SAME");
		   sprintf(name,"%s_t0_pulse_%s.png",inname.Data(),sgn?"neg":"pos");
		   c1->SaveAs(name);
		   for (channel=0;channel<640;channel++)
		   {
		   outshapefile[sgn] << channel << "\t";
		   outshapefile[sgn] << pulseProfile[sgn]->GetNbinsX() << "\t";
		   for (int i=0; i<pulseProfile[sgn]->GetNbinsX(); i++)
		   {
		   outshapefile[sgn] << pulseProfile[sgn]->GetBinCenter(i+1) << "\t";
		   outshapefile[sgn] << pulseProfile[sgn]->GetBinContent(i+1) << "\t";
		   outshapefile[sgn] << pulseProfile[sgn]->GetBinError(i+1) << "\t";
		   }
		   outshapefile[sgn] << endl;
		   }
		   outshapefile[sgn].close();
		   */

		T0_A[sgn]->GetYaxis()->SetRangeUser(minA[sgn],maxA[sgn]);
		T0_A[sgn]->GetXaxis()->SetRangeUser(minT0[sgn]-5.0,maxT0[sgn]+5.0);
		T0_A[sgn]->Draw("colz");
		sprintf(name,"%s_t0_T0_A_%s.png",inname.Data(),sgn?"neg":"pos");
		c1->SaveAs(name);

		tempArray = new int[T0_A[sgn]->GetNbinsY()];
		int count_integral = 0;
		for (int i=0;i<T0_A[sgn]->GetNbinsX();i++)
		{
			for (int j=0;j<T0_A[sgn]->GetNbinsY();j++) tempArray[j] = (int) T0_A[sgn]->GetBinContent(i+1,j+1);

			double center, spread;
			int count;
			doStats(T0_A[sgn]->GetNbinsY(),0,T0_A[sgn]->GetNbinsY()-1,tempArray,count,center,spread);
			count_integral+=count;
			if (count>0)
			{
				outdistfile[sgn] << T0_A[sgn]->GetXaxis()->GetBinCenter(1)+i*T0_A[sgn]->GetXaxis()->GetBinWidth(1) << "\t";
				outdistfile[sgn] << T0_A[sgn]->GetYaxis()->GetBinCenter(1)+center*T0_A[sgn]->GetYaxis()->GetBinWidth(1) << "\t";
				outdistfile[sgn] << count_integral << "\t";
				outdistfile[sgn] << endl;
			}
		}
		outdistfile[sgn].close();





		histChiProb[sgn]->Draw("colz");
		sprintf(name,"%s_t0_chiprob_%s.png",inname.Data(),sgn?"neg":"pos");
		c1->SaveAs(name);

		histT0_2d[sgn]->GetYaxis()->SetRangeUser(minT0[sgn]-5.0,maxT0[sgn]+5.0);
		histT0_2d[sgn]->Draw("colz");
		sprintf(name,"%s_t0_T0_hist_%s.png",inname.Data(),sgn?"neg":"pos");
		c1->SaveAs(name);

		histA_2d[sgn]->GetYaxis()->SetRangeUser(minA[sgn],maxA[sgn]);
		histA_2d[sgn]->Draw("colz");
		sprintf(name,"%s_t0_A_hist_%s.png",inname.Data(),sgn?"neg":"pos");
		c1->SaveAs(name);

		sprintf(name,"T0_graph_%s",sgn?"neg":"pos");
		sprintf(name2,"%s_t0_T0_%s.png",inname.Data(),sgn?"neg":"pos");
		sprintf(title,"Mean fitted T0, %s pulses;Channel;T0 [ns]",sgn?"negative":"positive");
		plotResults(title, name, name2, nChan[sgn], grChan[sgn], grT0[sgn], c1);

		sprintf(name,"T0_err_%s",sgn?"neg":"pos");
		sprintf(name2,"T0_spread_%s",sgn?"neg":"pos");
		sprintf(filename,"%s_t0_T0_sigma_%s.png",inname.Data(),sgn?"neg":"pos");
		sprintf(title,"Error and spread in fitted T0, %s pulses;Channel;T0 error [ns]",sgn?"negative":"positive");
		plotResults2(title, name, name2, filename, nChan[sgn], grChan[sgn], grT0_err[sgn], grT0_sigma[sgn], c1);

		sprintf(name,"A_graph_%s",sgn?"neg":"pos");
		sprintf(name2,"%s_t0_A_%s.png",inname.Data(),sgn?"neg":"pos");
		sprintf(title,"Mean fitted amplitude, %s pulses;Channel;Amplitude [ADC counts]",sgn?"negative":"positive");
		plotResults(title, name, name2, nChan[sgn], grChan[sgn], grA[sgn], c1);

		sprintf(name,"A_err_%s",sgn?"neg":"pos");
		sprintf(name2,"A_spread_%s",sgn?"neg":"pos");
		sprintf(filename,"%s_t0_A_sigma_%s.png",inname.Data(),sgn?"neg":"pos");
		sprintf(title,"Error and spread in fitted amplitude, %s pulses;Channel;Amplitude error [ADC counts]",sgn?"negative":"positive");
		plotResults2(title, name, name2, filename, nChan[sgn], grChan[sgn], grA_err[sgn], grA_sigma[sgn], c1);
	}

	// Start X-Windows
	//	theApp.Run();

	// Close file
	if (print_fit_status) fitfile.close();
	outfile[0].close();
	outfile[1].close();
	return(0);
}

