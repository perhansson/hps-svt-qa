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
#include <Data.h>
#include <DataRead.h>
#include <DataReadEvio.h>
#include <unistd.h>
#include <Filter.h>

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
using namespace std;

//#define corr1 599
//#define corr1 604
//#define corr2 500

// Process the data
// Pass root file to open as first and only arg.
int main ( int argc, char **argv ) {
	bool read_temp = true;
	int hybrid_type = 0;
	bool evio_format = false;
	bool no_gui = false;
	bool print_data = false;
	bool use_filter = false;
	int use_fpga = -1;
	int use_hybrid = -1;
	int use_apv = -1;
	int num_events = -1;
	int c;

	int n_coeffs = 10;
	int n_delay = 1;

	Filter filt_;
	int coefOrder[] = {0,5,1,6,2,7,3,8,4,9};

	DataRead        *dataRead;
	DevboardEvent    event;
	DevboardSample   *sample;
	uint            x;
	uint            y;
	int            eventCount;
	int runCount;
	TString inname;
	TString outdir;
	char title[200];
	char name[200];


	TH1S            *hist[7][3][5][35];
	double          histMin[7][3][5][35];
	double          histMax[7][3][5][35];
	double          pedestal[7][3][5];
	int          pedestalCount[7][3][5];
	double          syncSize[7][3][5];
	int          syncCount[7][3][5];
	double          meanVal[35];
	double          plotX[35];
	ifstream        inFile;
	string          inLine;
	stringstream    inLineStream;
	string          inValue;
	stringstream    inValueStream;

	while ((c = getopt(argc,argv,"ho:nt:H:F:A:e:Egdf:")) !=-1)
		switch (c)
		{
			case 'h':
				printf("-h: print this help\n");
				printf("-o: use specified output filename\n");
				printf("-n: physical channel numbering\n");
				printf("-t: hybrid type (1 for old test run hybrid, 2 for new 2014 hybrid)\n");
				printf("-F: use only specified FPGA\n");
				printf("-H: use only specified hybrid\n");
				printf("-A: use only specified APV\n");
				printf("-e: stop after specified number of events\n");
				printf("-E: use EVIO file format\n");
				printf("-g: no GUI\n");
				printf("-d: no printing of fits\n");
				printf("-f: read and apply filter\n");
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
			case 'g':
				no_gui = true;
				break;
			case 'd':
				print_data = true;
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
			case 'A':
				use_apv = atoi(optarg);
				break;
			case 'e':
				num_events = atoi(optarg);
				break;
			case 'E':
				evio_format = true;
				break;
			case 'f':
				use_filter = true;
				int fpga,hyb,apv,chan;
				for(fpga=0; fpga < 7; fpga++) {
					for(hyb=0; hyb <3; hyb++) {
						for(apv=0; apv <5; apv++) {
							filt_.filterData[fpga][hyb][apv][0] = 1;
							for(chan=1; chan <10; chan++) {
								filt_.filterData[fpga][hyb][apv][chan] = 0;
							}
						}
					}
				}

				{
					ifstream is;
					string line;
					is.open(optarg);

					//read filter ID
					while (getline(is,line))
					{
						if (line[0]!='#') break;
					}
					line.copy(filt_.filterId,filt_.IdLength);

					while (getline(is, line))
					{
						if (line[0]=='#') continue;
						//printf("%s\n",line.c_str());
						istringstream iss(line);
						if (!(iss >> fpga >> hyb >> apv)) { return 1; } // error
						for (int i = 0;i<filt_.CoefCount;i++)
							//if (!(iss >> filt_.filterData[fpga][hyb][apv][coefOrder[i]])) {return 1;};
							if (!(iss >> filt_.filterData[fpga][hyb][apv][i])) {return 1;};
					}
					is.close();
				}
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

	TCanvas         *c1;
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat("emrou");
	gStyle->SetPalette(1,0);
	gStyle->SetStatW(0.2);                
	gStyle->SetStatH(0.1);                
	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetPadLeftMargin(0.15);

	// Start X11 view
	//   TApplication theApp("App",NULL,NULL);

	// Root file is the first and only arg
	if ( argc-optind != 1 ) {
		cout << "Usage: meeg_baseline data_file\n";
		return(1);
	}

	for (int fpga = 0;fpga<7;fpga++)
		for (int hyb = 0;hyb<3;hyb++)
			for (int apv = 0;apv<5;apv++)
			{
				pedestal[fpga][hyb][apv] = 0;
				pedestalCount[fpga][hyb][apv] = 0;
				syncSize[fpga][hyb][apv] = 0;
				syncCount[fpga][hyb][apv] = 0;
				for (x=0; x < 35; x++) {
					sprintf(title,"sample_%d_F%d_H%d_A%d",x,fpga,hyb,apv);
					hist[fpga][hyb][apv][x] = new TH1S(title,title,10000,-0.5,1.5);
					histMin[fpga][hyb][apv][x] = 16384;
					histMax[fpga][hyb][apv][x] = -16384;
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
	ofstream outfile;
	//if (print_data)
	cout << "Writing calibration to " << inname+".filter" << endl;
	outfile.open(inname+".filter");
	outfile << inname << endl;

	//if (print_data)
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
	//if (print_data)
	cout << "Writing configuration to " <<outdir<<confname << endl;
	outconfig.open(outdir+confname);

	dataRead->next(&event);
	dataRead->dumpConfig(outconfig);
	outconfig.close();
	//dataRead.dumpStatus();

	runCount = atoi(dataRead->getConfig("RunCount").c_str());
	/*
	   double *moving_yi = new double[runCount];
	   double *moving_yi2 = new double[runCount];
	   */


	// Process each event
	eventCount = 0;

	do {
		int fpga = event.fpgaAddress();
		//printf("fpga %d\n",event.fpgaAddress());
		if (event.fpgaAddress()==7) 
		{
			//printf("not a data event\n");
			continue;
		}
		if (use_fpga!=-1 && fpga!=use_fpga) continue;
		//cout<<"  fpga #"<<event.fpgaAddress()<<"; number of samples = "<<event.count()<<endl;
		//if (print_data)
		if (eventCount%1000==0) printf("Event %d\n",eventCount);
		if (num_events!=-1 && eventCount >= num_events) break;
		if (read_temp && !event.isTiFrame()) for (uint i=0;i<4;i++) {
			printf("Event %d, temperature #%d: %f\n",eventCount,i,event.temperature(i,hybrid_type==1));
			read_temp = false;
		}
		for ( y=0; y < 6; y++ ) {
			for (int block=0; block < event.count()/128; block++) {
				int idx = 10000;
				double syncValue;
				int filterBuffer[10];
				int filterPtr = -1;
				for (int channel=0; channel < 128; channel++) {
					x = channel;
					// Get sample
					sample  = event.sample(block*128+channel);
					int hyb = sample->hybrid();
					int apv = sample->apv();
					if (use_hybrid!=-1 && sample->hybrid()!=use_hybrid) continue;
					if (use_apv!=-1 && sample->apv()!=use_apv) continue;
					//printf("hybrid %d\n",sample->hybrid());


					double adcValue = (double) (sample->value(y) & 0x3FFF);
					int adcValid = sample->value(y) >> 15;

					if (adcValid)
					{
						if (++filterPtr==10) filterPtr = 0;
						filterBuffer[filterPtr] = adcValue;

						if (use_filter) {
							if (x<10) continue;
							else {
								adcValue = 0;
								for (int i=0;i<10;i++) {
									//pedValue += filt_.filterData[fpga][hyb][apv][i]*filterBuffer[(filterPtr+i)%10];
									adcValue += filt_.filterData[fpga][hyb][apv][i]*filterBuffer[(filterPtr-i+10)%10];
								}
							}
						}
						double pedValue = adcValue;
						pedValue -= pedestal[fpga][hyb][apv];

						//printf("fpga: %d, x: %d, val: %x, adcValid: %d\n",event.fpgaAddress(),x, adcValue, adcValid);

						//if ( idx==10000 && adcValue > 0x2000 ) idx = 0;
						if ( pedValue > 0x2000 ) {
							syncValue = pedValue;

							double delta = syncValue-syncSize[fpga][hyb][apv];
							syncCount[fpga][hyb][apv]++;
							if (syncCount[fpga][hyb][apv]==1)
							{
								syncSize[fpga][hyb][apv] = syncValue;
							}
							else
							{
								syncSize[fpga][hyb][apv] += delta/syncCount[fpga][hyb][apv];
							}
							idx = 0;
							//printf("sync! event: %d, fpga: %d, hybrid: %d, apv: %d, samples: %d, x: %d, val: %x, adcValid: %d\n",eventCount,event.fpgaAddress(),sample->hybrid(),sample->apv(),y, x, adcValue, adcValid);
						}

						if ( idx < 35 ) {
							//cout << "Fill idx=" << dec << idx << " value=0x" << hex << value << endl;
							hist[fpga][hyb][apv][idx]->Fill(pedValue/syncValue);
							if ( pedValue < histMin[fpga][hyb][apv][idx] ) histMin[fpga][hyb][apv][idx] = pedValue;
							if ( pedValue > histMax[fpga][hyb][apv][idx] ) histMax[fpga][hyb][apv][idx] = pedValue;
						}
						if (idx > 20 && idx < 32) {
							pedestalCount[fpga][hyb][apv]++;
							double delta = adcValue-pedestal[fpga][hyb][apv];
							if (pedestalCount[fpga][hyb][apv]==1)
							{
								pedestal[fpga][hyb][apv] = adcValue;
							}
							else
							{
								pedestal[fpga][hyb][apv] += delta/pedestalCount[fpga][hyb][apv];
							}
						}
						idx++;
					}
				}
			}
		}
		eventCount++;

	} while ( dataRead->next(&event));
	dataRead->close();

	if (eventCount != runCount)
	{
		if (print_data)
			printf("ERROR: events read = %d, runCount = %d\n",eventCount, runCount);
	}



	for (int fpga = 0;fpga<7;fpga++)
	{
		if (use_fpga!=-1 && fpga!=use_fpga) continue;
		for (int hyb = 0;hyb<3;hyb++)
		{
			if (use_hybrid!=-1 && hyb!=use_hybrid) continue;
			for (int apv = 0;apv<5;apv++)
			{
				if (use_apv!=-1 && apv!=use_apv) continue;
				if (hist[fpga][hyb][apv][0]->GetEntries()<10) continue;
				if (!no_gui) {
					c1 = new TCanvas("c1","c1",1200,900);
					c1->Divide(7,5,0.0125,0.0125);
				}
				printf("F%d, H%d, A%d: pedestal %f, sync %f\n",fpga,hyb,apv,pedestal[fpga][hyb][apv],syncSize[fpga][hyb][apv]);
				double mean[35], rms[35];
				for (x=0; x < 35; x++) 
				{
					mean[x] = hist[fpga][hyb][apv][x]->GetMean();
					rms[x] = hist[fpga][hyb][apv][x]->GetRMS();
					//	for (int i=-16384;i<16384;i++) hist[x]->Fill(i);
				}

				for (x=0; x < 35; x++) {
					TF1 *gaus = new TF1("gaus"+x,"gaus",-16384,16384);
					gaus->FixParameter(0,0.0);
					gaus->SetParameter(1,mean[x]);
					gaus->SetParameter(2,rms[x]);
					if (no_gui)
						hist[fpga][hyb][apv][x]->Fit(gaus,"QL0","",mean[x]-3*rms[x],mean[x]+3*rms[x]);
					else {
						c1->cd(x+1);
						hist[fpga][hyb][apv][x]->Fit(gaus,"QL","",mean[x]-3*rms[x],mean[x]+3*rms[x]);
						//	hist[fpga][hyb][apv][x]->GetXaxis()->SetRangeUser(histMin[fpga][hyb][apv][x]-100,histMax[fpga][hyb][apv][x]+100);
						hist[fpga][hyb][apv][x]->GetXaxis()->SetRangeUser(mean[x]-5*rms[x],mean[x]+5*rms[x]);
						hist[fpga][hyb][apv][x]->Draw();
					}
					if (print_data)
						printf("mean %f, RMS %f, mode %f, fitted mean %f, fitted RMS %f\n",mean[x],rms[x],hist[fpga][hyb][apv][x]->GetBinCenter(hist[fpga][hyb][apv][x]->GetMaximumBin()),gaus->GetParameter(1),gaus->GetParameter(2));
					//meanVal[x]  = hist[x]->GetFunction("gaus")->GetParameter(1);
					//meanVal[x]  = hist[x]->GetMean() / 16383.0;
					plotX[x] = x;
					
					meanVal[(x+n_delay)%35]  = gaus->GetParameter(1);
					if (x==0) meanVal[(x+n_delay)%35]  = mean[x];
					
					//meanVal[(x+n_delay)%35]  = hist[fpga][hyb][apv][x]->GetMean();
					//meanVal[(x+n_delay)%35]  = hist[fpga][hyb][apv][x]->GetBinCenter(hist[fpga][hyb][apv][x]->GetMaximumBin());
					//meanVal[x]  = hist[fpga][hyb][apv][x]->GetBinCenter(hist[fpga][hyb][apv][x]->GetMaximumBin())+pedestal[fpga][hyb][apv]/syncSize[fpga][hyb][apv];
				}

				double avg = 0;
				for (x=22; x < 32; x++) {
					avg += meanVal[x];
				}
				avg/=10;

				double adjVal[35];
				for (x=0; x < 35; x++) {
					adjVal[x] = meanVal[x];
					//adjVal[x] = (meanVal[x]-avg)/(meanVal[0]-avg);
					//if (x!=n_delay) adjVal[x]*=-1; //invert the coefficients
				}

				double data[2*n_coeffs];
				gsl_fft_complex_wavetable * wavetable;
				gsl_fft_complex_workspace * workspace;

				for (int i = 0; i < n_coeffs; i++)
				{
					//REAL(data,(i+5)%n) = meanVal[i];
					REAL(data,i) = adjVal[i];
					IMAG(data,i) = 0.0;
					//					printf ("%d: %e %e\n", i, REAL(data,i), IMAG(data,i));
				}

				wavetable = gsl_fft_complex_wavetable_alloc (n_coeffs);
				workspace = gsl_fft_complex_workspace_alloc (n_coeffs);
				gsl_fft_complex_forward (data, 1, n_coeffs, wavetable, workspace);
				for (int i = 0; i < n_coeffs; i++)
				{
					//printf ("%d: %e %e\n", i, REAL(data,i), IMAG(data,i));
				}

				/*
				double last_filter[2*n];
				if (use_filter) {
					for (int i = 0; i < n; i++)
					{
						//REAL(data,(i+5)%n) = meanVal[i];
						REAL(last_filter,i) = filt_.filterData[fpga][hyb][apv][i];
						IMAG(last_filter,i) = 0.0;
						//					printf ("%d: %e %e\n", i, REAL(data,i), IMAG(data,i));
					}
					gsl_fft_complex_forward (last_filter, 1, n, wavetable, workspace);
				}
				*/


				for (int i = 0; i < n_coeffs; i++)
				{
					gsl_complex temp = gsl_complex_rect (REAL(data,i), IMAG(data,i));
					temp = gsl_complex_inverse(temp);
					temp = gsl_complex_mul(temp,gsl_complex_polar(1,-2*n_delay*i*M_PI*2.0/n_coeffs));
					/*
					   if (use_filter) {
					   gsl_complex temp2 = gsl_complex_rect (REAL(last_filter,i), IMAG(last_filter,i));
					   temp = gsl_complex_mul(temp,temp2);
					   }
					   */
					REAL(data,i) = GSL_REAL(temp);
					IMAG(data,i) = GSL_IMAG(temp);
					//printf ("%d: %e %e\n", i, REAL(data,i), IMAG(data,i));
				}

				//REAL(data,0) = 1.0;

				gsl_fft_complex_inverse (data, 1, n_coeffs, wavetable, workspace);
				for (int i = 0; i < n_coeffs; i++)
				{
					//					printf ("%d: %e %e\n", i, REAL(data,i), IMAG(data,i));
					adjVal[i]  = REAL(data,i);
					//meanVal[i]  = REAL(data,(i+5)%n);
				}

				/*
				for (x=0; x < 35; x++) {
					adjVal[x] = meanVal[x];
					if (x!=n_delay) adjVal[x]*=-1; //invert the coefficients
				}
				*/

				double          sum = 0;
				for (x=0; x < n_coeffs; x++) {
					sum += adjVal[x];
					if (print_data)
						printf("Idx=%d value=%f adj=%f\n",x,meanVal[x],adjVal[x]);
				}
				//if (print_data)
				cout << "Sum=" << sum << endl;

				outfile << fpga << "\t" << hyb << "\t" << apv;
				for (x=0; x < n_coeffs; x++) {
					outfile << "\t" << adjVal[x];
				}
				outfile << endl;

				if (!no_gui)
				{
					c1->cd();
					sprintf(name,"%s_fits_F%d_H%d_A%d.png",inname.Data(),fpga,hyb,apv);
					c1->SaveAs(name);
					delete c1;
					c1 = new TCanvas("c1","c1");
					c1->cd();
					TGraph          *plot;
					TGraph          *plot2;
					TMultiGraph          *mg = new TMultiGraph();
					plot = new TGraph(n_coeffs,plotX,adjVal);
					plot2 = new TGraph(35,plotX,meanVal);
					plot2->SetMarkerColor(2);
					mg->Add(plot);
					mg->Add(plot2);
					mg->Draw("a*");
					sprintf(name,"%s_coeffs_F%d_H%d_A%d.png",inname.Data(),fpga,hyb,apv);
					c1->SaveAs(name);
					delete plot;
					delete plot2;
					delete mg;
					delete c1;
				}
			}
		}
	}

	// Close file
	outfile.close();
	delete dataRead;
	return(0);
}

