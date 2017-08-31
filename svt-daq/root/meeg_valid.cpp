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
#include <meeg_utils.hh>
#include <stdarg.h>
#include <DevboardEvent.h>
#include <DevboardSample.h>
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
	bool read_temp = true;
	int hybrid_type = 0;
	bool evio_format = false;
	int fpga = -1;
	int hybrid = -1;
	int num_events = -1;
	int c;

	DataRead        *dataRead;
	DevboardEvent    event;
	DevboardSample   *sample;
	uint            x;
	uint            y;
	int            eventCount;
	int runCount;

	while ((c = getopt(argc,argv,"ht:H:F:e:Ed")) !=-1)
		switch (c)
		{
			case 'h':
				printf("-h: print this help\n");
				printf("-d: turn on debug\n");
				printf("-t: hybrid type (1 for old test run hybrid, 2 for new 2014 hybrid)\n");
				printf("-F: use only specified FPGA\n");
				printf("-H: use only specified hybrid\n");
				printf("-e: stop after specified number of events\n");
				printf("-E: use EVIO file format\n");
				return(0);
				break;
			case 't':
				hybrid_type = atoi(optarg);
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

	// Root file is the first and only arg
	if ( argc-optind < 1 ) {
		cout << "Usage: meeg_valid data_file\n";
		return(1);
	}


	while (optind<argc)
	{
		cout << "Reading data file " <<argv[optind] << endl;
		// Attempt to open data file
		if ( ! dataRead->open(argv[optind]) ) return(2);


		dataRead->next(&event);
		runCount = atoi(dataRead->getConfig("RunCount").c_str());
		int max_count = runCount==0 ? 10000 : runCount;

		read_temp = true;

		if(debug) printf("%d events expected from config\n",runCount);


		// Process each event
		eventCount = 0;

		do {
			if(debug) printf("Event %d\n",eventCount);

			if(debug) printf("fpga %d\n",event.fpgaAddress());
			if (event.fpgaAddress()==7) 
			{
				printf("not a data event\n");
				continue;
			}
			if (fpga!=-1 && ((int)event.fpgaAddress())!=fpga) continue;
			if(debug) cout<<"  fpga #"<<event.fpgaAddress()<<"; number of samples = "<<event.count() << " is" <<endl;
			if(debug) cout<<"  fpga #"<<event.fpgaAddress()<<"; number of samples = "<<event.count()<<endl;
			//if (eventCount%1000==0) printf("Event %d\n",eventCount);
			if (num_events!=-1 && eventCount >= num_events) break;
			if (read_temp && !event.isTiFrame()) {
				//printf("Temperatures, event %d:\t",eventCount);
				for (uint i=0;i<4;i++) {
					printf("temperature #%d: %f\t",i,event.temperature(i,hybrid_type==1));
				}
				//printf("\n");
				read_temp = false;
			}

			for (x=0; x < event.count(); x++) {
				// Get sample
				sample  = event.sample(x);
				if(debug) printf("event %d\tx=%d\tF%d H%d A%d channel %d, samples:\t%d\t%d\t%d\t%d\t%d\t%d\n",eventCount,x,event.fpgaAddress(),sample->hybrid(),sample->apv(),sample->channel(),sample->value(0),sample->value(1),sample->value(2),sample->value(3),sample->value(4),sample->value(5));
				if (hybrid!=-1 && ((int)sample->hybrid())!=hybrid) continue;
				//printf("hybrid %d\n",sample->hybrid());

				bool bad_event = false;
				if (sample->apv()<0 || sample->apv()>4) {
					cout << "Apv " << dec << sample->apv() << " out of range" << endl;
				}
				if (sample->channel()<0 || sample->channel()>127) {
					cout << "Channel " << dec << sample->channel() << " out of range" << endl;
				}

				// Filter APVs
				if ( eventCount >= 20 ) {
					for ( y=0; y < 6; y++ ) if (sample->value(y)==0) {
						printf("sample is zero\n");
						bad_event = true;
					}
				}
				if (bad_event) continue;
			}
			if (eventCount<max_count) {
			}
			eventCount++;

		} while ( dataRead->next(&event));
		dataRead->close();
		printf("events read = %d, runCount = %d\n",eventCount, runCount);
		optind++;
	}

	delete dataRead;
	return(0);
}

