//-----------------------------------------------------------------------------
// File          : DataRead.h
// Author        : Ryan Herbst  <rherbst@slac.stanford.edu>
// Created       : 04/12/2011
// adapted by mgraham on 4/2/2012 to read in evio files
// Project       : General Purpose
//-----------------------------------------------------------------------------
// Description :
// Read data & configuration from disk
//-----------------------------------------------------------------------------
// Copyright (c) 2011 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 04/12/2011: created
//-----------------------------------------------------------------------------
#ifndef __DATA_READ_EVIO_H__
#define __DATA_READ_EVIO_H__

#include <evio.h>
#include <string>
#include <map>
#include <vector>
#include <sys/types.h>
#include <DataRead.h>
#include <Data.h>
using namespace std;
#define MAXEVIOBUF   1000000

// Define variable holder
typedef map<string,string> VariableHolder;

//! Class to contain generic register data.
class DataReadEvio : public DataRead {
	bool debug_;

	int maxbuf ;
	Data *fpga_banks[16];
	int fpga_count, fpga_it;
	int svt_bank_min,svt_bank_range;
    
    bool is_engrun;
    int svt_data_tag;
    int svt_ti_data_tag;
    int svt_ti_data_size;
    int svt_config_tag;


	void parse_event( unsigned int *buf);
	void parse_eventBank( unsigned int *buf,int bank_length);
	void parse_SVTBank( unsigned int *buf,int bank_length);
	void parse_ECalBank( unsigned int *buf,int bank_length);
	void parse_ECalCompositeData( unsigned int *buf,int bank_length);

	void eventInfo(unsigned int *buf);

	int getFragType(int type);

	enum {
		BANK = 0,
		SEGMENT,
		TAGSEGMENT,
		UINT32,
		COMPOSITE,
		CHARSTAR8,
		UNEXPECTED,
	};
	int fragment_offset[3];

	int evtType;
	unsigned short evtTag;
	unsigned char evtNum;

	public:

	//! Constructor
	DataReadEvio ( );

	//! Deconstructor
	~DataReadEvio ( );

	void set_engrun(bool engrun);
	void set_bank_num(int bank_num);
	void set_bank_range(int bank_num);

	bool open ( string file, bool compressed = false );

	void close();

	//! Get next data record
	/*! 
	 * Returns true on success
	 * \param data Data object to store data
	 */
	bool  next ( Data *data);

	//! Get next data record & create new data object
	/*! 
	 * Returns NULL on failure
	 */
	//      Data *next ( );

};
#endif
