//-----------------------------------------------------------------------------
// File          : readfile_test.cpp
// Author        : Per Hansson Adrian  <phansson@slac.stanford.edu>
// Created       : 08/25/2017
// Project       : HPS
//-----------------------------------------------------------------------------
// Description :
// File to test reading binary data file.
//-----------------------------------------------------------------------------
// Copyright (c) 2009 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 08/25/2017: created
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
#include <sys/stat.h>
#include <fcntl.h>

using namespace std;


void xmlParse(int* fd, uint size, char* data) {
  char         *buff;
  uint         mySize;
  uint         myType;
  
  // Decode size
  myType = (size >> 28) & 0xF;
  mySize = (size & 0x0FFFFFFF);

   // Read file
   buff = (char *) malloc(mySize+1);
   size_t bytes_read;
   bytes_read = ::read(*fd, buff, mySize);
   cout << "DataRead::xmlParse read " << bytes_read << " bytes, asked for " << mySize << endl;
   if (((uint)bytes_read) != mySize) {
     cout << "DataRead::xmlParse -> Read error!" << endl;
     return;
   } else {
     cout << "DataRead::xmlParse -> Read ok into buff!" << endl;       
   }
   free(buff);
   
}


// Process the data
// Pass root file to open as first and only arg.
int main ( int argc, char **argv ) {

  int fd;
  size_t bytes_read;
  size_t buf_size = 4;
  //char buf[buf_size];
  uint buf;
  uint ibuf;
  

   // input file is the first and only arg
   if ( argc != 2 ) {
      cout << "Usage: data_file\n";
      return(1);
   }

   // Attempt to open data file

   cout << "Try to open data file: " << argv[1] << endl;

   fd = open(argv[1],O_RDONLY);

   cout << "fd: " << fd << endl;

   if( fd < 0) {

     cout << "Failed to open file, result: " << fd << endl;
     //cout << "errno " << errno << endl;

     return 1;
   } 

   cout << "Successfully opened file" << endl;
   
   // Attempt to read from the file

   bool readOK = true;
   bool found = false;
   char* shBuff = NULL;

   do {
     
     bytes_read = read(fd, &buf, buf_size);
     
     cout << "Read " << bytes_read << " bytes" << endl;
     
     cout << "\"" << buf << "\"" << endl;
     
     cout << "fd: " << fd << endl; 
     
     if(bytes_read != 4) {
       cout << "only read " << bytes_read << " bytes" << endl;
       break;
     } 
     
     uint frame_type;
     frame_type = (buf >> 28) & 0xF;
     
     cout << "frame_type: " << frame_type << endl;
     
     
     // Frame type
     switch ( frame_type ) {
       
       // Data
     case Data::RawData : found = true; break;
       
       // Configuration
     case Data::XmlConfig : xmlParse(&fd, buf,shBuff); break;
       
       // Status
     case Data::XmlStatus : xmlParse(&fd, buf,shBuff); break;
       
       // Start
     case Data::XmlRunStart : xmlParse(&fd, buf,shBuff); break;
       
       // Stop
     case Data::XmlRunStop : xmlParse(&fd, buf,shBuff); break;
       
       // Time
     case Data::XmlRunTime : xmlParse(&fd, buf,shBuff); break;
       
       // Unknown
     default: 
       cout << "DataRead::next -> Unknown data type 0x" 
            << hex << setw(8) << setfill('0') << ((buf >> 28) & 0xF) << " skipping." << endl;
       break;
     }
     
   } while( ! found);
   
   
   cout << "found raw data frame of size " << buf << endl;
   
   DevboardEvent event;
   bool dataReadStatus;
   dataReadStatus = event.read(fd, buf);
   cout << "FPGA: 0x"  << hex << event.fpgaAddress() << endl;
   cout << "Sample count: " << event.count() << endl;
   for(uint x=0; x!=event.count(); x++) {
     DevboardSample *sample = event.sample(x);
     cout << dec << "x: " << x
          << " fpga: " << sample->fpgaAddress()
          << " hyb: " << sample->hybrid()
          << " apv: " << sample->apv()
          << " channel: " << sample->channel();
     for(int ival=0;ival!=6; ++ival) {
       cout << " " << (sample->value(ival) & 0x3FFF);
     }
     cout << endl;
     
   }
     


   // Close file
   close(fd);
   
   return 0;
}

