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

   bytes_read = read(fd, &buf, buf_size);

   cout << "Read " << bytes_read << " bytes" << endl;
   
   cout << "\"" << buf << "\"" << endl;

   cout << "fd: " << fd << endl; 

   




   // Close file
   close(fd);
   
   return 0;
}

