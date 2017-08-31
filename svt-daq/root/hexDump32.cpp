#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;

int main ( int argc, char **argv ) {
   uint val;  
   int fd;
 
   if (argc != 2 ) {
     cout << "Usage: " << argv[0] << " data_file" << endl;
     return 1;
   }

   cout << "dumping file: " << argv[1] << endl;

   fd = open ( argv[1], 0 );
   
   if( fd < 0) {
     cout << "Failed to open file " << argv[1] << endl;
   }
   
   while ( read(fd,&val,4) == 4 ) cout << "0x" << setfill('0') << setw(8) << hex << val << endl;
   
} 
