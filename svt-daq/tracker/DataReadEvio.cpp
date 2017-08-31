//-----------------------------------------------------------------------------
// File          : DataReadEvio.cpp
// Author        : Ryan Herbst  <rherbst@slac.stanford.edu>
// Created       : 04/12/2011
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

#include <DataReadEvio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <iostream>
#include <iomanip>
using namespace std;

// Constructor
DataReadEvio::DataReadEvio ( ) {
	debug_=false;
	//debug_=true;
	fd_ = -1;
	maxbuf=MAXEVIOBUF;
	fragment_offset[0] = 2;//BANK
	fragment_offset[1]=1;//SEGMENT
	fragment_offset[2]=1;//TAGSEGMENT
	fpga_count = 0;
	fpga_it = 0;
	svt_bank_min = 3;
	svt_bank_range = 1;
    svt_data_tag = 3;
    svt_ti_data_tag = 57610;
    svt_ti_data_size = 4;
    svt_config_tag = 57614;
    is_engrun = false;
}

// Deconstructor
DataReadEvio::~DataReadEvio ( ) { }

void DataReadEvio::set_engrun(bool engrun) {
	is_engrun = engrun;
    if (engrun) {
        svt_bank_min = 51;
        svt_bank_range = 16;
    } else {
        svt_bank_min = 3;
        svt_bank_range = 1;
    }
}

void DataReadEvio::set_bank_num(int bank_num) {
    svt_bank_min = bank_num;
}

void DataReadEvio::set_bank_range(int bank_num) {
    svt_bank_range = bank_num;
}

// Open file
bool DataReadEvio::open ( string file, bool compressed ) {
    int status;
    char * filename = (char *) malloc((file.size()+1)*sizeof(char));
    strcpy(filename,file.c_str());
    if((status=evOpen(filename,"r",&fd_))!=0) {
        //      printf("\n ?Unable to open file %s, status=%d\n\n",file,status);
        cout<<"Unable to open file "<<file<<", status="<<status<<endl;
        return(false);
    }
    return(true);
}

void DataReadEvio::close () {
    evClose(fd_);
    fd_ = -1;
}

bool DataReadEvio::next(Data *data) {
    if(debug_)printf("fpga_count %d, fpga_it %d\n",fpga_count, fpga_it);
    if (fpga_count>fpga_it)
    {
        if(debug_)printf("pulling a bank out of cache\n");
        Data *source_data = fpga_banks[fpga_it++];
        data->copy(source_data->data(),source_data->size());
        if (fpga_it==fpga_count)
        {
            for (int i=0;i<fpga_count;i++)
            {
                delete fpga_banks[i];
            }
            fpga_count = 0;
            fpga_it = 0;
        }
        return true;
    }
    if(debug_)printf("reading an EVIO event\n");
    bool nodata = true;
    int nevents=0;
    int status;
    if ( fd_ < 0 ) {
        cout<<"DataReadEvio::next error fd_<0...no file open"<<endl;
        return(false);
    }

    do{    
        unsigned int *buf = (unsigned int*)malloc(maxbuf*sizeof(unsigned int));
        status = evRead(fd_,buf,maxbuf);
        if(status==S_SUCCESS){
            nevents++;
            //  here, get the offset and the length of the SVT data in the buffer (buf)
            //  do this by stepping through the banks until get an SVT bank
            //  then hand make a "Data" object 
            eventInfo(buf);
            if(debug_)printf("evtTag = %d\n",evtTag);
            if(evtTag>=32 || evtTag<16){
                parse_event(buf);
                //fpga_it = fpga_banks.begin();
                nodata=false;  
                free(buf);
            }else{
                if(evtTag==20)return(false); //this is the end of data
                //otherwise, just skip it. 
                cout<<"Not a data event...skipping"<<endl;
                free(buf);
            }
        } else if (status==EOF)
        {
            cout << "end of file" << endl;
            free(buf);
            return(false);
        }else{
            cout<<"oops...broke trying to evRead; error code "<<status<<endl;

            free(buf);
            return(false);
        }
    }  while(nodata);
    if(debug_)cout<<"Found a data event"<<endl;
    return(next(data));
}

void DataReadEvio::eventInfo(unsigned int *buf) {
    evtTag         = (buf[1]>>16)&0xffff;
    evtType        = (buf[1]>>8)&0x3f;
    evtNum         = buf[1]&0xff;
}

void DataReadEvio::parse_event(unsigned int *buf) {
    int length,type, padding=0;
    unsigned short tag;
    unsigned short num;
    length      = buf[0]+1;
    tag         = (buf[1]>>16)&0xffff;
    padding     = (buf[1]>>14)&0x3;
    type        = (buf[1]>>8)&0x3f;
    num         = buf[1]&0xff;

    if(debug_)
    {
        cout<<"Spitting out event header data: "<<endl;
        cout<<"length: "<<length<<", tag: "<<tag<<endl;
        cout<<"padding: "<<padding<<", type: "<<type<<", num: "<<num<<endl;
    }
    int fragType = getFragType(type);
    if (fragType!=BANK) printf("data type of event should be BANK but was %d\n",type);


    parse_eventBank(&buf[fragment_offset[fragType]],length-fragment_offset[fragType]);
}

void DataReadEvio::parse_eventBank(unsigned int *buf, int bank_length) {
    int ptr = 0;
    int length,type, padding=0;
    unsigned short tag;
    unsigned short num;
    while (ptr<bank_length) {
        if (debug_) printf("ptr = %d, bank_length = %d\n",ptr,bank_length);
        length      = buf[ptr]+1;
        tag         = (buf[ptr+1]>>16)&0xffff;
        type        = (buf[ptr+1]>>8)&0x3f;
        padding     = (buf[ptr+1]>>14)&0x3;
        num         = buf[ptr+1]&0xff;
        if(debug_)
        {
            cout<<"Spitting out event bank header data: "<<endl;
            cout<<"length: "<<length<<", tag: "<<tag<<endl;
            cout<<"padding: "<<padding<<", type: "<<type<<", num: "<<num<<endl;
        }
        int fragType = getFragType(type);
        if (fragType==BANK)
        {
            if (tag>=svt_bank_min && tag<svt_bank_min+svt_bank_range) {
                if (debug_) printf("found SVT bank, tag %d\n",tag);
                parse_SVTBank(&buf[ptr+2],length-2);
            }
            /*else switch (tag) {
              case 1:
              if (debug_) printf("ECal top bank\n");
              parse_ECalBank(&buf[ptr+2],length-2);
              break;
              case 2:
              if (debug_) printf("ECal bottom bank\n");
              parse_ECalBank(&buf[ptr+2],length-2);
              break;
              default:
              printf("Unexpected bank tag %d\n",tag);
              break;
              }*/
        }
        else
            if (debug_) printf("looking for event bank of type BANK but found %d\n",type);
        ptr+=length;
    }
}

void DataReadEvio::parse_SVTBank(unsigned int *buf, int bank_length) {
    int ptr = 0;
    int length,type, padding=0;
    unsigned short tag;
    unsigned short num;
    // Pointer to the data structure that will be passed out
    Data* tb = NULL; 
    // Pointer to the TI data for this bank
    uint* tiDataPtr = NULL;
    uint tiDataLen = 0;
    uint debug_local = 0;
    if(debug_local){
      cout<<"====\nparse SVT bank: bank length: "<<bank_length<<endl;
    }

    while (ptr<bank_length) {
        if (debug_) printf("ptr = %d, bank_length = %d\n",ptr,bank_length);
        length      = buf[ptr]+1;
        tag         = (buf[ptr+1]>>16)&0xffff;
        type        = (buf[ptr+1]>>8)&0x3f;
        padding     = (buf[ptr+1]>>14)&0x3;
        num         = buf[ptr+1]&0xff;
        if(debug_local)
        {
            cout<<"Spitting out SVT bank header data: "<<endl;
            cout<<"length: "<<length<<", tag: "<<tag<<endl;
            cout<<"padding: "<<padding<<", type: "<<type<<", num: "<<num<<endl;
        }
        int fragType = getFragType(type);
        
        
        if (fragType==UINT32 && (!is_engrun || tag==svt_data_tag))
        {
            if(debug_) {
              printf("Got data\n");
            }
            tb=new Data();
            if (!is_engrun){
                uint *data_  = (uint *)malloc((length-1) * sizeof(uint));
                memcpy(data_+1,&buf[ptr+2],(length-2)*sizeof(uint));
                data_[0] = tag;
                if (tag==7) data_[0]+=0x80000000;
                tb->copy(data_,length-1);
                free(data_);
            } else {
                tb->copy(&buf[ptr+2],length-2);
            }
            fpga_banks[fpga_count++] = tb;
        }
        else if (fragType==UINT32 && (!is_engrun || tag==svt_ti_data_tag))  {
          tiDataLen = length-2;
          if(is_engrun && tiDataLen!=svt_ti_data_size) {
            printf("DataReadEvio: Error: TI data length is wrong!? (%d)\n",tiDataLen);
            exit(1);            
          }
          if(debug_local) printf("Got TI data bank with data %d words\n",tiDataLen);
          tiDataPtr = (uint*) malloc(tiDataLen * sizeof(uint));
          if(debug_local) printf("Allocated TI data space at %p\n",tiDataPtr);
          memcpy(tiDataPtr,&buf[ptr+2],tiDataLen*sizeof(uint));
          if(debug_local) printf("Done memcpy\n");
        }
        else if (fragType==CHARSTAR8 && tag==svt_config_tag)  {
          printf("Found config/status bank %d in length(%lu %lu)\n",length,sizeof(uint),sizeof(char));
           char* str = (char*) malloc((length-2)*sizeof(uint)*sizeof(uint)/sizeof(char));
           memcpy(str,&buf[ptr+2],(length-2)*sizeof(uint)*sizeof(uint)/sizeof(char));
           printf("Done copy\n");
           printf("\"%s\"\n",str);
           printf("Done\n");
           free(str);
        } 
        else if (debug_)
            printf("data type of SVT bank should be UINT32 but was %d\n",type);
        ptr+=length;
    }

    // Create a new Data structure that includes the TI data if it's found
    // way to complicated but whatever
    if( tb!=NULL ) {
      if(tiDataPtr!=NULL) {
        if(debug_local) printf("Copy the TI data to the data object");
        uint len = tb->size()+tiDataLen;
        uint *d  = (uint *)malloc(len * sizeof(uint));
        if(debug_local) printf("Allocated %d for the new data\n", len);
        memcpy(d,tb->data(),tb->size()*sizeof(uint));
        memcpy(d+tb->size(),tiDataPtr,tiDataLen*sizeof(uint));
        if(debug_local) printf("Done memcpy\n");
        tb->copy(d,len);
        if(debug_local) printf("Done copy the TI data to the data object");
        
      } else {
        if(is_engrun) {
          if(debug_local) printf("Inject empty TI data");
          uint len = tb->size()+svt_ti_data_size;
          uint *d  = (uint *)malloc(len * sizeof(uint));
          if(debug_local) printf("Allocated %d for the new data\n", len);
          memcpy(d,tb->data(),tb->size()*sizeof(uint));
          memset(d+tb->size(),0xff,svt_ti_data_size*sizeof(uint));
          if(debug_local) printf("Done memcpy\n");
          tb->copy(d,len);
          if(debug_local) printf("Done copy the fake TI data to the data object");        
        }
      }
    } else {
      if(debug_local) printf("No tb object built. free TI data if it's there\n");
      if(tiDataPtr!=NULL) {
        free(tiDataPtr);
      }
    }

    if(debug_local) {
      cout<<"\n DONE parsing SVT bank\n===="<<endl;
    }
}

void DataReadEvio::parse_ECalBank(unsigned int *buf, int bank_length) {
    int ptr = 0;
    int length,type, padding=0;
    unsigned short tag;
    unsigned short num;
    while (ptr<bank_length) {
        if (debug_) printf("ptr = %d, bank_length = %d\n",ptr,bank_length);
        length      = buf[ptr]+1;
        tag         = (buf[ptr+1]>>16)&0xffff;
        type        = (buf[ptr+1]>>8)&0x3f;
        padding     = (buf[ptr+1]>>14)&0x3;
        num         = buf[ptr+1]&0xff;
        if(debug_)
        {
            cout<<"Spitting out ECal bank header data: "<<endl;
            cout<<"length: "<<length<<", tag: "<<tag<<endl;
            cout<<"padding: "<<padding<<", type: "<<type<<", num: "<<num<<endl;
        }
        int fragType = getFragType(type);

        if (fragType==COMPOSITE)
        {
            parse_ECalCompositeData(&buf[ptr+2],length-2);
        }
        else
            printf("data type of ECal bank should be COMPOSITE but was %d\n",type);
        ptr+=length;
    }
}

void DataReadEvio::parse_ECalCompositeData(unsigned int *buf, int bank_length) {
    int ptr = 0;
    int length,type, padding=0;
    unsigned short tag;
    unsigned short num;
    int fragType; 

    if (debug_) printf("ptr = %d, bank_length = %d\n",ptr,bank_length);
    tag         = (buf[ptr]>>20)&0xfff;
    type        = (buf[ptr]>>16)&0xf;
    length      = (buf[ptr]&0xffff) + 1;
    fragType = getFragType(type);
    if(debug_)
    {
        cout<<"Spitting out ECal composite data sub-TAGSEGMENT header data: "<<endl;
        cout<<"length: "<<length<<", tag: "<<tag<<", type: "<<type<<endl;
        for (int i=1;i<length;i++)
        {
            char c;
            //printf("%x\n",buf[ptr+i]);
            c = (buf[ptr+i])&0xff;
            printf("%c",c);
            c = (buf[ptr+i]>>8)&0xff;
            printf("%c",c);
            c = (buf[ptr+i]>>16)&0xff;
            printf("%c",c);
            c = (buf[ptr+i]>>24)&0xff;
            printf("%c",c);
        }
        printf("\n");
    }

    ptr+=length;
    if (debug_) printf("ptr = %d, bank_length = %d\n",ptr,bank_length);
    length      = buf[ptr]+1;
    tag         = (buf[ptr+1]>>16)&0xffff;
    type        = (buf[ptr+1]>>8)&0x3f;
    padding     = (buf[ptr+1]>>14)&0x3;
    num         = buf[ptr+1]&0xff;
    if(debug_)
    {
        cout<<"Spitting out ECal composite data sub-BANK header data: "<<endl;
        cout<<"length: "<<length<<", tag: "<<tag<<endl;
        cout<<"padding: "<<padding<<", type: "<<type<<", num: "<<num<<endl;
    }
}

int DataReadEvio::getFragType(int type){
    switch (type) {
        case 0x1:
            if(debug_)cout<<"found UINT32 data"<<endl;
            return(UINT32);
            break;
        case 0x3:
            if(debug_)cout<<"found CHARSTAR8 data"<<endl;
            return(CHARSTAR8);
            break;
        case 0xe:
        case 0x10:
            if(debug_)cout<<"found BANK/ALSOBANK data"<<endl;
            return(BANK);
            break;
        case 0xf:
            if(debug_)cout<<"found COMPOSITE data"<<endl;
            return(COMPOSITE);
            break;
        case 0xd:
        case 0x20:
            if(debug_)cout<<"found SEGMENT/ALSOSEGMENT data"<<endl;
            return(SEGMENT);
            break;
        case 0xc:
            if(debug_)cout<<"found TAGSEGMENT data"<<endl;
            return(TAGSEGMENT);
            break;
        default:
            printf("Unexpected fragment type: %d",type);
            return(UNEXPECTED);
            break;
    }
}
