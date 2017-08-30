
/* include files */
#include <evio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <expat.h>
#include <unistd.h>

#define MAXEVIOBUF   10000000

static int maxbuf         = MAXEVIOBUF;
static unsigned int evtTag;
static unsigned int evtType;
static unsigned int evtNum;
static unsigned int evtPad;

static char* eviofilename= (char*) "data0.evio";
static bool debug=false;
enum {
   BANK = 0,
   SEGMENT,
   TAGSEGMENT,
   UINT32,
   COMPOSITE,
   CHARSTAR8,
   UNEXPECTED,
};

typedef enum {
   DPM0=11,
   DPM1=12,
   DPM2=13,
   DPM3=14,
   DPM4=15,
   DPM5=16,
   DPM6=17,   
   DPM7=19,
} RCETAGS;



static int fragment_offset[] = {2 /*BANK*/, 1 /* SEGMENT */, 1 /* TAGSEGMENT */};


void eventInfo(unsigned int* buf);
int getFragType(int type);
void check_lfsr(unsigned int *buf, int length);
void parse_eventBank(unsigned int *buf, int bank_length);
void parse_event(unsigned int *buf);
unsigned int flfsr32(unsigned int input);
unsigned short int flfsr16(unsigned short int input);
unsigned short int lfsr16(unsigned short int input);


int main (int argc, char **argv) {

   int status, handle, nevent;
   int c;
   while( (c = getopt(argc,argv,"f:d")) != -1) {
      switch (c) 
      {
         case 'f':
            eviofilename = optarg;
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
   }
   
   
  

  /* open evio (binary) output file */
   if((status=evOpen(eviofilename,"r",&handle))!=0) {
      printf("Unable to open evio output file %s, status=%d\n\n",eviofilename,status);
      exit(1);
   } else {
      printf("Opened evio output file %s, status=%d\n\n",eviofilename,status);     
   }
   
   /* event loop */
   nevent=0;
   status=0;
   while((status==0) ) {
      ++nevent;      
      unsigned int *buf = (unsigned int*)malloc(maxbuf*sizeof(unsigned int));
      if((buf==NULL)) {
         int sz=maxbuf*sizeof(unsigned int);
         printf("\n   *** Unable to allocate buffers ***\n\n");
         printf("\n buf size=%d bytes, addr=0x%p \n",sz,buf);
         exit(1);
      }      
      status=evRead(handle,buf,maxbuf);
      eventInfo(buf);
      if(debug) {
         printf("event %d\n",nevent);
         printf("%d %d %d %d\n",evtTag,evtPad,evtType,evtNum);
      }
      
      if(evtTag==1) {
         if(debug) printf("found data event\n");
         parse_event(buf);
         free(buf);
      }
      else {
         if(evtTag==20) {
            if(debug) printf("end of data evtTag found\n");
            free(buf);
            return false;
         }
         else {
            printf("not a data event. skip evtTag=%d\n",evtTag);
            free(buf);
         }
      }      
   }
   if(status!=EOF) printf("\n   *** error reading file, status is: 0x%x ***\n\n",status);
   else printf("EOF: status=%d (EFO=%d)\n\n",status, EOF);
   
   
   /* done */
   evClose(handle);
  exit(0);
  
}


void eventInfo(unsigned int* buf) {
//    evtTag         = (buf[1]>>16)&0xffff;
//    evtType        = (buf[1]>>8)&0x3f;
//    evtNum         = buf[1]&0xff;
//    evtTag         = 0;
  evtTag         = (buf[1]>>16)&0xffff;
  evtType        = (buf[1]>>8)&0x3f;
  evtPad         = (buf[1]>>14)&0x3;
  evtNum         = (buf[1])&0xff; 
}


int getFragType(int type){
	switch (type) {
		case 0x1:
			if(debug)printf("found UINT32 data\n");
			return(UINT32);
			break;
		case 0x3:
			if(debug)printf("found CHARSTAR8 data\n");
			return(CHARSTAR8);
			break;
		case 0xe:
		case 0x10:
			if(debug)printf("found BANK/ALSOBANK data\n");
			return(BANK);
			break;
		case 0xf:
			if(debug)printf("found COMPOSITE data\n");
			return(COMPOSITE);
			break;
		case 0xd:
		case 0x20:
			if(debug)printf("found SEGMENT/ALSOSEGMENT data\n");
			return(SEGMENT);
			break;
		case 0xc:
			if(debug)printf("found TAGSEGMENT data\n");
			return(TAGSEGMENT);
			break;
		default:
           printf("Unexpected fragment type: %d (0x%x)",type,type);
			return(UNEXPECTED);
			break;
	}
}


void parse_DPMBank(unsigned int *buf, int bank_length) {
   if(debug) printf("parse_DPMBank with length %i\n",bank_length);
	int ptr = 0;
	int length,type,word,fragType,seed, padding=0;
	unsigned short tag;
	unsigned short num;
    unsigned int expected;
    bool valid;
    unsigned int word_errors = 0;
    bool debug_pred = false;
	while (ptr<bank_length) {
		if (debug) printf("ptr = %d, bank_length = %d\n",ptr,bank_length);
		length      = buf[ptr]+1;
		tag         = (buf[ptr+1]>>16)&0xffff;
		type        = (buf[ptr+1]>>8)&0x3f;
		padding     = (buf[ptr+1]>>14)&0x3;
		num         = buf[ptr+1]&0xff;
        if(debug)
        {
           printf("Spitting out DPM bank header data: \n"); 
           printf("length: %i, tag: %d\n",length,tag);
           printf("padding: %i, type: %i, num: %d\n", padding, type, num);
        }
		fragType = getFragType(type);        
		if (fragType==UINT32)
		{
           seed = 5;
           expected = buf[seed];
           if(debug && debug_pred) printf("seed: 0x%x\n",expected);
           for(word=2; word<length; ++word) {
              if(debug && !debug_pred) {
                 printf("0x%x\t",buf[word]);
                 if((word-2)%5==0) printf("\n");
              }
              if(word<seed) {
                 //do nothing
              }
              else if(word==seed) {
                 // set the seed
                 expected = buf[word];            
              } else {
                 //predict
                 expected = flfsr32(expected);
                 if(debug && debug_pred) printf(" data: 0x%x pred: 0x%x -> ", buf[word], expected);
                 if(expected!=buf[word]) {
                    valid = false;
                    word_errors++;
                    if(debug && debug_pred) printf(" NO\n");
                 } else {
                    valid = true;
                    if(debug && debug_pred) printf(" YES\n");
                 }
              }
              if(debug && debug_pred) printf("\n");
           }
        }
		else {
           printf("data type of DPM bank should be UINT32 but was %d??\n",type);
           exit(1);
        }
		ptr+=length;
	}
    printf("\nword errors in this bank: %d\n",word_errors);
    
}


void check_lfsr(unsigned int *buf, int length) {
   if(debug) printf("parse_DPM with length %i\n",length);	   
   //uint *data_  = (uint *)malloc((length) * sizeof(uint));
   //memcpy(data_,buf,length*sizeof(uint));
   unsigned int expected;
   int word;
   bool valid;
   unsigned int word_errors = 0;
   // seed 
   expected = buf[0];
   for(word=1; word<length; ++word) {
      if(debug) printf("input 0x%x\t->\t0x%x?\n",expected,buf[word]);
      expected =  flfsr32(expected);
      if(debug) printf("predicted\t0x%x",expected);
      if(expected!=buf[word]) {
         valid = false;
         word_errors++;
         if(debug) printf("\tNO\n"); 
      } else {
         valid = true;
         if(debug) printf("\tYES\n"); 
      }
      //if(debug) {
      //   if(word%5!=0) printf("0x%x\t",buf[word]);
      //   else         printf("0x%x\n",buf[word]);
      //}
   }
   //free(data_);

}


unsigned int flfsr32(unsigned int input) {
/*  SsiPrbsTx.vhd
 -- PRBS Config
      PRBS_SEED_SIZE_G           : natural range 32 to 128    := 32;
      PRBS_TAPS_G                : NaturalArray               := (0 => 31, 1 => 6, 2 => 2, 3 => 1);
*/
   //unsigned int bit = ( (input >> 31) ^ (input >> 6) ^ (input >> 2) ^ (input >> 1) ) & 1;
   //unsigned int next = (input >> 1) | (bit << 31);
   unsigned int bit = ( (input >> 31) ^ (input >> 6) ^ (input >> 2) ^ (input >> 1 ) ) & 1; // &1 at the end picks out least sign bit
   unsigned int next = (input << 1) | bit; // shift into least sign bit
   return next;
}

unsigned short int lfsr16(unsigned short int input) {
   unsigned short int lfsr = input;
   lfsr = (lfsr >> 1) ^ (-(lfsr & 1u) & 0xB400u);
   return lfsr;
}

unsigned short int flfsr16(unsigned short int input) {
   unsigned short int bit = ( (input >> 0) ^ (input >> 2) ^ (input >> 3) ^ (input >> 5) ) & 1;
   return (input >> 1) | (bit << 15);
}


void parse_eventBank(unsigned int *buf, int bank_length) {
   if(debug) printf("parse_eventBank\n");
	int ptr = 0;
	int length,type, padding=0;
	unsigned short tag;
	unsigned short num;
	while (ptr<bank_length) {
       if (debug) printf("ptr = %d, bank_length = %d\n",ptr,bank_length);
       length      = buf[ptr]+1;
       tag         = (buf[ptr+1]>>16)&0xffff;
       type        = (buf[ptr+1]>>8)&0x3f;
       padding     = (buf[ptr+1]>>14)&0x3;
       num         = buf[ptr+1]&0xff;
       if(debug)
       {
          printf("Spitting out event bank header data: \n"); 
          printf("length: %i, tag: %d\n",length,tag);
          printf("padding: %i, type: %i, num: %d\n", padding, type, num);
       }
       int fragType = getFragType(type);
       if (fragType==UINT32)
       { 
          switch (tag) {
             case DPM0:
                if(debug) printf("DPM0 bank %d (0x%x)\n",tag,tag);
                parse_DPMBank(&buf[ptr+2], length-2);
                break;
             case DPM1:
                if(debug) printf("DPM1 bank %d (0x%x)\n",tag,tag);
                parse_DPMBank(&buf[ptr+2], length-2);
                break;
             case DPM2:
                if(debug) printf("DPM2 bank %d (0x%x)\n",tag,tag);
                parse_DPMBank(&buf[ptr+2], length-2);
                break;
             case DPM3:
                if(debug) printf("DPM3 bank %d (0x%x)\n",tag,tag);
                parse_DPMBank(&buf[ptr+2], length-2);
                break;
             case DPM4:
                if(debug) printf("DPM4 bank %d (0x%x)\n",tag,tag);
                parse_DPMBank(&buf[ptr+2], length-2);
                break;
             case DPM5:
                if(debug) printf("DPM5 bank %d (0x%x)\n",tag,tag);
                parse_DPMBank(&buf[ptr+2], length-2);
                break;
             case DPM6:
                if(debug) printf("DPM6 bank %d (0x%x)\n",tag,tag);
                parse_DPMBank(&buf[ptr+2], length-2);
                break;
             case DPM7:
                if(debug) printf("DPM7 bank %d (0x%x)\n",tag,tag);
                parse_DPMBank(&buf[ptr+2], length-2);
                break;
             default:
                if(debug) printf("Unexpected bank tag %d (0x%x)\n",tag,tag);
                break;
          }
       }
       else
          if (debug) printf("data type of event bank should be UINT32 but was %d (0x%x)\n",type,type);
       ptr+=length;
	}
}
void parse_event(unsigned int *buf) {
   int length,type, padding=0;
   unsigned short tag;
   unsigned short num;
   length      = buf[0]+1;
   tag         = (buf[1]>>16)&0xffff;
   padding     = (buf[1]>>14)&0x3;
   type        = (buf[1]>>8)&0x3f;
   num         = buf[1]&0xff;
   
   if(debug)
   {
      printf("Spitting out event header data: \nlength: %i, tag: %d\npadding: %i, type: %i, num: %d\n",length,tag,padding,type,num);
   }
   int fragType = getFragType(type);
   if (fragType!=BANK) printf("data type of event should be BANK but was %d\n",type);
   else {
      if(debug) printf("this is a BANK type\n");
   }
   
   parse_eventBank(&buf[fragment_offset[fragType]],length-fragment_offset[fragType]);
}
