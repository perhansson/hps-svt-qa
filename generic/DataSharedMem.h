//-----------------------------------------------------------------------------
// File          : DataSharedMem.h
// Author        : Ryan Herbst  <rherbst@slac.stanford.edu>
// Created       : 03/29/2013
// Project       : General Purpose
//-----------------------------------------------------------------------------
// Description :
// Shared memory for live display
//-----------------------------------------------------------------------------
// Copyright (c) 2013 by SLAC. All rights reserved.
// Proprietary and confidential to SLAC.
//-----------------------------------------------------------------------------
// Modification history :
// 01/11/2013: created
//-----------------------------------------------------------------------------
#ifndef __DATA_SHARED_MEM_H__
#define __DATA_SHARED_MEM_H__
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>   
#include <string.h>   
#include <unistd.h>
#include <stdio.h>

#define DATA_BUFF_SIZE 2097152
#define DATA_BUFF_CNT  20
#define DATA_NAME_SIZE 200

typedef struct {

   unsigned int wrAddr;
   unsigned int wrCount;
   unsigned int wrAddrLast;
   unsigned int wrCountLast;
   char         sharedName[DATA_NAME_SIZE];
   char         buffer[DATA_BUFF_CNT][DATA_BUFF_SIZE];

} DataSharedMemory;

// Open and map shared memory
inline int dataSharedOpenAndMap ( DataSharedMemory **ptr, const char *system, unsigned int id, int uid=-1 ) {
   int           smemFd;
   char          shmName[200];
   int           lid;

   // ID to use?
   if ( uid == -1 ) lid = getuid();
   else lid = uid;

   // Generate shared memory
   sprintf(shmName,"data_shared.%i.%s.%i",lid,system,id);

   // Attempt to open existing shared memory
   if ( (smemFd = shm_open(shmName, O_RDWR, (S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH)) ) < 0 ) {

      // Otherwise open and create shared memory
      if ( (smemFd = shm_open(shmName, (O_CREAT | O_RDWR), (S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH)) ) < 0 ) return(-1);

      // Force permissions regardless of umask
      fchmod(smemFd, (S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH));
    
      // Set the size of the shared memory segment
      ftruncate(smemFd, sizeof(DataSharedMemory));
   }

   // Map the shared memory
   if((*ptr = (DataSharedMemory *)mmap(0, sizeof(DataSharedMemory),
              (PROT_READ | PROT_WRITE), MAP_SHARED, smemFd, 0)) == MAP_FAILED) return(-2);

   // Store name
   strcpy((*ptr)->sharedName,shmName);

   return(smemFd);
}

// Close shared memory
inline void dataSharedClose ( DataSharedMemory *ptr ) {
   char shmName[200];

   // Get shared name
   strcpy(shmName,ptr->sharedName);

   // Unlink
   shm_unlink(shmName);
}

// Init data structure, called by ControlServer
inline void dataSharedInit ( DataSharedMemory *ptr ) {
   ptr->wrAddr      = 0;
   ptr->wrAddrLast  = 0;
   ptr->wrCount     = 0;
   ptr->wrCountLast = 0;
}

// Write to shared buffer
inline void dataSharedWrite ( DataSharedMemory *ptr, unsigned int flag, const char *data, unsigned int count ) {
   if ( (count+1) < DATA_BUFF_SIZE ) {
      ptr->wrAddrLast  = ptr->wrAddr;
      ptr->wrCountLast = ptr->wrCount;

      *((unsigned int *)ptr->buffer[ptr->wrAddr]) = flag;

      memcpy(&(ptr->buffer[ptr->wrAddr][4]),data,count);

      ptr->wrCount++;
      ptr->wrAddr++;
      if ( ptr->wrAddr == DATA_BUFF_CNT ) ptr->wrAddr = 0;
   }
}

// Read from shared buffer
inline int dataSharedRead ( DataSharedMemory *ptr, unsigned int *rdAddr, unsigned int *rdCount, unsigned int *flag, char **data ) {

   // Detect if reader is ahead of writer or too far behind writer
   if ( (*rdCount == 0 && ptr->wrCount != 0) || ( *rdCount > ptr->wrCount ) || ( (ptr->wrCount - *rdCount) >= (DATA_BUFF_CNT/2) ) ) {
      printf("Adjusting pointers. wrAddr=%i, wrCount=%i, wrAddrLast=%i, wrCountLast=%i, rdAddr=%i, rdCount=%i\n",ptr->wrAddr,ptr->wrCount,
         ptr->wrAddrLast,ptr->wrCountLast,*rdAddr,*rdCount);
      *rdAddr  = ptr->wrAddrLast;
      *rdCount = ptr->wrCountLast;
   }

   if ( *rdCount == ptr->wrCount ) return(0);

   *flag = *((unsigned int *)ptr->buffer[*rdAddr]);
   *data = &(ptr->buffer[*rdAddr][4]);

   (*rdCount)++;
   (*rdAddr)++;
   if ( *rdAddr == DATA_BUFF_CNT ) *rdAddr = 0;
   return(1);
}

#endif

