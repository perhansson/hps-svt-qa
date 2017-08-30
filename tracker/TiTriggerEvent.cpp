#include "TiTriggerEvent.h"

TiTriggerEvent::TiTriggerEvent() : TriggerEvent() {}

TiTriggerEvent::~TiTriggerEvent() {}


// Get sample count
uint TiTriggerEvent::count ( ) {
   if (eventCodeMatch()) {
      return((size_- (kHeadSize + kTailSize + _tiDataSize) ) / sampleSize_);
   } else {
      return 0;
   }
}

// Get sample at index
void TiTriggerEvent::sample (uint index, TrackerSample* sample) {
   if ( index < count() ) {
      sample->setData(&(data_[kHeadSize + (index*sampleSize_)]));
   }
}

// Get TI data
unsigned long TiTriggerEvent::timeStamp() {
  uint w3 = data_[size_ - 1]; 
  uint w2 = data_[size_ - 2]; 
  //printf("w0 w1 w2 w3:  0x%x 0x%x 0x%x 0x%x\n",w0, w1, w2, w3);
  unsigned long t = 0;
  t = w2;
  unsigned long tu = w3 & 0x0000ffff;
  t = (tu << 32) | t; 
  return t;
}

// Get TI event number
unsigned long TiTriggerEvent::tiEventNumber() {
  uint w3 = data_[size_ - 1]; 
  uint w1 = data_[size_ - 3]; 
  unsigned long t = 0;
  t = w1;
  unsigned long tu = w3 & 0xffff0000;
  t = (tu << 32) | t; 
  return t;

}
