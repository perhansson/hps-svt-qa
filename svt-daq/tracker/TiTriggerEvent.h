#ifndef __TI_TRIGGER_EVENT_H__
#define __TI_TRIGGER_EVENT_H__

#include "TriggerEvent.h"


//! Trigger Event Container Class
class TiTriggerEvent : public TriggerEvent {

  
 public:
  
  TiTriggerEvent();
  ~TiTriggerEvent();
  
  // overrides
  uint count();
  void sample (uint index, TrackerSample* sample);

  unsigned long timeStamp();
  unsigned long tiEventNumber();
  
 private:

  static const uint _tiDataSize   = 4;
  
};


#endif
