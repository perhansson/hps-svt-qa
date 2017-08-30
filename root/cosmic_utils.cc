#include "cosmic_utils.hh"

double cosmic_utils::getMaxSampleOfThree(double* samples) {

  // find maximum among three consecutive samples
  double peak_val = -1.;
  double y0,y1,y2;
  for ( int i=1; i<4; ++i ) {
    y0 = samples[i-1];
    y1 = samples[i];
    y2 = samples[i+1];
    if ( y1 > y0 && y1 > y2) {
      if ( y1 > peak_val ) {
	peak_val = y1;
      }
    }
  }
  return peak_val;
}

int cosmic_utils::getNrOfConsecutiveSamplesAboveThresh(double* samples,double thresh) {
  
  unsigned int c = 0;
  unsigned int cmax = 0;
  bool prevAbove = true;
  for ( unsigned int i=0; i<6; ++i ) {
    if ( samples[i] > thresh ) {
      if ( prevAbove ) {
	++c;
      }
      prevAbove = true;
    } else {
      if ( c > cmax) cmax = c;      
      prevAbove = false;
      c=0;
    }
  }
  return cmax;
}
