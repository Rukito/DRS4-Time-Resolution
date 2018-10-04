#include <TGraph.h>
#include <iostream>
#include <TF1.h>
#include <TFitResult.h>
#include <TMath.h>

Int_t Create_temp(const Double_t *Taxis, Double_t *shape, const Double_t *waveform, const Double_t *time,
                  const Double_t frac, const Int_t nbins, const Double_t overflow) {

 //----- Find peak index ---------------------
  TGraph gr1(nbins, time, waveform);
  Int_t index_max = TMath::LocMin(nbins, waveform);

  //------ Check if signal is below threshold, overflow or not signal at all (no true peak) -------------
  if(index_max < 30 || index_max > nbins-30 || waveform[index_max]<=overflow || waveform[index_max]>-0.15) {return 1;}

  //------ Find signal maximum ----------------------------------------------
  TFitResultPtr r1 = gr1.Fit("pol2","QSF","",time[index_max-20], time[index_max+20]);
  Double_t p2 = r1->Parameter(2);
  Double_t p1 = r1->Parameter(1);
  Double_t p0 = r1->Parameter(0);
  Double_t max = - (p1*p1 - 4*p0*p2) / (4*p2);

  //------ Check if fitted extremum is overflow ------------------------------
  Double_t threshold = frac*max;
  if(max<=overflow || waveform[index_max]>threshold || threshold>=0) {return 1;}

  //-------- Find where peak == threshold ------------------------------------
  Int_t index_thr=0;
  for(Int_t ind=index_max; ind>10; ind--){
    if(waveform[ind]<=threshold && waveform[ind-1]>threshold) {index_thr = ind; break;}
  }
  if (index_thr == 0) {return 1;}

  //------ Signal shape added to the array ---------------------------------
  int i=0, j=0;
  Double_t time_beg = time[index_max] - 50.;
  if (time_beg < .0) j = int( abs(floor(time_beg/0.2)) );
  Double_t a=0, b=0;

  for (i; i<nbins; i++) {
    if (j >= nbins) {break;}
    while (Taxis[j] < (time[i]-time_beg) && j<nbins) {
      shape[j] += Taxis[j] * a + b;
      j++;
    }
    if ( ((time[i]-time_beg) <= Taxis[j] ) && ((time[i+1]-time_beg) > Taxis[j]) ) {
      a = (waveform[i+1] - waveform[i]) / (time[i+1] - time[i]);
      b = waveform[i] - a * (time[i]-time_beg);
      shape[j] += Taxis[j] * a + b;
      j++;
    }
  }

  return 0;
}
