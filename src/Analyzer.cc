//AUTHOR: Tomasz Rudnicki

#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TF1.h>
#include <CFD.cc>

const Int_t nbins = 1024;
const Int_t nchns = 2;
const Double_t overflow = -0.499;
const Double_t frac = 0.2;

TGraph* template_;
double fit(double* x, double* pars);

int main(int argc, char** argv){

  //------- Usage check ---------------------------------------
  if(argc < 2){
    std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
    return 1;
  }

  //------- Open file and load tree with data -------------------
  TFile *f = new TFile(argv[1]);
  TTree *t1 = (TTree*)f->Get("Waveforms");

  //----- Set branches ------------------------------------------
  Double_t waveform[nchns][nbins], time[nchns][nbins];
  t1->SetBranchAddress("waveform1", &waveform[0]);
  t1->SetBranchAddress("waveform2", &waveform[1]);
  t1->SetBranchAddress("time1", &time[0]);
  t1->SetBranchAddress("time2", &time[1]);

  //------- Open root files for results ------------------
  TFile file1("../root_files/TemplateResults.root", "RECREATE");

  //------ Create template shape of waveform -------------
  std::cout << "\n Shape Creation\n" << std::endl;

  Double_t Taxis[nbins], shape[nchns][nbins];
  for(int i=0; i<nbins; i++) {Taxis[i]=i*0.2; shape[0][i] = 0; shape[1][i] = 0;}

  Int_t nentries = t1->GetEntries(), counter1 = 0, counter2 = 0;
  for(int nentry=0; nentry<nentries; nentry++){
    t1->GetEntry(nentry);
    if (Create_temp(Taxis, shape[0], waveform[0], time[0], frac, nbins, overflow) == 0) counter1++;
    if (Create_temp(Taxis, shape[1], waveform[1], time[1], frac, nbins, overflow) == 0) counter2++;
  }
  for(int i=0; i<nbins; i++) {shape[0][i] = shape[0][i]/counter1; shape[1][i] = shape[1][i]/counter2;}

  //------ Prepare file for data saving ------------------
  std::cout << "\n Data Analysis\n" << std::endl;

  TTree *t2 = new TTree("Resolution", "Resolution");
  Double_t dt1_2;
  t2->Branch("dt1_2", &dt1_2, "dt1_2/D");

  //------ Loop over entries & channels ---------------------
  for(Int_t nentry=0; nentry<nentries; nentry++){
    t1->GetEntry(nentry);
    Double_t t0[nchns];
    Bool_t breaker = false;
    for(Int_t i=0; i<nchns; i++){

      //------- Waveform correctness checker --------------------
      Int_t current_max_ind = TMath::LocMin(nbins, waveform[i]);
      if (waveform[i][current_max_ind] <= overflow
         || waveform[i][current_max_ind] > -0.15
         || current_max_ind < 20 || current_max_ind > nbins-20) {breaker = true;}
      if(breaker) continue;

      //------ Fit template to impuls ---------------------------
      TGraph gr1 (nbins, time[i], waveform[i]);
      TGraph *temporary = new TGraph(nbins, Taxis, shape[i]);
      template_ = temporary;

      TF1 f("fit", fit, time[i][current_max_ind]-6, time[i][current_max_ind]+3, 2);
      f.SetParameter(0, time[i][current_max_ind]-50.);
      gr1.Fit(&f, "Q");
      t0[i] = f.GetParameter(0);
    }

    //------ Save result dt if waveform classified as correct ----------
    dt1_2 = t0[0]-t0[1];
    if(!breaker) t2->Fill();
  }

  //--- Clear and finish ------------------------------------
  t2->Write();
  delete t2;
  delete t1;
  file1.Close();
  f->Close();

  return 0;
}

double fit(double* x, double* pars) {

  const double xx = x[0]-pars[0];
  Double_t yy = template_->Eval(xx);
  const double y1 = pars[1]*yy;
  return y1;
}

