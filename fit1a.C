#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TROOT.h"


// fit1.C
// entries is the number of random samples filled into the histogram
void fit1a(int experiments=1000, int entries=1000, bool save=false) {
   //Simple histogram fitting examples
  gROOT->Reset();  // useful to reset ROOT to a cleaner state

  TFile *tf=0;
  if (save) tf=new TFile("histo.root","recreate");

  TH1F *randomHist1 = new TH1F("randomHist1", "Random Histogram;x;frequency", 100, 0, 100);
  TRandom2 *generator=new TRandom2(0);  // parameter == seed, 0->use clock

  TH1F *rChiSquare = new TH1F("rChiSquare", "Reduced Chi2;x;frequency", 100, 0, 2);
  TGraph *chiSquareProb = new TGraph();
  chiSquareProb->SetTitle("Chi2 p-value;Reduced Chi2;p-value");
  chiSquareProb->SetMarkerStyle(21);
  TH1F *distMean = new TH1F("distMean", "Distribution mean;x;frequency", 100, 49, 51);
  TH1F *distMeanError = new TH1F("distMeanError", "Error of Mean;x;frequency", 100, 0, 0.2);

  for(int j = 0; j < experiments; j++)
  {
    for (int i=0 ; i<entries ; i++){
      randomHist1->Fill(generator->Gaus(50,10)); // params: mean, sigma
    }
    // simple fits may be performed automatically
    // gStyle->SetOptFit(111);  // show reduced chi2 and params
    gStyle->SetOptFit(1111); // show reduced chi2, probability, and params
    randomHist1->Fit("gaus");  
    randomHist1->DrawCopy("e");

    TF1 *fitfunc = randomHist1->GetFunction("gaus");
    rChiSquare->Fill(fitfunc->GetChisquare() / fitfunc->GetNDF());
    chiSquareProb->SetPoint(j, fitfunc->GetChisquare() / fitfunc->GetNDF(), fitfunc->GetProb());
    distMean->Fill(fitfunc->GetParameter(1));
    distMeanError->Fill(fitfunc->GetParError(1));
  }

  TCanvas *c = new TCanvas("c", "Plots", 1000, 800);
  c->Divide(2,2);

  c->cd(1);
  rChiSquare->Draw();

  c->cd(2);
  distMean->Draw();

  c->cd(3);
  chiSquareProb->Draw("AP");

  c->cd(4);
  distMeanError->Draw();

  c->Update();

  if (save) {
    tf->Write();
    tf->Close();
  }
  cout << "Use .q to exit root" << endl;
}







