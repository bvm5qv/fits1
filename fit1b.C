#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TROOT.h"


// entries is the number of random samples filled into the histogram
void fit1b(int experiments=1000, int entries=10, bool save=false) {
   //Simple histogram fitting examples
  gROOT->Reset();  // useful to reset ROOT to a cleaner state

  TFile *tf=0;
  if (save) tf=new TFile("histo.root","recreate");

  TH1F *randomHist1 = new TH1F("randomHist1", "Random Histogram1;x;frequency", 100, 0, 100);
  TH1F *randomHist2 = new TH1F("randomHist2", "Random Histogram2;x;frequency", 100, 0, 100);
  TRandom2 *generator=new TRandom2(0);  // parameter == seed, 0->use clock

  TH1F *Chi2Mean = new TH1F("Chi2Mean", "Chi2 Distribution Mean;x;frequency", 100, 48, 52);
  TH1F *NLLMean = new TH1F("NLLMean", "NLL Distribution Mean;x;frequency", 100, 48, 52);
  TH1F *Chi2Sigma = new TH1F("Chi2Sigma", "Chi2 Distribution Sigma;Val;frequency", 100, 8, 12);
  TH1F *NLLSigma = new TH1F("NLLSigma", "NLL Distribution Sigma;Val;frequency", 100, 8, 12);

  for(int j = 0; j < experiments; j++)
  {
    for (int i=0 ; i<entries ; i++){
      randomHist1->Fill(generator->Gaus(50,10)); // params: mean, sigma
    }
    gStyle->SetOptFit(1111); // show reduced chi2, probability, and params

    randomHist1->Fit("gaus");
    randomHist1->DrawCopy("e");
    TF1 *fitfunc1 = randomHist1->GetFunction("gaus");
    Chi2Mean->Fill(fitfunc1->GetParameter(1));
    Chi2Sigma->Fill(fitfunc1->GetParameter(2));

    randomHist1->Fit("gaus", "L");
	randomHist1->DrawCopy("e");
    TF1 *fitfunc2 = randomHist1->GetFunction("gaus");
    NLLMean->Fill(fitfunc2->GetParameter(1));
    NLLSigma->Fill(fitfunc2->GetParameter(2));
  }

  TCanvas *c1 = new TCanvas("c1", "Mean Plots", 1000, 800);
  c1->Divide(2,1);
  c1->cd(1);
  Chi2Mean->Draw();
  c1->cd(2);
  NLLMean->Draw();
  c1->Update();

  TCanvas *c2 = new TCanvas("c2", "Sigma Plots", 1000, 800);
  c2->Divide(2,1);
  c2->cd(1);
  Chi2Sigma->Draw();
  c2->cd(2);
  NLLSigma->Draw();
  c2->Update();

  if (save) {
    tf->Write();
    tf->Close();
  }
  cout << "Use .q to exit root" << endl;
}







