#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TROOT.h"

#include <iostream>
#include <algorithm>
#include <vector>

void histofit1a()
{
  TFile *f = TFile::Open("histo25.root");
  TH1F *h = (TH1F*)f->Get("randomHist1");
  
  TF1 *fitFcn = new TF1("fitFcn", "gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  h->Fit(fitFcn, "L");  // Likelihood fit
  double best_mean = fitFcn->GetParameter(1);
  double best_sigma = fitFcn->GetParameter(2);
  double best_ampl = fitFcn->GetParameter(0);
  
  std::vector<double> means, nll_values;
  for (double m = best_mean - 3*best_sigma; m <= best_mean + 3*best_sigma; m += 0.1*best_sigma) {
      double NLL = 0;
      for (int i = 1; i <= h->GetNbinsX(); ++i) {
          double n = h->GetBinContent(i);
          double x = h->GetBinCenter(i);
          double mu = best_ampl * exp(-0.5 * pow((x - m)/best_sigma, 2)) * h->GetBinWidth(i);
          if (mu > 0 && n > 0) NLL += -2 * (n * log(mu) - mu);
      }
      means.push_back(m);
      nll_values.push_back(NLL);
  }
  TGraph *gNLL = new TGraph(means.size(), &means[0], &nll_values[0]);
  gNLL->SetTitle("-2ln(L) vs mean; Mean; -2ln(L)");
  gNLL->SetLineColor(kBlue);
  gNLL->Draw("AL");

    
  TFile *f2 = TFile::Open("histo1k.root");
  TH1F *h2 = (TH1F*)f2->Get("randomHist1");
  
  TF1 *fitChi2 = new TF1("fitChi2", "gaus", h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax());
  h2->Fit(fitChi2);   // Default = Chi2 fit
  double mean_chi2 = fitChi2->GetParameter(1);
  double sigma_chi2 = fitChi2->GetParameter(2);
  double ampl_chi2 = fitChi2->GetParameter(0);
  
  std::vector<double> means2, chi2_values;
  for (double m = mean_chi2 - 3*sigma_chi2; m <= mean_chi2 + 3*sigma_chi2; m += 0.1*sigma_chi2) {
      double chi2 = 0;
      for (int i = 1; i <= h2->GetNbinsX(); ++i) {
          double n = h2->GetBinContent(i);
          double err = h2->GetBinError(i);
          double x = h2->GetBinCenter(i);
          double mu = ampl_chi2 * exp(-0.5 * pow((x - m)/sigma_chi2, 2));
          if (err > 0) chi2 += pow((n - mu)/err, 2);
      }
      means2.push_back(m);
      chi2_values.push_back(chi2);
  }
  TGraph *gChi2 = new TGraph(means2.size(), &means2[0], &chi2_values[0]);
  gChi2->SetTitle("Chi2 vs mean; Mean; Chi2");
  gChi2->SetLineColor(kGreen+2);
  gChi2->Draw("AL");

  TCanvas *c1 = new TCanvas("c1", "NLL and Chi2 Comparison", 800, 600);
  gNLL->SetLineColor(kBlue);
  gNLL->SetLineWidth(2);
  gNLL->SetTitle("Comparison of -2ln(L) and Chi2 Contours;Mean;Value");
  gNLL->Draw("AL");  // "A" = draw axes, "L" = draw line

  gChi2->SetLineColor(kRed);
  gChi2->SetLineWidth(2);
  gChi2->Draw("L SAME");

  double nll_min = *std::min_element(nll_values.begin(), nll_values.end());
  double chi2_min = *std::min_element(chi2_values.begin(), chi2_values.end()); 
  TLine *l_nll1 = new TLine(gNLL->GetXaxis()->GetXmin(), nll_min + 1, gNLL->GetXaxis()->GetXmax(), nll_min + 1);
  l_nll1->SetLineStyle(2);
  l_nll1->SetLineColor(kBlue);
  l_nll1->Draw("same");
  TLine *l_chi1 = new TLine(gChi2->GetXaxis()->GetXmin(), chi2_min + 1, gChi2->GetXaxis()->GetXmax(), chi2_min + 1);
  l_chi1->SetLineStyle(2);
  l_chi1->SetLineColor(kRed);
  l_chi1->Draw("same");

  auto legend = new TLegend(0.65, 0.75, 0.95, 0.88);
  legend->AddEntry(gNLL, "-2ln(L) (histo25.root)", "l");
  legend->AddEntry(gChi2, "Chi2 (histo1k.root)", "l");
  legend->Draw();
    
  std::cout << "NLL par error: " << fitFcn->GetParError(1) << '\n';
  std::cout << "Chi2 par error: " << fitChi2->GetParError(1) << '\n';

}


