#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TROOT.h"

#include <iostream>
#include <vector>

void histofit()
{
  TFile *f = TFile::Open("histo25.root");
  TH1F *h = (TH1F*)f->Get("randomHist1");
  h->Draw();

  TF1 *fitFcn = new TF1("fitFcn", "gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  h->Fit(fitFcn, "L");

double NLL_data = 0;
for (int i = 1; i <= h->GetNbinsX(); ++i) {
    double n = h->GetBinContent(i);
    double x = h->GetBinCenter(i);
    double mu = fitFcn->Eval(x) * h->GetBinWidth(i); // expected count
    if (mu > 0 && n > 0) NLL_data += -2 * (n * log(mu) - mu);
}
std::cout << "NLL for data = " << NLL_data << endl;

int Ntoys = 1000;
TH1F *hToy = (TH1F*)h->Clone("hToy");
hToy->Reset();

std::vector<double> NLL_toys;
for (int t = 0; t < Ntoys; ++t) {
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double mu = fitFcn->Eval(h->GetBinCenter(i)) * h->GetBinWidth(i);
        int n = gRandom->Poisson(mu);
        hToy->SetBinContent(i, n);
    }

    double NLL_toy = 0;
    for (int i = 1; i <= hToy->GetNbinsX(); ++i) {
        double n = hToy->GetBinContent(i);
        double mu = fitFcn->Eval(hToy->GetBinCenter(i)) * hToy->GetBinWidth(i);
        if (mu > 0 && n > 0) NLL_toy += -2 * (n * log(mu) - mu);
    }
    NLL_toys.push_back(NLL_toy);
}

TH1F *hNLL = new TH1F("hNLL", "NLL from pseudo-experiments; -2logL; Entries", 50, 
                       TMath::MinElement(NLL_toys.size(), &NLL_toys[0]) - 1,
                       TMath::MaxElement(NLL_toys.size(), &NLL_toys[0]) + 1);

for (auto val : NLL_toys) hNLL->Fill(val);

hNLL->Draw();
TLine *line = new TLine(NLL_data, 0, NLL_data, hNLL->GetMaximum());
line->SetLineColor(kRed);
line->SetLineWidth(2);
line->Draw("same");

int nAbove = 0;
for (auto val : NLL_toys) if (val > NLL_data) nAbove++;
double p_value = (double)nAbove / NLL_toys.size();
cout << "p-value = " << p_value << endl;


}


