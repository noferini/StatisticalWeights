#include "tools.C"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

void testweight(){
  float par[2];
  StatisticalWeightsGaus manager;
  par[0] = 0; par[1] = 1;
  manager.AddSpecies("uno",par);
  par[0] = 2; par[1] = 1;
  manager.AddSpecies("due",par);
  par[0] = -1; par[1] = 1;
  manager.AddSpecies("tre",par);
  manager.Init();

  StatisticalWeights manager2;
  manager2.AddSpecies("unobis","TMath::Gaus(x,0,1)");
  manager2.AddSpecies("duebis","TMath::Gaus(x,2,1)");
  manager2.AddSpecies("trebis","TMath::Gaus(x,-1,1)");
  manager2.Init();

  Float_t xmin = manager.GetXmin();
  Float_t xmax = manager.GetXmax();
  Float_t xbin = manager.GetXbin();
  int nbin = (xmax-xmin)/xbin;

  TH1F *h1 = new TH1F("h1","",nbin,xmin,xmax);
  TH1F *h2 = new TH1F("h2","",nbin,xmin,xmax);
  TH1F *h3 = new TH1F("h3","",nbin,xmin,xmax);
  h2->SetLineColor(2);
  h3->SetLineColor(4);

  for(int i=1;i <= nbin;i++){
    h1->SetBinContent(i,manager.GetWeight(0,h1->GetBinCenter(i)));
    h2->SetBinContent(i,manager.GetWeight(1,h2->GetBinCenter(i)));
    h3->SetBinContent(i,manager.GetWeight(2,h3->GetBinCenter(i)));
  }

  new TCanvas;
  h1->Draw();
  h2->Draw("SAME");
  h3->Draw("SAME");
  manager.GetResponse(0)->Draw("SAME");
  manager.GetResponse(1)->Draw("SAME");
  manager.GetResponse(2)->Draw("SAME");


  // fast simulation
  Int_t nparticle[3];
  nparticle[0] = 100000;
  nparticle[1] = 30000;
  nparticle[2] = 50000;

  float nreco[3];
  nreco[0] = 0;
  nreco[1] = 0;
  nreco[2] = 0;
  float nreco2[3];
  nreco2[0] = 0;
  nreco2[1] = 0;
  nreco2[2] = 0;

  TH1F *h1or = new TH1F("h1or","",nbin,xmin,xmax);
  TH1F *h2or = new TH1F("h2or","",nbin,xmin,xmax);
  TH1F *h3or = new TH1F("h3or","",nbin,xmin,xmax);
  TH1F *htot = new TH1F("htot","",nbin,xmin,xmax);
  h1or->SetLineColor(2);
  h2or->SetLineColor(4);
  h3or->SetLineColor(6);

  for(int isp=0; isp < manager.GetNspecies(); isp++){
    manager.GetResponse(isp)->SetNpx(10000);
    for(int j=0; j < nparticle[isp]; j++){
      float x = manager.GetResponse(isp)->GetRandom(xmin,xmax);
      nreco[0] += manager.GetWeight(0,x);
      nreco[1] += manager.GetWeight(1,x);
      nreco[2] += manager.GetWeight(2,x);
      nreco2[0] += manager2.GetWeight(0,x);
      nreco2[1] += manager2.GetWeight(1,x);
      nreco2[2] += manager2.GetWeight(2,x);
      if(isp==0) h1or->Fill(x);
      if(isp==1) h2or->Fill(x);
      if(isp==2) h3or->Fill(x);
      htot->Fill(x);

    }
  }

  for(int isp=0; isp < manager.GetNspecies(); isp++){
    printf("%d) genareted = %d - reconstructed = %f (%f)\n",isp,nparticle[isp],nreco[isp],nreco2[isp]);
  }

  new TCanvas();
  htot->Draw();
  h1or->Draw("SAME");
  h2or->Draw("SAME");
  h3or->Draw("SAME");

  TF1 *f1 = new TF1("f1","gaus",xmin,xmax);
  TF1 *f2 = new TF1("f2","gaus",xmin,xmax);
  TF1 *f3 = new TF1("f3","gaus",xmin,xmax);

  f1->SetParameter(0,1);
  f1->SetParameter(1,0);
  f1->SetParameter(2,1);
  f2->SetParameter(0,1);
  f2->SetParameter(1,2);
  f2->SetParameter(2,1);
  f3->SetParameter(0,1);
  f3->SetParameter(1,-1);
  f3->SetParameter(2,1);

  f1->SetParameter(0,nreco[0]/f1->Integral(-10,10)*h1or->GetBinWidth(1));
  f2->SetParameter(0,nreco[1]/f2->Integral(-10,10)*h2or->GetBinWidth(1));
  f3->SetParameter(0,nreco[2]/f3->Integral(-10,10)*h3or->GetBinWidth(1));

  f1->Draw("SAME");
  f2->Draw("SAME");
  f3->Draw("SAME");
  f1->SetLineColor(1);
  f2->SetLineColor(1);
  f3->SetLineColor(1);

}
