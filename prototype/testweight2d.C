#include "tools.C"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF2.h"

void testweight2d(){
  StatisticalWeights *manager = new StatisticalWeights();
  manager->AddSpecies("uno","TMath::Gaus(x,-1,1)");
  manager->AddSpecies("due","TMath::Gaus(x,1,1)");

  manager->AddSignal(0,"TMath::Gaus(x,-0.5,1)");
  manager->AddSignal(1,"TMath::Gaus(x,0.5,1)");
  manager->Init();

  float par[2];
  StatisticalWeightsGaus managerGaus;
  par[0] = -1; par[1] = 1;
  managerGaus.AddSpecies("unoG",par);
  par[0] = -0.5; par[1] = 1;
  managerGaus.AddSignal(0,par);

  par[0] = 1; par[1] = 1;
  managerGaus.AddSpecies("dueG",par);
  par[0] = 0.5; par[1] = 1;
  managerGaus.AddSignal(1,par);

  managerGaus.Init();

  manager = &managerGaus;

  Float_t xmin = manager->GetXmin();
  Float_t xmax = manager->GetXmax();
  Float_t xbin = manager->GetXbin();
  int nbin = (xmax-xmin)/xbin;

  TH2F *h1 = new TH2F("h1","",nbin,xmin,xmax,nbin,xmin,xmax);
  TH2F *h2 = new TH2F("h2","",nbin,xmin,xmax,nbin,xmin,xmax);
  h1->SetLineColor(2);
  h2->SetLineColor(4);

  float x[2];
  for(int i=1;i <= nbin;i++){
    for(int j=1;j <= nbin;j++){
      x[0] = h1->GetXaxis()->GetBinCenter(i);
      x[1] = h1->GetYaxis()->GetBinCenter(j);
      h1->SetBinContent(i,j,manager->GetWeight(0,x));
      h2->SetBinContent(i,j,manager->GetWeight(1,x));
    }
  }

  TCanvas *c1 = new TCanvas;
  c1->Divide(2,1);
  c1->cd(1);
  h1->Draw("lego");
  c1->cd(2);
  h2->Draw("lego");

  TCanvas *c2 = new TCanvas;
  c2->Divide(2,1);
  c2->cd(1);
  manager->GetResponse(0,0)->Draw("");
  manager->GetResponse(1,0)->Draw("SAME");
  manager->GetResponse(0,0)->SetLineColor(2);
  manager->GetResponse(1,0)->SetLineColor(4);
  c2->cd(2);
  manager->GetResponse(0,1)->Draw("");
  manager->GetResponse(1,1)->Draw("SAME");
  manager->GetResponse(0,1)->SetLineColor(2);
  manager->GetResponse(1,1)->SetLineColor(4);


  // fast simulation
  Int_t nparticle[2];
   nparticle[0] = 100000;
   nparticle[1] = 30000;

   float nreco[2];
   nreco[0] = 0;
   nreco[1] = 0;

   TH2F *h1or = new TH2F("h1or","",nbin,xmin,xmax,nbin,xmin,xmax);
   TH2F *h2or = new TH2F("h2or","",nbin,xmin,xmax,nbin,xmin,xmax);
   TH2F *htot = new TH2F("htot","",nbin,xmin,xmax,nbin,xmin,xmax);
   h1or->SetLineColor(2);
   h2or->SetLineColor(4);

   for(int isp=0; isp < manager->GetNspecies(); isp++){
     for(int j=0; j < nparticle[isp]; j++){
       manager->SetPar(0,1,-1);
       x[0] = manager->GetResponse(isp,0)->GetRandom(xmin,xmax);
       x[1] = manager->GetResponse(isp,1)->GetRandom(xmin,xmax);
       nreco[0] += manager->GetWeight(0,x);
       nreco[1] += manager->GetWeight(1,x);

       if(isp==0) h1or->Fill(x[0],x[1]);
       if(isp==1) h2or->Fill(x[0],x[1]);
       htot->Fill(x[0],x[1]);

     }
   }

   for(int isp=0; isp < manager->GetNspecies(); isp++){
     printf("%d) genareted = %d - reconstructed = %f\n",isp,nparticle[isp],nreco[isp]);
   }

   TCanvas *c3 = new TCanvas();
   c3->Divide(3,1);
   c3->cd(1);
   htot->Draw("lego");
   c3->cd(2);
   h1or->Draw("lego");
   c3->cd(3);
   h2or->Draw("lego");

   h1or->SetMaximum(700);
   h2or->SetMaximum(700);
   htot->SetMaximum(700);

   TF2 *f1 = new TF2("f1","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",xmin,xmax,xmin,xmax);
   TF2 *f2 = new TF2("f2","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",xmin,xmax,xmin,xmax);

   f1->SetParameter(0,1);
   f1->SetParameter(1,-1);
   f1->SetParameter(2,1);
   f1->SetParameter(3,-0.5);
   f1->SetParameter(4,1);
   f2->SetParameter(0,1);
   f2->SetParameter(1,1);
   f2->SetParameter(2,1);
   f2->SetParameter(3,0.5);
   f2->SetParameter(4,1);

   f1->SetParameter(0,nreco[0]/f1->Integral(-10,10,-10,10)*h1or->GetXaxis()->GetBinWidth(1)*h1or->GetYaxis()->GetBinWidth(1));
   f2->SetParameter(0,nreco[1]/f2->Integral(-10,10,-10,10)*h2or->GetXaxis()->GetBinWidth(1)*h2or->GetYaxis()->GetBinWidth(1));

   TCanvas *c4 = new TCanvas();
   c4->Divide(2,1);
   c4->cd(1);
   f1->Draw("SURF");
   f1->SetMaximum(700);
   c4->cd(2);
   f2->Draw("SURF");
   f2->SetMaximum(700);
   f1->SetLineColor(2);
   f2->SetLineColor(4);
}
