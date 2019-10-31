#include "tools.C"
void makeperf(){
  float par[2];
  StatisticalWeightsGaus manager;
  manager.SetRange(-1000,1000);

  par[0]=-100, par[1] = 20;
  manager.AddSpecies("pion",par);
  par[0]=0, par[1] = 20;
  manager.AddSpecies("kaon",par);
  par[0]=100, par[1] = 20;
  manager.AddSpecies("proton",par);
  par[0]=-200, par[1] = 20;
  manager.AddSpecies("electron",par);

  manager.Init();

  // load priors
  TH2F *hPriorEl,*hPriorPi,*hPriorKa,*hPriorPr;
  TFile *fPriors = TFile::Open("priors.root");
  if(! fPriors){
    printf("Set Default priors\n");
    hPriorEl = new TH2F("hPriorsEl","",20,0,10,20,-1.2,1.2);
    hPriorPi = new TH2F("hPriorsPi","",20,0,10,20,-1.2,1.2);
    hPriorKa = new TH2F("hPriorsKa","",20,0,10,20,-1.2,1.2);
    hPriorPr = new TH2F("hPriorsPr","",20,0,10,20,-1.2,1.2);
    for(int i=1;i <= 20;i++){
      for(int j=1;j <= 20;j++){
	hPriorEl->SetBinContent(i,j,1);
	hPriorPi->SetBinContent(i,j,1);
	hPriorKa->SetBinContent(i,j,1);
	hPriorPr->SetBinContent(i,j,1);
      }
    }
  }
  else{
    hPriorEl = (TH2F *) fPriors->Get("hPriorsEl");
    hPriorPi = (TH2F *) fPriors->Get("hPriorsPi");
    hPriorKa = (TH2F *) fPriors->Get("hPriorsKa");
    hPriorPr = (TH2F *) fPriors->Get("hPriorsPr");
  }
  hPriorEl->SetName("oldElpr");
  hPriorPi->SetName("oldPipr");
  hPriorKa->SetName("oldKapr");
  hPriorPr->SetName("oldPrpr");

  TH2F *hPtBayesEl,*hPtBayesPi,*hPtBayesKa,*hPtBayesPr;
  hPtBayesEl = new TH2F("hPriorsEl","",20,0,10,20,-1.2,1.2);
  hPtBayesPi = new TH2F("hPriorsPi","",20,0,10,20,-1.2,1.2);
  hPtBayesKa = new TH2F("hPriorsKa","",20,0,10,20,-1.2,1.2);
  hPtBayesPr = new TH2F("hPriorsPr","",20,0,10,20,-1.2,1.2);

  TH2F *hPtTrueEl,*hPtTruePi,*hPtTrueKa,*hPtTruePr;
  hPtTrueEl = new TH2F("hTrueEl","",20,0,10,20,-1.2,1.2);
  hPtTruePi = new TH2F("hTruePi","",20,0,10,20,-1.2,1.2);
  hPtTrueKa = new TH2F("hTrueKa","",20,0,10,20,-1.2,1.2);
  hPtTruePr = new TH2F("hTruePr","",20,0,10,20,-1.2,1.2);

  gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");

  int pdgEl = 11;  // electron ID
  int pdgPi = 211; // pion     ID
  int pdgKa = 321; // kaon     ID
  int pdgPr = 2212;// proton   ID

  float mEl = 0.000511; // electron mass in GeV/c^2
  float mPi = 0.139570; // pion     mass in GeV/c^2
  float mKa = 0.493600; // kaon     mass in GeV/c^2
  float mPr = 0.938270; // proton   mass in GeV/c^2

  TFile *f = new TFile("my_output.root");
  TTree *t = (TTree *) f->Get("tree");

  int nev = t->GetEntries();

  TH1F *hbtof = new TH1F("hbtof",";#beta",101,0,1.01);
  TH2F *htof2D = new TH2F("htof2D",";p;#beta",100,0,10,150,0,1.5);
  TH2F *hdeltat2Dpi = new TH2F("hdeltat2Dpi",";p;t - t_{exp}^{K} (ps)",100,0,10,200,-1000,1000);
  TH2F *hdeltat2Dka = new TH2F("hdeltat2Dka",";p;t - t_{exp}^{K} (ps)",100,0,10,200,-1000,1000);
  TH2F *hdeltat2Dpr = new TH2F("hdeltat2Dpr",";p;t - t_{exp}^{K} (ps)",100,0,10,200,-1000,1000);

  TH2F *hmaskCorr = new TH2F("hmaskCorr",";det1;det2",5,0,5,5,0,5);
  hmaskCorr->Fill("btof","btof",0);
  hmaskCorr->Fill("etof","etof",0);
  hmaskCorr->Fill("htof1","htof1",0);
  hmaskCorr->Fill("htof2","htof2",0);
  hmaskCorr->Fill("htof2p","htof2p",0);


  // set branch addresses
  // particles
  std::vector<TParticle> part;
  std::vector<TParticle> *apart = &part;  
  t->SetBranchAddress("particle",&apart);
  // mask with tracking flags
  std::vector<int> mask;
  std::vector<int> *amask = &mask;  
  t->SetBranchAddress("trackMask",&amask);
  // btof infos
  std::vector<float> btof_tof;
  std::vector<float> *abtof_tof = &btof_tof;  
  t->SetBranchAddress("btof_tof",&abtof_tof);
  std::vector<float> btof_texp;
  std::vector<float> *abtof_texp = &btof_texp;  
  t->SetBranchAddress("btof_texp",&abtof_texp);
  // etof infos
  std::vector<float> etof_tof;
  std::vector<float> *aetof_tof = &etof_tof;  
  t->SetBranchAddress("etof_tof",&aetof_tof);
  std::vector<float> etof_texp;
  std::vector<float> *aetof_texp = &etof_texp;  
  t->SetBranchAddress("etof_texp",&aetof_texp);
  // htof1 infos
  std::vector<float> htof1_tof;
  std::vector<float> *ahtof1_tof = &htof1_tof;  
  t->SetBranchAddress("htof1_tof",&ahtof1_tof);
  std::vector<float> htof1_texp;
  std::vector<float> *ahtof1_texp = &htof1_texp;  
  t->SetBranchAddress("htof1_texp",&ahtof1_texp);
  // htof2 infos
  std::vector<float> htof2_tof;
  std::vector<float> *ahtof2_tof = &htof2_tof;  
  t->SetBranchAddress("htof2_tof",&ahtof2_tof);
  std::vector<float> htof2_texp;
  std::vector<float> *ahtof2_texp = &htof2_texp;  
  t->SetBranchAddress("htof2_texp",&ahtof2_texp);
  // htof2p infos
  std::vector<float> htof2p_tof;
  std::vector<float> *ahtof2p_tof = &htof2p_tof;  
  t->SetBranchAddress("htof2p_tof",&ahtof2p_tof);
  std::vector<float> htof2p_texp;
  std::vector<float> *ahtof2p_texp = &htof2p_texp;  
  t->SetBranchAddress("htof2p_texp",&ahtof2p_texp);

  float p,pt,betaEl,betaPi,betaKa,betaPr,meanEl,meanPi,meanPr,weightPi,weightKa,weightPr,weightEl,time,expEl,expPi,expKa,expPr,distVert,distVertT,eta,phi;
  int particleID;
  float bayesEl,bayesPi,bayesKa,bayesPr;

  TFile *fout =  new TFile("results.root","RECREATE");
  TTree *tOut = new TTree("tOut","tOut");
  tOut->Branch("p",&p,"p/F");
  tOut->Branch("pt",&pt,"pt/F");
  tOut->Branch("eta",&eta,"eta/F");
  tOut->Branch("phi",&phi,"phi/F");
  tOut->Branch("expEl",&expEl,"expEl/F");
  tOut->Branch("expPi",&expPi,"expPi/F");
  tOut->Branch("expKa",&expKa,"expKa/F");
  tOut->Branch("expPr",&expPr,"expPr/F");
  tOut->Branch("weightEl",&weightEl,"weightEl/F");
  tOut->Branch("weightPi",&weightPi,"weightPi/F");
  tOut->Branch("weightKa",&weightKa,"weightKa/F");
  tOut->Branch("weightPr",&weightPr,"weightPr/F");
  tOut->Branch("time",&time,"time/F");
  tOut->Branch("pdg",&particleID,"pdg/I");
  tOut->Branch("distVertT",&distVertT,"distVertT/F");
  tOut->Branch("bayesEl",&bayesEl,"bayesEl/F");
  tOut->Branch("bayesPi",&bayesPi,"bayesPi/F");
  tOut->Branch("bayesKa",&bayesKa,"bayesKa/F");
  tOut->Branch("bayesPr",&bayesPr,"bayesPr/F");

  for(int i=0; i < nev;i ++){ // loop over events
    t->GetEvent(i);

    //    printf("ev=%d) npart = %ld\n",i,part.size());
    for(int j=0;j < part.size();j++){ // loop over particles
      if(mask[j] == 0) continue; // skip particles without a reconstructed track

      particleID = part[j].GetPdgCode(); // get particle identity (MC truth)
      p = part[j].P();                   // get total momentum (MC truth)
      pt = part[j].Pt();                 // get transverse momentum (MC truth)
      eta = part[j].Eta();               // get pseudorapidity (MC truth)
      phi = part[j].Phi();               // get aziumthal angle (MC truth)

      float ptin = pt;
      float etain = eta;
      if(ptin > 9.9999) ptin = 9.9999;
      if(etain > 1.19999) etain = 1.19999;
      if(etain < -1.1999) etain = -1.1999;
      
      int binX = hPriorPi->GetXaxis()->FindBin(ptin);
      int binY = hPriorPi->GetYaxis()->FindBin(etain);

      manager.SetPriors(0,hPriorPi->GetBinContent(binX,binY));
      manager.SetPriors(1,hPriorKa->GetBinContent(binX,binY));
      manager.SetPriors(2,hPriorPr->GetBinContent(binX,binY));
      manager.SetPriors(3,hPriorEl->GetBinContent(binX,binY));

      distVert = part[j].R();
      distVertT = sqrt(part[j].Vx()*part[j].Vx() + part[j].Vy()*part[j].Vy());
      if(distVertT > 3) continue;
      // compute expected velocity for each mass hypotesis
      betaEl = p/sqrt(mEl*mEl + p*p);
      betaPi = p/sqrt(mPi*mPi + p*p);
      betaKa = p/sqrt(mKa*mKa + p*p);
      betaPr = p/sqrt(mPr*mPr + p*p);

      //      if(abs(particleID) != pdgKa) continue; // select only true kaons

      if(mask[j] & 1){ // btof has a track fot this particle

	btof_texp[j] += distVert*3.3356410;
	expEl = btof_texp[j]/betaEl;
	expPi = btof_texp[j]/betaPi;
	expKa = btof_texp[j]/betaKa;
	expPr = btof_texp[j]/betaPr;
	
	meanEl = expEl-expKa;//btof_texp[j]*(1./betaEl - 1./betaKa);
	meanPi = expPi-expKa;//btof_texp[j]*(1./betaPi - 1./betaKa);
	meanPr = expPr-expKa;//btof_texp[j]*(1./betaPr - 1./betaKa);
	manager.SetPar(0,1,meanPi);
	manager.SetPar(1,1,0);
	manager.SetPar(2,1,meanPr);
	manager.SetPar(3,1,meanEl);
	if(meanPi>-10) manager.SetPar(1,1,10000);
	if(meanPr<10) manager.SetPar(2,1,20000);
	if(meanPi-meanEl<10) manager.SetPar(3,1,-10000);
	weightEl =  manager.GetWeight(3,btof_tof[j] - btof_texp[j]/betaKa);
	weightPi =  manager.GetWeight(0,btof_tof[j] - btof_texp[j]/betaKa);
	weightKa =  manager.GetWeight(1,btof_tof[j] - btof_texp[j]/betaKa);
	weightPr =  manager.GetWeight(2,btof_tof[j] - btof_texp[j]/betaKa);
	time = btof_tof[j];

	bayesEl = manager.GetBayesWeight(3,btof_tof[j] - btof_texp[j]/betaKa);
	bayesPi = manager.GetBayesWeight(0,btof_tof[j] - btof_texp[j]/betaKa);
	bayesKa = manager.GetBayesWeight(1,btof_tof[j] - btof_texp[j]/betaKa);
	bayesPr = manager.GetBayesWeight(2,btof_tof[j] - btof_texp[j]/betaKa);

	tOut->Fill();

	if(meanPi-meanEl > 10) hPtBayesEl->Fill(pt,eta,bayesEl);
	hPtBayesPi->Fill(pt,eta,bayesPi);
	if(meanPi < -10) hPtBayesKa->Fill(pt,eta,bayesKa);
	if(meanPr > 10) hPtBayesPr->Fill(pt,eta,bayesPr);

	if(abs(particleID)==pdgPi) hPtTruePi->Fill(pt,eta);
	else if(abs(particleID)==pdgKa) hPtTrueKa->Fill(pt,eta);
	else if(abs(particleID)==pdgPr) hPtTruePr->Fill(pt,eta);
	else if(abs(particleID)==pdgEl) hPtTrueEl->Fill(pt,eta);

	hbtof->Fill(btof_texp[j]/btof_tof[j]);
	htof2D->Fill(p,btof_texp[j]/btof_tof[j]);
	hdeltat2Dpi->Fill(p,btof_tof[j] - btof_texp[j]/betaPi);
	hdeltat2Dka->Fill(p,btof_tof[j] - btof_texp[j]/betaKa);
	hdeltat2Dpr->Fill(p,btof_tof[j] - btof_texp[j]/betaPr);
	hmaskCorr->Fill("btof","btof",1);

	if(mask[j] & 2) hmaskCorr->Fill("btof","etof",1);
	if(mask[j] & 4) hmaskCorr->Fill("btof","htof1",1);
	if(mask[j] & 8) hmaskCorr->Fill("btof","htof2",1);
	if(mask[j] & 16) hmaskCorr->Fill("btof","htof2p",1);
      }
      if(mask[j] & 2){ // etof has a track fot this particle
	htof2D->Fill(p,etof_texp[j]/etof_tof[j]);
	//	hdeltat2Dka->Fill(p,etof_tof[j] - etof_texp[j]/betaKa);
	hmaskCorr->Fill("etof","etof",1);

	if(mask[j] & 4) hmaskCorr->Fill("etof","htof1",1);
	if(mask[j] & 8) hmaskCorr->Fill("etof","htof2",1);
	if(mask[j] & 16) hmaskCorr->Fill("etof","htof2p",1);
      }
      if(mask[j] & 4){ // htof1 has a track fot this particle
	htof2D->Fill(p,htof1_texp[j]/htof1_tof[j]);
	//	hdeltat2Dka->Fill(p,htof1_tof[j] - htof1_texp[j]/betaKa);
	hmaskCorr->Fill("htof1","htof1",1);

	if(mask[j] & 8) hmaskCorr->Fill("htof1","htof2",1);
	if(mask[j] & 16) hmaskCorr->Fill("htof1","htof2p",1);
      }
      if(mask[j] & 8){ // htof2 has a track fot this particle
	htof2D->Fill(p,htof2_texp[j]/htof2_tof[j]);
	//	hdeltat2Dka->Fill(p,htof2_tof[j] - htof2_texp[j]/betaKa);
	hmaskCorr->Fill("htof2","htof2",1);

	if(mask[j] & 16) hmaskCorr->Fill("htof2","htof2p",1);
      }
      if(mask[j] & 16){ // htof2p has a track fot this particle
	htof2D->Fill(p,htof2p_texp[j]/htof2p_tof[j]);
	//	hdeltat2Dka->Fill(p,htof2p_tof[j] - htof2_texp[j]/betaKa);
	hmaskCorr->Fill("htof2p","htof2p",1);
      }
    }
  }
  TCanvas *cbtof = new TCanvas("cbtof","cbtof");
  hbtof->Draw();

  TCanvas *ctof2D = new TCanvas("ctof2D","ctof2D");
  htof2D->Draw("colz");

  TCanvas *cdeltat2Dka = new TCanvas("cdeltat2Dka","cdeltat2Dka");
  hdeltat2Dka->Draw("colz");

  TCanvas *ccorrMask = new TCanvas("ccorrMask","ccorrMAsk");
  hmaskCorr->Draw("box");

  fout->cd();
  tOut->Write();
  fout->Close();

  fout = new TFile("newpriors.root","RECREATE");
  hPtBayesEl->Write();
  hPtBayesPi->Write();
  hPtBayesKa->Write();
  hPtBayesPr->Write();
  hPtBayesEl->ProjectionX()->Write();
  hPtBayesPi->ProjectionX()->Write();
  hPtBayesKa->ProjectionX()->Write();
  hPtBayesPr->ProjectionX()->Write();
  hPtTrueEl->Write();
  hPtTruePi->Write();
  hPtTrueKa->Write();
  hPtTruePr->Write();
  hPtTrueEl->ProjectionX()->Write();
  hPtTruePi->ProjectionX()->Write();
  hPtTrueKa->ProjectionX()->Write();
  hPtTruePr->ProjectionX()->Write();
  fout->Close();
}
