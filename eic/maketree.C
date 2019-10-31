/*
Prints complete input particle arborescence for the first 10 events. Useful for debugging purposes.
root -l examples/Example5.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif


//------------------------------------------------------------------------------

void maketree(const char *inputFile="delphes.root")
{
  gSystem->Load("libDelphes");
  gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchBTOF = treeReader->UseBranch("bTOF");
  TClonesArray *branchETOF = treeReader->UseBranch("eTOF");
  TClonesArray *branchHTOF1 = treeReader->UseBranch("hTOF1");
  TClonesArray *branchHTOF2 = treeReader->UseBranch("hTOF2");
  TClonesArray *branchHTOF2p = treeReader->UseBranch("hTOF2p");

  // output structure
  TFile *fout = new TFile("my_output.root","RECREATE");
  TTree *treeOut = new TTree("tree","tree");
  std::vector<TParticle> part;
  std::vector<int> trackMask;
  treeOut->Branch("particle",&part);
  treeOut->Branch("trackMask",&trackMask);

  std::vector<float> btof_tof;
  std::vector<float> btof_texp;
  treeOut->Branch("btof_tof",&btof_tof);
  treeOut->Branch("btof_texp",&btof_texp);

  std::vector<float> etof_tof;
  std::vector<float> etof_texp;
  treeOut->Branch("etof_tof",&etof_tof);
  treeOut->Branch("etof_texp",&etof_texp);

  std::vector<float> htof1_tof;
  std::vector<float> htof1_texp;
  treeOut->Branch("htof1_tof",&htof1_tof);
  treeOut->Branch("htof1_texp",&htof1_texp);

  std::vector<float> htof2_tof;
  std::vector<float> htof2_texp;
  treeOut->Branch("htof2_tof",&htof2_tof);
  treeOut->Branch("htof2_texp",&htof2_texp);

  std::vector<float> htof2p_tof;
  std::vector<float> htof2p_texp;
  treeOut->Branch("htof2p_tof",&htof2p_tof);
  treeOut->Branch("htof2p_texp",&htof2p_texp);

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    part.clear();
    trackMask.clear();
    btof_tof.clear();
    btof_texp.clear();    
    etof_tof.clear();
    etof_texp.clear();
    htof1_tof.clear();
    htof1_texp.clear();
    htof2_tof.clear();
    htof2_texp.clear();
    htof2p_tof.clear();
    htof2p_texp.clear();

//    if(entry>10) break;
    
//    cout<<"" <<endl;
//    cout<<"--------- New Event ---------" <<endl;
//    cout<<"" <<endl;
 
    double tof,texp;
    // loop over all input particles in the event
    for(Int_t i=0; i < branchParticle->GetEntriesFast(); i++)
    {    
     GenParticle *gen = (GenParticle*) branchParticle->At(i);     
     gen->IsPU = -1;
     if(gen->PID == 0 || gen->PID > 10000) continue;
     gen->IsPU = part.size();
     part.emplace_back(gen->PID, gen->Status, gen->M1, gen->M2, gen->D1, gen->D2, gen->Px, gen->Py, gen->Pz, gen->E, gen->X, gen->Y, gen->Z, gen->T);
     trackMask.emplace_back(0);
     btof_tof.emplace_back(-1);
     btof_texp.emplace_back(-1);
     etof_tof.emplace_back(-1);
     etof_texp.emplace_back(-1);
     htof1_tof.emplace_back(-1);
     htof1_texp.emplace_back(-1);
     htof2_tof.emplace_back(-1);
     htof2_texp.emplace_back(-1);
     htof2p_tof.emplace_back(-1);
     htof2p_texp.emplace_back(-1);
//     status = gen->Status;
//     pid = gen->PID;
//     px = gen->Px;
//     py = gen->Py;
//     pz = gen->Pz;
//     cout<<"N: "<<i<<", St: "<<gen->Status<<", PID: "<<gen->PID<<", E: "<<gen->E<<", Px: "<<gen->Px<<", Py: "<<gen->Py<<", Pz: "<<gen->Pz<<", M: "<<gen->Mass<<", M1: "<<gen->M1<<", M2: "<<gen->M2<<", D1: "<<gen->D1<<", D2: "<<gen->D2<<endl;
    }
    for(Int_t i=0; i < branchBTOF->GetEntriesFast(); i++)
    {
      Track *tr = (Track*) branchBTOF->At(i);
      auto particle = (GenParticle *)tr->Particle.GetObject();
      tof = tr->TOuter * 1.e12; // in ps
      texp = tr->L * 3.3356410; // in ps
      if(particle->IsPU==-1){
        printf("no particle for btof %f %f\n",tof,texp);
        continue;
      }
      trackMask[particle->IsPU] += 1;
      btof_tof[particle->IsPU] = tof;
      btof_texp[particle->IsPU] = texp;
    }
    for(Int_t i=0; i < branchETOF->GetEntriesFast(); i++)
    {
      Track *tr = (Track*) branchETOF->At(i);
      auto particle = (GenParticle *)tr->Particle.GetObject();
      tof = tr->TOuter * 1.e12; // in ps
      texp = tr->L * 3.3356410; // in ps
      if(particle->IsPU==-1){
        printf("no particle for etof %f %f\n",tof,texp);
        continue;
      }
      trackMask[particle->IsPU] += 2;
      etof_tof[particle->IsPU] = tof;
      etof_texp[particle->IsPU] = texp;
    }
    for(Int_t i=0; i < branchHTOF1->GetEntriesFast(); i++)
    {
      Track *tr = (Track*) branchHTOF1->At(i);
      auto particle = (GenParticle *)tr->Particle.GetObject();
      tof = tr->TOuter * 1.e12; // in ps
      texp = tr->L * 3.3356410; // in ps
      if(particle->IsPU==-1){
        printf("no particle for htof1 %f %f\n",tof,texp);
        continue;
      }
      trackMask[particle->IsPU] += 4;
      htof1_tof[particle->IsPU] = tof;
      htof1_texp[particle->IsPU] = texp;
    }
    for(Int_t i=0; i < branchHTOF2->GetEntriesFast(); i++)
    {
      Track *tr = (Track*) branchHTOF2->At(i);
      auto particle = (GenParticle *)tr->Particle.GetObject();
      tof = tr->TOuter * 1.e12; // in ps
      texp = tr->L * 3.3356410; // in ps
      if(particle->IsPU==-1){
        printf("no particle for htof2 %f %f\n",tof,texp);
        continue;
      }
      trackMask[particle->IsPU] += 8;
      htof2_tof[particle->IsPU] = tof;
      htof2_texp[particle->IsPU] = texp;
    }
    for(Int_t i=0; i < branchHTOF2p->GetEntriesFast(); i++)
    {
      Track *tr = (Track*) branchHTOF2p->At(i);
      auto particle = (GenParticle *)tr->Particle.GetObject();
      tof = tr->TOuter * 1.e12; // in ps
      texp = tr->L * 3.3356410; // in ps
      if(particle->IsPU==-1){
        printf("no particle for htof2p %f %f\n",tof,texp);
        continue;
      }
      trackMask[particle->IsPU] += 16;
      htof2p_tof[particle->IsPU] = tof;
      htof2p_texp[particle->IsPU] = texp;
    }

    treeOut->Fill();
  }

  treeOut->Write();
  fout->Close();
}
