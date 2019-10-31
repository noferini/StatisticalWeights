#include "TF1.h"
#include "TMatrix.h"
#include "TMath.h"

class StatisticalWeights{
public:
  int GetNspecies() const {return fNspecies;}
  void AddSpecies(const char *name,const char *formula, float priors=1);
  virtual void AddSpecies(const char *name,float *par, float priors=1);
  virtual void AddSignal(int species, const char* formula);
  virtual void AddSignal(int species,float *par){printf("Parameteric add signal not vali for base class\n");}

  TF1* GetResponse(int i,int j=0) {return fResponse[i][j];} 
  const char* GetName(int i) const {if(!fName[i]) return "";return fName[i]->Data();} 

  void SetRange(Float_t xmin, Float_t xmax) {fXmin=xmin; fXmax = xmax;
    for(int i=0; i < fNspecies;i++)   fResponse[i][0]->SetRange(fXmin,fXmax); }
  void SetXbin(Float_t xbin) {fXbin=xbin;}
  Float_t GetXmin() const {return fXmin;}
  Float_t GetXmax() const {return fXmax;}
  Float_t GetXbin() const {return fXbin;}

  virtual void SetPar(int isp, int par, float value, int isig=0);

  virtual void Init();
  virtual Float_t Eval(Int_t i, Float_t x, Int_t isig=0);

  Float_t GetWeight(Int_t i, Float_t x);
  Float_t GetWeight(Int_t i, Float_t *x);
  virtual Float_t GetBayesWeight(Int_t i, Float_t x);

  void SetPriors(Int_t i, Float_t x){if(i >= fNspecies) return; fPriors[i] = x;}

  virtual Double_t GetScalarProduct(Int_t i, Int_t j);

protected:
  TMatrix *fMat = NULL;
  TMatrix *fMatInv = NULL;
  Bool_t fToBeInitialized = kTRUE;
  Bool_t fRedoMatrix = kTRUE;
  int fNspecies=0;
  int fNsignals=0;
  int fSignalFilled[100];
  void SetResponse(int i, int j, TF1 *f){fResponse[i][j] = f;}
private:
  Float_t fXmin = -5;
  Float_t fXmax =  5;
  Float_t fXbin = 0.2;

  TF1 *fResponse[100][10];
  TString *fName[100];
  float fPriors[100];
};

class StatisticalWeightsGaus : public StatisticalWeights{
public:
  void AddSignal(int species, const char* formula){printf("Invalid AddSignal with formula for Gaussian derived class\n");}
  void AddSignal(int species,float *par);
  Double_t GetScalarProduct(Int_t i, Int_t j);
  void SetPar(int isp, int par, float value, int isig=0);
  void Init();
  Float_t Eval(Int_t i, Float_t x, Int_t isig=0) {return TMath::Gaus(x,fMean[i][isig],fSigma[i][isig])*0.39894228*fSigmaInv[i][isig];}
private:
  float fMean[100][100];
  float fSigma[100][100];
  float fSigmaInv[100][100];
};

class StatisticalManager{
};


void StatisticalWeights::AddSpecies(const char *name,const char *formula, float priors){
  if(fNspecies == 100) return;
  
  fName[fNspecies] = new TString(name);
  fSignalFilled[fNspecies]=0;
  fPriors[fNspecies]=priors;
  fNspecies++;
  fToBeInitialized = kTRUE;
  fRedoMatrix = kTRUE;

  AddSignal(fNspecies-1,formula);

  if(! fResponse[fNspecies-1][0])
    fNspecies--;
  
}

void StatisticalWeights::AddSignal(int species, const char* formula){
  if(species >= fNspecies || fSignalFilled[species]==100) return;
  fResponse[species][fSignalFilled[species]] = new TF1(Form("f_%s_%d",fName[species]->Data(),fSignalFilled[species]),Form("(%s)*[0]",formula));

  fResponse[species][fSignalFilled[species]]->SetParameter(0,1);
  fResponse[species][fSignalFilled[species]]->SetRange(fXmin,fXmax);
  fSignalFilled[species]++;
}

void StatisticalWeights::Init(){
  if(fRedoMatrix){
    if(fMat) delete fMat;
    if(fMatInv) delete fMatInv;
    
    fMat = new TMatrix(fNspecies,fNspecies);
    fMatInv = new TMatrix(fNspecies,fNspecies);

    for(Int_t i=0; i < fNspecies; i++){
      for(int k=0; k < fSignalFilled[i];k++){
	fResponse[i][k]->SetRange(fXmin,fXmax);
	int np = (fXmax - fXmin)/fXbin;
	fResponse[i][k]->SetNpx(np);
	fResponse[i][k]->SetParameter(0,1);
	fResponse[i][k]->SetParameter(0,1./fResponse[i][k]->Integral(fXmin,fXmax));
      }
    }
  }

  for(Int_t i=0; i < fNspecies; i++){
    for(Int_t j=i; j < fNspecies; j++){
      (*fMat)[i][j] = this->GetScalarProduct(i,j);
      (*fMat)[j][i] = (*fMat)[i][j];
      (*fMatInv)[i][j] = (*fMat)[i][j];
      (*fMatInv)[j][i] = (*fMat)[i][j];
    }
  }

  fMatInv->Invert();

  if(fRedoMatrix){
    fMat->Print();
    fMatInv->Print();
  }
  fToBeInitialized = kFALSE;
  fRedoMatrix = kFALSE;
}


void StatisticalWeightsGaus::Init(){
  if(fRedoMatrix){
    if(fMat) delete fMat;
    if(fMatInv) delete fMatInv;
    
    fMat = new TMatrix(fNspecies,fNspecies);
    fMatInv = new TMatrix(fNspecies,fNspecies);
  }


  for(Int_t i=0; i < fNspecies; i++){
    for(Int_t j=i; j < fNspecies; j++){
      (*fMat)[i][j] = this->GetScalarProduct(i,j);
      (*fMat)[j][i] = (*fMat)[i][j];
      (*fMatInv)[i][j] = (*fMat)[i][j];
      (*fMatInv)[j][i] = (*fMat)[i][j];
    }
  }

  fMatInv->Invert();

  if(fRedoMatrix){
    fMat->Print();
    fMatInv->Print();
  }
  fToBeInitialized = kFALSE;
  fRedoMatrix = kFALSE;
}

Double_t StatisticalWeights::GetScalarProduct(Int_t i, Int_t j){
  if(i >= fNspecies || j >= fNspecies) return 0;

  double res=1;
  for(int k=0;k < fSignalFilled[i];k++){
    Float_t x = fXmin + fXbin*0.5;
    double integral=0;
    while(x < fXmax){
      integral += Eval(i,x,k)*Eval(j,x,k);
      x += fXbin;
    }
    res *= integral*fXbin;
  }

  return res;
}

Double_t StatisticalWeightsGaus::GetScalarProduct(Int_t i, Int_t j){
  if(i >= GetNspecies() || j >= GetNspecies()) return 0;

  float res=1;
  for(int k=0;k < fSignalFilled[i];k++){
    float sigmainv = sqrt(2./(fSigma[i][k]*fSigma[i][k]+fSigma[j][k]*fSigma[j][k]));
    float nsigma = (fMean[i][k] - fMean[j][k])*sigmainv;
    res *= TMath::Exp(-nsigma*nsigma*0.25)*0.28209479*sigmainv;
  }
  return res;
}

Float_t StatisticalWeights::GetWeight(Int_t i, Float_t x){
  if(i >= fNspecies) return 0;
  if(fToBeInitialized) Init();

  Float_t weight = 0;
  for(Int_t j=0; j < fNspecies; j++){
    //    printf("%f %f\n",(*fMatInv)[i][j], Eval(j,x));
    weight += (*fMatInv)[i][j] * Eval(j,x);
  }

  return weight;
}

Float_t StatisticalWeights::GetWeight(Int_t i, Float_t *x){
  if(i >= fNspecies) return 0;
  if(fToBeInitialized) Init();
  Float_t weight = 0;
  for(Int_t j=0; j < fNspecies; j++){
    float val=1;
    for(int k=0;k < fSignalFilled[j];k++)
      val *= Eval(j,x[k],k);
    
    weight += (*fMatInv)[i][j] * val;
  }

  return weight;
}

Float_t StatisticalWeights::GetBayesWeight(Int_t i, Float_t x){
  if(i >= fNspecies) return 0;
  Float_t sum = 0;
  for(Int_t j=0; j < fNspecies; j++){
    sum += fPriors[j] * this->Eval(j,x);
  }
  if(sum==0) return 0;

  return fPriors[i] * this->Eval(i,x)/sum;
}

Float_t StatisticalWeights::Eval(Int_t i, Float_t x, Int_t isig){
  if(i >= fNspecies) return 0;
  if(isig >= fSignalFilled[i]) return 0;
  return fResponse[i][isig]->Eval(x);
}

void StatisticalWeights::SetPar(int isp, int par, float value, int isig){
  if(isp >= fNspecies) return;
  fResponse[isp][isig]->SetParameter(par,value);
  fToBeInitialized = kTRUE;
}

void StatisticalWeightsGaus::SetPar(int isp, int par, float value, int isig){
  if(isp >= fNspecies) return;
  if(par == 1) fMean[isp][isig] = value;
  else if(par == 2){
    fSigma[isp][isig] = value;
    fSigmaInv[isp][isig] = 1./value;
  }
  // GetResponse(isp,isig)->SetParameter(par,value);
  fToBeInitialized = kTRUE;
}

void StatisticalWeightsGaus::AddSignal(int species,float *par){
  if(species >= fNspecies || fSignalFilled[species]==100) return;

  SetResponse(species,fSignalFilled[species], new TF1(Form("f_%s_%d",GetName(species),fSignalFilled[species]),"TMath::Gaus(x,[1],[2])"));
  GetResponse(species,fSignalFilled[species])->SetParameter(0,1);
  GetResponse(species,fSignalFilled[species])->SetParameter(1,par[0]);
  GetResponse(species,fSignalFilled[species])->SetParameter(2,par[1]);
  GetResponse(species,fSignalFilled[species])->SetRange(GetXmin(),GetXmax());

  fMean[species][fSignalFilled[species]] = par[0];
  fSigma[species][fSignalFilled[species]] = par[1];
  fSigmaInv[species][fSignalFilled[species]] = 1./par[1];

  fSignalFilled[species]++;

}

void StatisticalWeights::AddSpecies(const char *name,float *par, float priors){
  if(fNspecies == 100) return;
  
  fName[fNspecies] = new TString(name);
  fSignalFilled[fNspecies]=0;

  fToBeInitialized = kTRUE;
  fRedoMatrix = kTRUE;
  fPriors[fNspecies]=priors;

  fNspecies++;
  AddSignal(fNspecies-1,par);
}
