#include "TF1.h"
#include "TMatrix.h"
#include "TMath.h"
#include "TH1F.h"

class StatisticalWeights{
public:
  int GetNspecies() const {return fNspecies;}
  void AddSpecies(const char *name,const char *formula, double priors=1);
  virtual void AddSpecies(const char *name,double *par, double priors=1);
  virtual void AddSpecies(const char *name,TH1F *h, double priors=1);
  virtual void AddSignal(int species, const char* formula);
  virtual void AddSignal(int species,double *par){printf("Parameteric add signal not vali for base class\n");}
  virtual void AddSignal(int species,TH1F *h){printf("TH1F add signal not vali for base class\n");}

  TF1* GetResponse(int i,int j=0) {return fResponse[i][j];} 
  const char* GetName(int i) const {if(!fName[i]) return "";return fName[i]->Data();} 

  void SetRange(Double_t xmin, Double_t xmax) {fXmin=xmin; fXmax = xmax;
    for(int i=0; i < fNspecies;i++)   fResponse[i][0]->SetRange(fXmin,fXmax); }
  void SetXbin(Double_t xbin) {fXbin=xbin;}
  Double_t GetXmin() const {return fXmin;}
  Double_t GetXmax() const {return fXmax;}
  Double_t GetXbin() const {return fXbin;}

  virtual void SetPar(int isp, int par, double value, int isig=0);

  virtual void Init();
  virtual Double_t Eval(Int_t i, Double_t x, Int_t isig=0);

  Double_t GetWeight(Int_t i, Double_t x);
  Double_t GetWeight(Int_t i, Double_t *x);
  Double_t GetBayesWeight(Int_t i, Double_t x);
  Double_t GetBayesWeight(Int_t i, Double_t *x);

  void SetPriors(Int_t i, Double_t x){if(i >= fNspecies) return; fPriors[i] = x;}

  virtual Double_t GetScalarProduct(Int_t i, Int_t j);

  void SetMaskSignal(int i, bool val) {fToBeInitialized = kTRUE; fMaskSignal[i] = val;}

  TMatrixD *GetMatrix() {return fMat;}

protected:
  TMatrixD *fMat = NULL;
  TMatrixD *fMatInv = NULL;
  Bool_t fToBeInitialized = kTRUE;
  Bool_t fRedoMatrix = kTRUE;
  int fNspecies=0;
  int fNsignals=0;
  int fSignalFilled[100];
  bool fMaskSignal[100];
  void SetResponse(int i, int j, TF1 *f){fResponse[i][j] = f;}
private:
  Double_t fXmin = -5;
  Double_t fXmax =  5;
  Double_t fXbin = 0.2;

  TF1 *fResponse[100][100];
  TString *fName[100];
  double fPriors[100];
};

class StatisticalWeightsGaus : public StatisticalWeights{
public:
  void AddSignal(int species, const char* formula){printf("Invalid AddSignal with formula for Gaussian derived class\n");}
  void AddSignal(int species,double *par);
  Double_t GetScalarProduct(Int_t i, Int_t j);
  void SetPar(int isp, int par, double value, int isig=0);
  void Init();
  Double_t Eval(Int_t i, Double_t x, Int_t isig=0) {return TMath::Gaus(x,fMean[i][isig],fSigma[i][isig])*0.39894228*fSigmaInv[i][isig];}
private:
  double fMean[100][100];
  double fSigma[100][100];
  double fSigmaInv[100][100];
};

class StatisticalWeightsHisto : public StatisticalWeights{
public:
  void AddSignal(int species, const char* formula){printf("Invalid AddSignal with formula for Histo derived class\n");}
  void AddSignal(int species,double *par) {printf("Invalid AddSignal with parameters for Histo derived class\n");}
  void AddSignal(int species, TH1F *h);
  Double_t GetScalarProduct(Int_t i, Int_t j);
  void Init();
  Double_t Eval(Int_t i, Double_t x, Int_t isig=0) {return fH[i][isig]->GetBinContent(fH[i][isig]->FindBin(x));}
private:
  TH1F *fH[100][100];
};

class StatisticalManager{
public:
  bool AddDiscriminator(StatisticalWeights *dis) {if(fNum < 4999){ fDiscriminators[fNum] = dis; fNum++; fRedoMatrix=kTRUE; fToBeInitialized=kTRUE; return 1;} else return 0;} 
  Double_t GetWeight(Int_t i, Double_t *x);
  void Init();
private:
  int fNum = 0;
  StatisticalWeights *fDiscriminators[5000];
  TMatrixD *fMat = NULL;
  double fDeterminant;
  TMatrixD *fMatInv = NULL;
  Bool_t fToBeInitialized = kTRUE;
  Bool_t fRedoMatrix = kTRUE;
};



void StatisticalWeightsHisto::AddSignal(int species, TH1F *h){
  if(species >= fNspecies || fSignalFilled[species]==100) return;
  h->Scale(1./h->Integral());
  fH[species][fSignalFilled[species]] = h;
  fMaskSignal[fSignalFilled[species]] = 1;
  fSignalFilled[species]++;
}

Double_t StatisticalWeightsHisto::GetScalarProduct(Int_t i, Int_t j){
  if(i >= fNspecies || j >= fNspecies) return 0;
  double sp=1;
  for(int k=0;k < fSignalFilled[i];k++){
    double sum=0;
    if(! fMaskSignal[k]) continue;
    for(int ii=1; ii <= fH[i][k]->GetNbinsX(); ii++){
      sum += fH[i][k]->GetBinContent(ii) * fH[j][k]->GetBinContent(ii);
    }
    sp *= sum;
  }
  return sp;
}

void StatisticalWeights::AddSpecies(const char *name,const char *formula, double priors){
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

void StatisticalWeights::AddSpecies(const char *name,TH1F *h, double priors){
  if(fNspecies == 100) return;
  
  fName[fNspecies] = new TString(name);
  fSignalFilled[fNspecies]=0;
  fPriors[fNspecies]=priors;
  fNspecies++;
  fToBeInitialized = kTRUE;
  fRedoMatrix = kTRUE;

  AddSignal(fNspecies-1,h);
}

void StatisticalWeights::AddSignal(int species, const char* formula){
  if(species >= fNspecies || fSignalFilled[species]==100) return;
  fResponse[species][fSignalFilled[species]] = new TF1(Form("f_%s_%d",fName[species]->Data(),fSignalFilled[species]),Form("(%s)*[0]",formula));

  fResponse[species][fSignalFilled[species]]->SetParameter(0,1);
  fResponse[species][fSignalFilled[species]]->SetRange(fXmin,fXmax);
  fMaskSignal[fSignalFilled[species]] = 1;
  fSignalFilled[species]++;
}

void StatisticalWeightsHisto::Init(){
  if(fRedoMatrix){
    if(fMat) delete fMat;
    if(fMatInv) delete fMatInv;
    
    fMat = new TMatrixD(fNspecies,fNspecies);
    fMatInv = new TMatrixD(fNspecies,fNspecies);
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

  if(fRedoMatrix && 0){
    fMat->Print();
    fMatInv->Print();
  }
  fToBeInitialized = kFALSE;
  fRedoMatrix = kFALSE;
}

void StatisticalWeights::Init(){
  if(fRedoMatrix){
    if(fMat) delete fMat;
    if(fMatInv) delete fMatInv;
    
    fMat = new TMatrixD(fNspecies,fNspecies);
    fMatInv = new TMatrixD(fNspecies,fNspecies);

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

  if(fRedoMatrix && 0){
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
    
    fMat = new TMatrixD(fNspecies,fNspecies);
    fMatInv = new TMatrixD(fNspecies,fNspecies);
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

  if(fRedoMatrix && 0){
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
    if(! fMaskSignal[k]) continue;
    Double_t x = fXmin + fXbin*0.5;
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

  double res=1;
  for(int k=0;k < fSignalFilled[i];k++){
    if(! fMaskSignal[k]) continue;
   double sigmainv = sqrt(2./(fSigma[i][k]*fSigma[i][k]+fSigma[j][k]*fSigma[j][k]));
    double nsigma = (fMean[i][k] - fMean[j][k])*sigmainv;
    res *= TMath::Exp(-nsigma*nsigma*0.25)*0.28209479*sigmainv;
  }
  return res;
}

Double_t StatisticalWeights::GetWeight(Int_t i, Double_t x){
  if(i >= fNspecies) return 0;
  if(fToBeInitialized) Init();

  Double_t weight = 0;
  for(Int_t j=0; j < fNspecies; j++){
    //    printf("%f %f\n",(*fMatInv)[i][j], Eval(j,x));
    weight += (*fMatInv)[i][j] * Eval(j,x);
  }

  return weight;
}

void StatisticalManager::Init(){
  if (fNum < 1 || fNum > 0) return;
  if(fRedoMatrix){
    fDeterminant = 1;
    if(fMat) delete fMat;
    if(fMatInv) delete fMatInv;
    fDiscriminators[0]->Init();
    TMatrixD *mm = fDiscriminators[0]->GetMatrix(); 
    fMat = new TMatrixD(fDiscriminators[0]->GetNspecies(),fDiscriminators[0]->GetNspecies());
    *fMat = *mm;
    fMatInv = new TMatrixD(fDiscriminators[0]->GetNspecies(),fDiscriminators[0]->GetNspecies());
    
   for(int i=1; i < fNum; i++){
      fDiscriminators[i]->Init();    

      if(fMat->Determinant() < 1E-1 || fMat->Determinant() > 1E1){
        double det = fMat->Determinant();
        fDeterminant *= det;
        *fMat *= 1./det;       
      }

      for(int ii=0; ii < fDiscriminators[i]->GetNspecies(); ii++)
	for(int jj=0; jj < fDiscriminators[i]->GetNspecies(); jj++)
	  (*fMat)[ii][jj] *= (*(fDiscriminators[i]->GetMatrix()))[ii][jj];
    } 
  }

  (*fMatInv) = (*fMat);
  fMatInv->Invert();
  
  if(fRedoMatrix){
    fMat->Print();
    fMatInv->Print();
  }
  fToBeInitialized = kFALSE;
  fRedoMatrix = kFALSE;
}

Double_t StatisticalManager::GetWeight(Int_t i, Double_t *x){
  if(fToBeInitialized) Init();

  Double_t weight = 0;

  if( fNum > 0 ){
    for(int iDis=0; iDis < fNum; iDis++){
      weight += fDiscriminators[iDis]->GetWeight(i, x[iDis]);
    }

    return weight/fNum;
  }

  for(Int_t j=0; j < fDiscriminators[0]->GetNspecies(); j++){
    double resp = 1;
    for(int iDis=0; iDis < fNum; iDis++){
      resp *= fDiscriminators[iDis]->Eval(j,x[iDis]);
    }
    weight += (*fMatInv)[i][j] * resp;
  }
//printf("det = %e %e\n",weight,fDeterminant);
  return weight/fDeterminant;
}

Double_t StatisticalWeights::GetWeight(Int_t i, Double_t *x){
  if(i >= fNspecies) return 0;
  if(fToBeInitialized) Init();
  Double_t weight = 0;
  for(Int_t j=0; j < fNspecies; j++){
    double val=1;
    for(int k=0;k < fSignalFilled[j];k++)
      if(fMaskSignal[k])
	val *= Eval(j,x[k],k);
    
    weight += (*fMatInv)[i][j] * val;
  }

  return weight;
}

Double_t StatisticalWeights::GetBayesWeight(Int_t i, Double_t x){
  if(i >= fNspecies) return 0;
  Double_t sum = 0;
  for(Int_t j=0; j < fNspecies; j++){  
    sum += fPriors[j] * this->Eval(j,x);
  }
  if(sum==0) return 0;

  return fPriors[i] * this->Eval(i,x)/sum;
}

Double_t StatisticalWeights::GetBayesWeight(Int_t i, Double_t *x){
  if(i >= fNspecies) return 0;
  Double_t sum = 0;
  double cval = 1;
  for(Int_t j=0; j < fNspecies; j++){
    double val = 1;
    for(int k=0;k < fSignalFilled[j];k++){
      if(! fMaskSignal[k]) continue;
      val *= this->Eval(j,x[k],k);     
      if(i==j) cval *= this->Eval(j,x[k],k);
    }
   
    sum += fPriors[j] * val;
  }
  if(sum==0) return 0;

  return fPriors[i] * cval/sum;
}

Double_t StatisticalWeights::Eval(Int_t i, Double_t x, Int_t isig){
  if(i >= fNspecies) return 0;
  if(isig >= fSignalFilled[i]) return 0;
  return fResponse[i][isig]->Eval(x);
}

void StatisticalWeights::SetPar(int isp, int par, double value, int isig){
  if(isp >= fNspecies) return;
  fResponse[isp][isig]->SetParameter(par,value);
  fToBeInitialized = kTRUE;
}

void StatisticalWeightsGaus::SetPar(int isp, int par, double value, int isig){
  if(isp >= fNspecies) return;
  if(par == 1) fMean[isp][isig] = value;
  else if(par == 2){
    fSigma[isp][isig] = value;
    fSigmaInv[isp][isig] = 1./value;
  }
  // GetResponse(isp,isig)->SetParameter(par,value);
  fToBeInitialized = kTRUE;
}

void StatisticalWeightsGaus::AddSignal(int species,double *par){
  if(species >= fNspecies || fSignalFilled[species]==100) return;

  SetResponse(species,fSignalFilled[species], new TF1(Form("f_%s_%d",GetName(species),fSignalFilled[species]),"TMath::Gaus(x,[1],[2])"));
  GetResponse(species,fSignalFilled[species])->SetParameter(0,1);
  GetResponse(species,fSignalFilled[species])->SetParameter(1,par[0]);
  GetResponse(species,fSignalFilled[species])->SetParameter(2,par[1]);
  GetResponse(species,fSignalFilled[species])->SetRange(GetXmin(),GetXmax());

  fMean[species][fSignalFilled[species]] = par[0];
  fSigma[species][fSignalFilled[species]] = par[1];
  fSigmaInv[species][fSignalFilled[species]] = 1./par[1];
  fMaskSignal[fSignalFilled[species]] = 1;

  fSignalFilled[species]++;

}

void StatisticalWeights::AddSpecies(const char *name,double *par, double priors){
  if(fNspecies == 100) return;
  
  fName[fNspecies] = new TString(name);
  fSignalFilled[fNspecies]=0;

  fToBeInitialized = kTRUE;
  fRedoMatrix = kTRUE;
  fPriors[fNspecies]=priors;

  fNspecies++;
  AddSignal(fNspecies-1,par);
}
