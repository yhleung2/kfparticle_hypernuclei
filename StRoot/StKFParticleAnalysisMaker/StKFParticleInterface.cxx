//TFG18n
#include "StKFParticleInterface.h"

#include "KFParticleTopoReconstructor.h"
#include "KFMCTrack.h"

#include "TMath.h"
#include "TArrayD.h"
#include "TH1F.h"
#include "TH2F.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoArrays.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"

#include "StBichsel/Bichsel.h"
#include "StBichsel/StdEdxModel.h"
#include "StProbPidTraits.h"
#include "StMuDSTMaker/COMMON/StMuBTofHit.h"
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"

#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StKFParticleInterface);
StKFParticleInterface *StKFParticleInterface::fgStKFParticleInterface = 0;
StKFParticleInterface::StKFParticleInterface(): 
  fKFParticleTopoReconstructor(0), fParticles(0), fParticlesPdg(0), fNHftHits(0),
  fCollectTrackHistograms(false), fCollectPIDHistograms(false),
  fStrictTofPID(true), fCleanKaonsWitTof(true), fdEdXMode(1), fTriggerMode(false),
  fChiPrimaryCut(18.6), fChiPrimaryMaxCut(2e4), fCleanLowPVTrackEvents(false), fUseHFTTracksOnly(false)
{
  fKFParticleTopoReconstructor = new KFParticleTopoReconstructor();
  fgStKFParticleInterface = this;
  // set default cuts
  SetPrimaryProbCut(0.0001); // 0.01% to consider primary track as a secondary;
}

StKFParticleInterface::~StKFParticleInterface()
{  
  if(fKFParticleTopoReconstructor) delete fKFParticleTopoReconstructor;
  fgStKFParticleInterface = 0;
}

void StKFParticleInterface::SetField(float field) 
{ 
  if(fKFParticleTopoReconstructor)
    fKFParticleTopoReconstructor->SetField(field); 
}
void StKFParticleInterface::SetBeamLine(KFParticle& p)
{ 
  if(fKFParticleTopoReconstructor)
    fKFParticleTopoReconstructor->SetBeamLine(p);
}

void StKFParticleInterface::InitParticles()
{ 
  fKFParticleTopoReconstructor->Init( fParticles, &fParticlesPdg, &fNHftHits );
  Int_t NPV =  fKFParticleTopoReconstructor->NPrimaryVertices();
  fKFParticleTopoReconstructor->GetKFParticleFinder()->Init(NPV);
  fKFParticleTopoReconstructor->FillPVIndices();
}

void StKFParticleInterface::ReconstructParticles()
{ 
  fKFParticleTopoReconstructor->SortTracks();
  fKFParticleTopoReconstructor->ReconstructParticles();
  
//   static int iEvent=0;
//   iEvent++;
//   std::cout << "Event " << iEvent << ": init " << fKFParticleTopoReconstructor->StatTime( 0 ) 
//             << " pv " << fKFParticleTopoReconstructor->StatTime( 1 )
//             << " sort " << fKFParticleTopoReconstructor->StatTime( 2 )
//             << " particles " << fKFParticleTopoReconstructor->StatTime( 3 ) <<  std::endl;
}

void StKFParticleInterface::ReconstructTopology()
{ 
  fKFParticleTopoReconstructor->Init( fParticles, &fParticlesPdg );
  fKFParticleTopoReconstructor->ReconstructPrimVertex(0);
  fKFParticleTopoReconstructor->SortTracks();
  fKFParticleTopoReconstructor->ReconstructParticles();
  
//   static int iEvent=0;
//   iEvent++;
//   std::cout << "Event " << iEvent << ": init " << fKFParticleTopoReconstructor->StatTime( 0 ) 
//             << " pv " << fKFParticleTopoReconstructor->StatTime( 1 )
//             << " sort " << fKFParticleTopoReconstructor->StatTime( 2 )
//             << " particles " << fKFParticleTopoReconstructor->StatTime( 3 ) <<  std::endl;
}

void StKFParticleInterface::AddPV(const KFVertex &pv, const vector<int> &tracks) { 
  fKFParticleTopoReconstructor->AddPV(pv, tracks);
  fKFParticleTopoReconstructor->FillPVIndices();
}
void StKFParticleInterface::CleanPV() {
  fKFParticleTopoReconstructor->CleanPV();
}

void StKFParticleInterface::AddPV(const KFVertex &pv) { 
  fKFParticleTopoReconstructor->AddPV(pv);
}

void StKFParticleInterface::AddParticle(const KFParticle &p) { 
  fKFParticleTopoReconstructor->AddParticle(p);
}

void StKFParticleInterface::AddCandidate(const KFParticle& candidate, int iPV) {
  fKFParticleTopoReconstructor->AddCandidate(candidate, iPV);
}
void StKFParticleInterface::AddDecayToReconstructionList(Int_t pdg) {
   fKFParticleTopoReconstructor->GetKFParticleFinder()->AddDecayToReconstructionList(pdg);
}
std::vector<KFParticle> const &StKFParticleInterface::GetParticles() const { return fKFParticleTopoReconstructor->GetParticles(); }
void StKFParticleInterface::RemoveParticle(const int iParticle) { fKFParticleTopoReconstructor->RemoveParticle(iParticle); }
const std::vector<KFParticle>* StKFParticleInterface::GetSecondaryCandidates() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetSecondaryCandidates();                           } // Get secondary particles with the mass constraint
const std::vector<KFParticle>& StKFParticleInterface::GetSecondaryK0() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetSecondaryK0();                           } // Get secondary particles with the mass constraint
const std::vector<KFParticle>& StKFParticleInterface::GetSecondaryLambda() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetSecondaryLambda();                           } // Get secondary particles with the mass constraint
const std::vector<KFParticle>& StKFParticleInterface::GetSecondaryAntiLambda() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetSecondaryAntiLambda();                           } // Get secondary particles with the mass constraint
const std::vector<KFParticle>& StKFParticleInterface::GetSecondaryGamma() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetSecondaryGamma();                           } // Get secondary particles with the mass constraint
const std::vector<KFParticle>& StKFParticleInterface::GetSecondaryPi0() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetSecondaryPi0();                           } // Get secondary particles with the mass constraint
const std::vector< std::vector<KFParticle> >* StKFParticleInterface::GetPrimaryCandidates() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetPrimaryCandidates();                } // Get primary particles with the mass constraint
const std::vector< std::vector<KFParticle> >* StKFParticleInterface::GetPrimaryTopoCandidates() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetPrimaryTopoCandidates();        } // Get primary particles with the topologigal constraint
const std::vector< std::vector<KFParticle> >* StKFParticleInterface::GetPrimaryTopoMassCandidates() const {return fKFParticleTopoReconstructor->GetKFParticleFinder()->GetPrimaryTopoMassCandidates();} // Get primary particles with the topologigal and mass constraint

KFParticleFinder*  StKFParticleInterface::GetKFParticleFinder() { return fKFParticleTopoReconstructor->GetKFParticleFinder(); }
void StKFParticleInterface::SetMaxDistanceBetweenParticlesCut(float cut) { GetKFParticleFinder()->SetMaxDistanceBetweenParticlesCut(cut); }
void StKFParticleInterface::SetLCut(float cut)                           { GetKFParticleFinder()->SetLCut(cut); }
void StKFParticleInterface::SetChiPrimaryCut2D(float cut)                { GetKFParticleFinder()->SetChiPrimaryCut2D(cut); }
void StKFParticleInterface::SetChi2Cut2D(float cut)                      { GetKFParticleFinder()->SetChi2Cut2D(cut); }
void StKFParticleInterface::SetLdLCut2D(float cut)                       { GetKFParticleFinder()->SetLdLCut2D(cut); }
void StKFParticleInterface::SetLdLCutXiOmega(float cut)                  { GetKFParticleFinder()->SetLdLCutXiOmega(cut); }
void StKFParticleInterface::SetChi2TopoCutXiOmega(float cut)             { GetKFParticleFinder()->SetChi2TopoCutXiOmega(cut); }
void StKFParticleInterface::SetChi2CutXiOmega(float cut)                 { GetKFParticleFinder()->SetChi2CutXiOmega(cut); }
void StKFParticleInterface::SetChi2TopoCutResonances(float cut)          { GetKFParticleFinder()->SetChi2TopoCutResonances(cut); }
void StKFParticleInterface::SetChi2CutResonances(float cut)              { GetKFParticleFinder()->SetChi2CutResonances(cut); }
void StKFParticleInterface::SetPtCutLMVM(float cut)                      { GetKFParticleFinder()->SetPtCutLMVM(cut); }
void StKFParticleInterface::SetPCutLMVM(float cut)                       { GetKFParticleFinder()->SetPCutLMVM(cut); }
void StKFParticleInterface::SetPtCutJPsi(float cut)                      { GetKFParticleFinder()->SetPtCutJPsi(cut); }
void StKFParticleInterface::SetPtCutCharm(float cut)                     { GetKFParticleFinder()->SetPtCutCharm(cut); }
void StKFParticleInterface::SetChiPrimaryCutCharm(float cut)             { GetKFParticleFinder()->SetChiPrimaryCutCharm(cut); }
void StKFParticleInterface::SetLdLCutCharmManybodyDecays(float cut)      { GetKFParticleFinder()->SetLdLCutCharmManybodyDecays(cut); }
void StKFParticleInterface::SetChi2TopoCutCharmManybodyDecays(float cut) { GetKFParticleFinder()->SetChi2TopoCutCharmManybodyDecays(cut); }
void StKFParticleInterface::SetChi2CutCharmManybodyDecays(float cut)     { GetKFParticleFinder()->SetChi2CutCharmManybodyDecays(cut); }
void StKFParticleInterface::SetLdLCutCharm2D(float cut)                  { GetKFParticleFinder()->SetLdLCutCharm2D(cut); }
void StKFParticleInterface::SetChi2TopoCutCharm2D(float cut)             { GetKFParticleFinder()->SetChi2TopoCutCharm2D(cut); }
void StKFParticleInterface::SetChi2CutCharm2D(float cut)                 { GetKFParticleFinder()->SetChi2CutCharm2D(cut); }
  
double StKFParticleInterface::InversedChi2Prob(double p, int ndf) const
{
  double epsilon = 1.e-14;
  double chi2Left = 0.f;
  double chi2Right = 10000.f;
  
  double probLeft = p - TMath::Prob(chi2Left, ndf);
  
  double chi2Centr = (chi2Left+chi2Right)/2.f;
  double probCentr = p - TMath::Prob( chi2Centr, ndf);
  
  while( TMath::Abs(chi2Right-chi2Centr)/chi2Centr > epsilon )
  {
    if(probCentr * probLeft > 0.f)
    {
      chi2Left = chi2Centr;
      probLeft = probCentr;
    }
    else
    {
      chi2Right = chi2Centr;
    }
    
    chi2Centr = (chi2Left+chi2Right)/2.f;
    probCentr = p - TMath::Prob( chi2Centr, ndf);
  }
  
  return chi2Centr;
}

void StKFParticleInterface::SetPrimaryProbCut(float prob)
{ 
  fKFParticleTopoReconstructor->SetChi2PrimaryCut( InversedChi2Prob(prob, 2) );
}

void StKFParticleInterface::CollectTrackHistograms()
{
  TDirectory *dirs[7] = {0};
  dirs[0] = TDirectory::CurrentDirectory(); assert(dirs[0]);
  dirs[0]->cd();
  if (! dirs[0]->GetDirectory("Tracks")) {
    dirs[0]->mkdir("Tracks");
  }
  dirs[1] = dirs[0]->GetDirectory("Tracks"); assert(dirs[1]);
  dirs[1]->cd();
  
  fTrackHistograms2D[0] = (TH2F *)   dirs[1]->Get("hdEdX");
  if (! fTrackHistograms2D[0]) fTrackHistograms2D[0] = new TH2F("hdEdX", "hdEdX", 2000, 0, 10, 3000, 0, 30);

  fTrackHistograms2D[1] = (TH2F *)   dirs[1]->Get("hdEdXPos");
  if (! fTrackHistograms2D[1]) fTrackHistograms2D[1] = new TH2F("hdEdXPos", "hdEdXPos", 2000, 0, 10, 3000, 0, 30);
  
  fTrackHistograms2D[2] = (TH2F *)   dirs[1]->Get("hdEdXNeg");
  if (! fTrackHistograms2D[2]) fTrackHistograms2D[2] = new TH2F("hdEdXNeg", "hdEdXNeg", 2000, 0, 10, 3000, 0, 30);
  
  fTrackHistograms2D[3] = (TH2F *)   dirs[1]->Get("hdEdXwithToF");
  if (! fTrackHistograms2D[3]) fTrackHistograms2D[3] = new TH2F("hdEdXwithToF", "hdEdXwithToF", 2000, 0, 10, 3000, 0, 30);
  
  fTrackHistograms2D[4] = (TH2F *)   dirs[1]->Get("hTofPID");
  if (! fTrackHistograms2D[4]) fTrackHistograms2D[4] = new TH2F("hTofPID", "hTofPID", 300, 0, 15, 1100, -1, 10);

  fTrackHistograms[0] = (TH1F *)   dirs[1]->Get("hNHFTHits");
  if (! fTrackHistograms[0]) fTrackHistograms[0] = new TH1F("hNHFTHits", "hNHFTHits",11, -0.5, 10.5);
  
  fTrackHistograms[1] = (TH1F *)   dirs[1]->Get("hPVError");
  if (! fTrackHistograms[1]) fTrackHistograms[1] = new TH1F("hPVError", "hPVError", 10000, 0, 1);

  fTrackHistograms2D[5] = (TH2F *)   dirs[1]->Get("hPVErrorVsNTracks");
  if (! fTrackHistograms2D[5]) fTrackHistograms2D[5] = new TH2F("hPVErrorVsNTracks", "hPVErrorVsNTracks", 5000, 0, 5000, 5000, 0, 0.5);

  fTrackHistograms2D[6] = (TH2F *)   dirs[1]->Get("hPVErrorVsNPVTracks");
  if (! fTrackHistograms2D[6]) fTrackHistograms2D[6] = new TH2F("hPVErrorVsNPVTracks", "hPVErrorVsNPVTracks", 5000, 0, 5000, 5000, 0, 0.5);
  
  dirs[0]->cd();
  
  fCollectTrackHistograms = true;
}

void StKFParticleInterface::CollectPIDHistograms()
{
  TDirectory *dirs[7] = {0};
  dirs[0] = TDirectory::CurrentDirectory(); assert(dirs[0]);
  dirs[0]->cd();
  if (! dirs[0]->GetDirectory("Tracks")) {
    dirs[0]->mkdir("Tracks");
  }
  dirs[1] = dirs[0]->GetDirectory("Tracks"); assert(dirs[1]);
  dirs[1]->cd();
  
  int pdgTrackHisto[NTrackHistoFolders] = { 11, -11, 13, -13, 211, -211, 321, -321, 2212, -2212, 
                                            1000010020, -1000010020, 1000010030, -1000010030, 1000020030, -1000020030, 1000020040, -1000020040 };
  TString trackFolderName[NTrackHistoFolders] = {"e-", "e+", "mu-", "mu+", "pi+", "pi-", "K+", "K-", "p", "p-", "d", "d-", "t", "t-", "He3", "He3-", "He4", "He4-"};
                    
  for(int iTrackHisto=0; iTrackHisto<NTrackHistoFolders; iTrackHisto++)
  {
    if (!dirs[1]->GetDirectory(trackFolderName[iTrackHisto].Data()))
      dirs[1]->mkdir(trackFolderName[iTrackHisto].Data());
    
    dirs[2] = dirs[1]->GetDirectory(trackFolderName[iTrackHisto].Data()); assert(dirs[2]);
    dirs[2]->cd();
    
    fTrackPdgToHistoIndex[ pdgTrackHisto[iTrackHisto] ] = iTrackHisto;
    
    fHistodEdXTracks[iTrackHisto] = (TH2F *)   dirs[2]->Get("hdEdX");
    if (! fHistodEdXTracks[iTrackHisto]) fHistodEdXTracks[iTrackHisto] = new TH2F("hdEdX", "hdEdX", 2000, 0, 10, 3000, 0, 30);

    fHistodEdXwithToFTracks[iTrackHisto] = (TH2F *)   dirs[2]->Get("hdEdXwithToF");
    if (! fHistodEdXwithToFTracks[iTrackHisto]) fHistodEdXwithToFTracks[iTrackHisto] = new TH2F("hdEdXwithToF", "hdEdXwithToF", 2000, 0, 10, 3000, 0, 30);
  
    fHistoTofPIDTracks[iTrackHisto] = (TH2F *)   dirs[2]->Get("hTofPID");
    if (! fHistoTofPIDTracks[iTrackHisto]) fHistoTofPIDTracks[iTrackHisto] = new TH2F("hTofPID", "hTofPID", 300, 0, 15, 1100, -1, 10);
  
    fHistoMomentumTracks[iTrackHisto] = (TH1F *)   dirs[2]->Get("hMomentum");
    if (! fHistoMomentumTracks[iTrackHisto]) fHistoMomentumTracks[iTrackHisto] = new TH1F("hMomentum", "hMomentum", 1000, 0, 10);
    
    fHistodEdXPull[iTrackHisto] = (TH2F *)   dirs[2]->Get("hdEdXPull");
    if (! fHistodEdXPull[iTrackHisto]) fHistodEdXPull[iTrackHisto] = new TH2F("hdEdXPull", "hdEdXPull", 2000, 0, 10, 600, -30, 30);
    
    fHistodEdXZ[iTrackHisto] = (TH2F *)   dirs[2]->Get("hdEdXZ");
    if (! fHistodEdXZ[iTrackHisto]) fHistodEdXZ[iTrackHisto] = new TH2F("hdEdXZ", "hdEdXZ", 2000, -5, 5, 280, -1, 6);
    
    dirs[1]->cd();
  }
  
  dirs[0]->cd();
  
  fCollectPIDHistograms = true;
}

bool StKFParticleInterface::IsGoodPV(const KFVertex& pv)
{
//  bool isGoodPV = (pv.X() > -0.3) && (pv.X() < -0.1) &&
//                  (pv.Y() > -0.27) && (pv.Y() < -0.1);
  bool isGoodPV = (pv.Z() > -70) && (pv.Z() < 70) && 
		  (sqrt(pv.X()*pv.X() + pv.Y()*pv.Y()) < 2);                 
  return isGoodPV;
}

bool StKFParticleInterface::GetTrack(const StDcaGeometry& dcaG, KFPTrack& track, int q, int index)
{
  Double_t xyzp[6], CovXyzp[21];
  dcaG.GetXYZ(xyzp,CovXyzp);
  
  bool goodTrack=1;
  for(int iPar=0; iPar<6; iPar++)
    goodTrack = goodTrack && finite(xyzp[iPar]);
  for(int iC=0; iC<21; iC++)
    goodTrack = goodTrack && finite(CovXyzp[iC]);
  goodTrack &= goodTrack && CovXyzp[0]  >=0.f && CovXyzp[0]  < 100.f;
  goodTrack &= goodTrack && CovXyzp[2]  >=0.f && CovXyzp[2]  < 100.f;
  goodTrack &= goodTrack && CovXyzp[5]  >=0.f && CovXyzp[5]  < 100.f;
  goodTrack &= goodTrack && CovXyzp[9]  >=0.f && CovXyzp[9]  < 1.f;
  goodTrack &= goodTrack && CovXyzp[14] >=0.f && CovXyzp[14] < 1.f;
  goodTrack &= goodTrack && CovXyzp[20] >=0.f && CovXyzp[20] < 1.f;
  if(!goodTrack) return false;
  
  track.SetParameters(xyzp);
  track.SetCovarianceMatrix(CovXyzp);
  track.SetNDF(1);
  //    track.SetChi2(GlobalTracks_mChiSqXY[k]);
  track.SetID(index);

  track.SetCharge(q);
  return true;
}

std::vector<int> StKFParticleInterface::GetTofPID(double m2, double p, int q, const int trackId)
{
  static const int order = 4;
  static const double parMean[6][order+1] = { { 0.02283190,-0.01482910, 0.01883130,-0.01824250, 0.00409811  }, //pi+
                                              { 0.24842500,-0.00699781,-0.00991387, 0.01327170,-0.00694824  }, //K+
                                              { 0.863211  , 0.0264171 ,-0.0230833 , 0.00239637, 0.000262309 }, //p
                                              { 0.0224095 ,-0.0123235 , 0.0145216 ,-0.0149944 , 0.00325952  }, //pi-
                                              { 0.250696  ,-0.0151308 , 0.00437457, 0.00516669,-0.00529184  }, //K-
                                              { 0.886912  ,-0.0298543 , 0.0449904 ,-0.0286879 , 0.00541963  }};//p-
  static const double parSigma[6][order+1] = { { 0.0112498,-0.0400571, 0.0733615,-0.0316505, 0.00629469 }, //pi+
                                               { 0.0154830,-0.0396312, 0.0719647,-0.0290683, 0.00637164 }, //K+
                                               { 0.114465 ,-0.287213 , 0.356536 ,-0.169257 , 0.0299844  }, //p
                                               { 0.0111682,-0.0394877, 0.0718342,-0.0302914, 0.00587317 }, //pi-
                                               { 0.0157322,-0.0402606, 0.0716639,-0.0272101, 0.00564467 }, //K-
                                               { 0.0899438,-0.211922 , 0.273122 ,-0.129597 , 0.0231844  }};//p-
  double pMax = 2.;
  double nSigmas[3];
  for(int iHypothesys = 0; iHypothesys<3; iHypothesys++)
  {
    double x = p;
    if(x>=pMax) x = pMax;
    
    int iSet = iHypothesys;
    if(q<0)
      iSet += 3;
    double mean = 0;
    for(int iTerm=0; iTerm<=order; iTerm++)
      mean += parMean[iSet][iTerm]*TMath::Power(x,iTerm);  
    
    double sigma = 0;
    for(int iTerm=0; iTerm<=order; iTerm++)
      sigma += parSigma[iSet][iTerm]*TMath::Power(x,iTerm);  
    
    nSigmas[iHypothesys] = fabs((m2 - mean)/sigma);
    fTrackPidTof[iHypothesys][trackId] = nSigmas[iHypothesys];
  }
  
  double minNSigma = nSigmas[0];
  int minHypothesis = 0;
  for(int iHypothesys=1; iHypothesys<3; iHypothesys++)
  {
    if(minNSigma > nSigmas[iHypothesys]) 
    {
      minNSigma = nSigmas[iHypothesys];
      minHypothesis = iHypothesys;
    }
  }

  int pdgHypothesis[3] = {211, 321, 2212};
  vector<int> tofPID;
  
  if(fStrictTofPID)
  {
    if(minNSigma < 3)
      tofPID.push_back(pdgHypothesis[minHypothesis]*q);
  }
  else
  {    
    for(int iHypothesys=0; iHypothesys<3; iHypothesys++)
      if(nSigmas[iHypothesys] < 3)
        tofPID.push_back(pdgHypothesis[iHypothesys]*q);
  }
  
  return tofPID;
}

std::vector<int> StKFParticleInterface::GetPID(double m2, double p, int q, double dEdX, double dEdXPull[7], bool isTofm2, const int trackId)
{
  vector<int> ToFPDG;
  if(isTofm2)
    ToFPDG = GetTofPID(m2, p, q, trackId);
  
  for(int iPdg=0; iPdg<3; iPdg++)
    fTrackPidTpc[iPdg][trackId] = dEdXPull[iPdg];
  
  vector<int> dEdXPDG;
  float nSigmaCut = 3.f; //TODO
  
  bool checkKTof = false;
  if(fCleanKaonsWitTof)
    checkKTof = (p > 0.5) && (p < 2.);
  bool checkKHasTof = 0;
  for(unsigned int iTofPDG=0; iTofPDG<ToFPDG.size(); iTofPDG++)
    if(abs(ToFPDG[iTofPDG]) == 321)
      checkKHasTof = 1;

  if(dEdXPull[0] < nSigmaCut)                                           dEdXPDG.push_back(211*q);  
  if(dEdXPull[1] < 2.f && ((checkKTof && checkKHasTof) || !checkKTof) ) dEdXPDG.push_back(321*q);
  if(dEdXPull[2] < nSigmaCut)                                           dEdXPDG.push_back(2212*q); 
      
  vector<int> totalPDG;
  if(!isTofm2)
    totalPDG = dEdXPDG;
  else
  {
    for(unsigned int iPDG=0; iPDG<dEdXPDG.size(); iPDG++)
      for(unsigned int iTofPDG=0; iTofPDG<ToFPDG.size(); iTofPDG++)
        if(dEdXPDG[iPDG] == ToFPDG[iTofPDG])
          totalPDG.push_back(ToFPDG[iTofPDG]);        
  }
  
  {
//    if(dEdXPull[5] < nSigmaCut) totalPDG.push_back(1000020030*q);
//    if(dEdXPull[6] < nSigmaCut) { totalPDG.push_back(1000020040*q); }
    
    if(dEdXPull[3] < nSigmaCut && dEdXPull[2] > nSigmaCut) 
      if( isTofm2 && (m2 > 2 && m2<6) ) //if( !isTofm2 || (isTofm2 && (m2 > 2 && m2<6)) )
        totalPDG.push_back(1000010020*q); 
    if(dEdXPull[4] < nSigmaCut && dEdXPull[3] > nSigmaCut) 
      if( isTofm2 && (m2 > 5) ) //if( !isTofm2 || (isTofm2 && (m2 > 5)) )
        totalPDG.push_back(1000010030*q);
   
    if(p>0.6 && p<1.8)
    {
      double lowerParameters[4] = {2.48051e+01,-1.40961e+00, 2.54327e-01, 2.07450e-01};
      double lowerHe3Bound = lowerParameters[0]*TMath::Power(p, lowerParameters[1] + lowerParameters[2]*log(p) + lowerParameters[3]*log(p)*log(p));
      
      double upperParameters[4] = {3.36169e+01,-1.33686e+00, 2.82856e-01, 1.91841e-01};
      double upperHe3Bound = upperParameters[0]*TMath::Power(p, upperParameters[1] + upperParameters[2]*log(p) + upperParameters[3]*log(p)*log(p));
      
      if(dEdX > lowerHe3Bound && dEdX < upperHe3Bound && dEdX > 9.) 
        if( !isTofm2 || (isTofm2 && (m2>1.) && (m2<5.) ) )
          totalPDG.push_back(1000020030*q);
      if(dEdX > upperHe3Bound && dEdX > 9.)               
        if( !isTofm2 || (isTofm2 && (m2>2.5) && (m2<6.) ) )
          totalPDG.push_back(1000020040*q);
    }
    else if(p>=1.8 && dEdX > 9.)
    {
      if(dEdXPull[5] < nSigmaCut && dEdXPull[4] > nSigmaCut) 
        if( !isTofm2 || (isTofm2 && (m2>1.) && (m2<5.) ) )
          totalPDG.push_back(1000020030*q);
      if(dEdXPull[6] < nSigmaCut && dEdXPull[5] > nSigmaCut) 
        if( !isTofm2 || (isTofm2 && (m2>2.5) && (m2<6.) ) )
          totalPDG.push_back(1000020040*q);
    }
  }
    
  if(totalPDG.size() == 0)
    totalPDG.push_back(-1);
  
  return totalPDG;
}

void StKFParticleInterface::AddTrackToParticleList(const KFPTrack& track, int nHftHitsInTrack, int index, const std::vector<int>& totalPDG, KFVertex& pv, 
  std::vector<int>& primaryTrackList, std::vector<int>& nHftHits, std::vector<int>& particlesPdg, std::vector<KFParticle>& particles, int& nPartSaved)
{
  for(unsigned int iPDG=0; iPDG<totalPDG.size(); iPDG++)
  {
    if( fTriggerMode && (nHftHitsInTrack < 3) ) continue;

    int pdg = totalPDG[iPDG];
    
    KFPTrack trackPDG = track;
    if(abs(pdg) == 1000020030 || abs(pdg) == 1000020040)
    {
      trackPDG.SetCharge( trackPDG.Charge()*2.f );
      trackPDG.SetPx( trackPDG.GetPx()*2.f );
      trackPDG.SetPy( trackPDG.GetPy()*2.f );
      trackPDG.SetPz( trackPDG.GetPz()*2.f );

      const int index2[9] = { 6,7,8, 10,11,12, 15,16,17 };
      for(int iIndex=0; iIndex<9; iIndex++)
      {
        const int iC = index2[iIndex];
//        cout<<iC<<endl;
        trackPDG.SetCovariance( iC, trackPDG.GetCovariance(iC)*2.f );
      }
      const int index4[6] = { 9, 13,14, 18,19,20 };
      for(int iIndex=0; iIndex<6; iIndex++)
      {
        const int iC = index4[iIndex];
        trackPDG.SetCovariance( iC, trackPDG.GetCovariance(iC)*4.f );
      }
    }
    
    nHftHits[nPartSaved] = nHftHitsInTrack;
   

    KFParticle particle(trackPDG, pdg);
    float chiPrim = particle.GetDeviationFromVertex(pv);
    if(chiPrim < fChiPrimaryCut)
    {
      if(fTriggerMode) continue;
      primaryTrackList.push_back(nPartSaved);
    }
    if(fTriggerMode && chiPrim > fChiPrimaryMaxCut) continue;
//     if( chiPrim > 1.e6 ) continue;
//
//    cout<<"index:"<<index<<" "<<totalPDG.size()<<endl;     
    particle.SetId(index);
    particles[nPartSaved] = particle;

    particlesPdg[nPartSaved] = pdg;
//    cout<<"nSaved:"<<nPartSaved<<endl;
    nPartSaved++;
  }
}

void StKFParticleInterface::FillPIDHistograms(StPicoTrack *gTrack, const std::vector<int>& pdgVector, const bool isTofm2, float m2tof)
{
  float momentum = gTrack->gPtot();
  for(unsigned int iPdg = 0; iPdg<pdgVector.size(); iPdg++)
  {
    int pdg = pdgVector[iPdg];
    const int iTrackHisto = fTrackPdgToHistoIndex[pdg];
    if( ! (iTrackHisto < 0 || iTrackHisto >= NTrackHistoFolders) )
    {
      fHistoMomentumTracks[iTrackHisto] -> Fill(momentum);
      fHistodEdXTracks[iTrackHisto] -> Fill(momentum, gTrack->dEdx());
      if(isTofm2)
      {
        fHistodEdXwithToFTracks[iTrackHisto] -> Fill(momentum, gTrack->dEdx());
        fHistoTofPIDTracks[iTrackHisto] -> Fill(momentum, m2tof);
        
        if(abs(pdg)==211)
        {
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(0.139570, fdEdXMode, 1));
          float betaGamma = TMath::Log10(momentum/0.139570);
          float z = gTrack->dEdxPull(0.139570, fdEdXMode, 1)*gTrack->dEdxError();
          fHistodEdXZ[iTrackHisto]->Fill(betaGamma, z);
          
          betaGamma = TMath::Log10(momentum/5.485799e-4);
          z = gTrack->nSigmaElectron()*gTrack->dEdxError();
          fHistodEdXZ[0]->Fill(betaGamma, z);
        }
        if(abs(pdg)==321)
        {
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(0.493677, fdEdXMode, 1));
          float betaGamma = TMath::Log10(momentum/0.493677);
          float z = gTrack->dEdxPull(0.493677, fdEdXMode, 1)*gTrack->dEdxError();
          fHistodEdXZ[iTrackHisto]->Fill(betaGamma, z);
        }
        if(abs(pdg)==2212)
        {
          fHistodEdXPull[iTrackHisto] -> Fill(momentum, gTrack->dEdxPull(0.938272, fdEdXMode, 1));
          float betaGamma = TMath::Log10(momentum/0.938272);
          float z = gTrack->dEdxPull(0.938272, fdEdXMode, 1)*gTrack->dEdxError();
          fHistodEdXZ[iTrackHisto]->Fill(betaGamma, z);
        }
      }
    }
  }
}

bool StKFParticleInterface::OpenCharmTrigger() 
{
  bool triggerDMesons = false;
  if(fKFParticleTopoReconstructor->NPrimaryVertices() == 0) return false;
    
  for(unsigned int iParticle=0; iParticle<GetParticles().size(); iParticle++)
  {
    KFParticle particle = GetParticles()[iParticle];
    
    if( abs(particle.GetPDG()) == 421 ||
        abs(particle.GetPDG()) == 429 || 
        abs(particle.GetPDG()) == 420 || 
        abs(particle.GetPDG()) == 411 || 
        abs(particle.GetPDG()) == 431 || 
        abs(particle.GetPDG()) == 4122 ||
        abs(particle.GetPDG()) == 426 )
    {
      KFParticleSIMD tempSIMDPart(particle);
      float_v l,dl;
      KFParticleSIMD pv(fKFParticleTopoReconstructor->GetPrimVertex());
      tempSIMDPart.GetDistanceToVertexLine(pv, l, dl);
      
      if(abs(particle.GetPDG()) == 411)
        triggerDMesons = (l[0] < 0.4);
      else    
        triggerDMesons = (l[0] < 0.2);
    }
  }
  
  return triggerDMesons;
}

void StKFParticleInterface::OpenCharmTriggerCompression(int nTracksTriggered, int nTracksInEvent, bool triggerDMesons) 
{
  static int nTriggeredEvents = 0;
  static int nTracksInEventTriggered = 0;
  static int nTracksInEventTotal = 0;
  static int nEvents = 0;
  nEvents++;
  nTracksInEventTotal += nTracksInEvent;
  if(triggerDMesons)
  {
    nTriggeredEvents++;
    nTracksInEventTriggered += nTracksTriggered;
    std::cout << "N Events " << nEvents << "    N triggered events " << nTriggeredEvents << "    ratio " << (double(nEvents)/double(nTriggeredEvents)) << std::endl;
    std::cout << "N Tracks " << nTracksInEventTotal << "    N triggered events " << nTracksInEventTriggered << "    ratio " << (double(nTracksInEventTotal)/double(nTracksInEventTriggered)) << std::endl;
  }
}

void StKFParticleInterface::ResizeTrackPidVectors(const int nTracks)
{
  for(int iHypothesis=0; iHypothesis<3; iHypothesis++)
  {
    fTrackPidTof[iHypothesis].clear();
    fTrackPidTof[iHypothesis].resize(nTracks, -1);
    
    fTrackPidTpc[iHypothesis].clear();
    fTrackPidTpc[iHypothesis].resize(nTracks, -1);
  }
}

bool StKFParticleInterface::ProcessEvent(StPicoDst* picoDst, std::vector<int>& triggeredTracks)
{
  triggeredTracks.resize(0);
  
  //read PV from pico Event
  KFVertex primaryVertex;
  vector<int> primaryTrackList;
    
  StPicoEvent* picoEvent = picoDst->event();
  if(!picoEvent) return 0;
  //cout<<"Event No:"<<picoEvent->eventId()<<endl;

  //b                                

  
  const TVector3 picoPV = picoEvent->primaryVertex();
  const TVector3 picoPVError = picoEvent->primaryVertexError();
  
  KFPVertex primVtx_tmp;
  primVtx_tmp.SetXYZ(picoPV.x(), picoPV.y(), picoPV.z());
  double dx = picoPVError.x();
  double dy = picoPVError.y();
  double dz = picoPVError.z();
//cout<<"dx:"<<dx<<" "<<dy<<" "<<dz<<endl;
  primVtx_tmp.SetCovarianceMatrix( dx*dx, 0, dy*dy, 0, 0, dz*dz );
  primaryVertex = KFVertex(primVtx_tmp);
   if(!IsGoodPV(primaryVertex)) return 0;
//   if(IsGoodPV(primaryVertex)){cout<<"Good PV!"<<endl; }

  
  Int_t nGlobalTracks = picoDst->numberOfTracks( );
  
  fParticles.resize(nGlobalTracks*7);
  fNHftHits.resize(nGlobalTracks*7);
  fParticlesPdg.resize(nGlobalTracks*7);
  int nPartSaved = 0;
  int nUsedTracks = 0;
  
  for (Int_t iTrack = 0; iTrack < nGlobalTracks; iTrack++) 
  {
    StPicoTrack *gTrack = picoDst->track(iTrack);
    if (! gTrack)            continue;
    if (! gTrack->charge())  continue;
    if (  gTrack->nHitsFit() < 15) continue;
    if (  gTrack->dEdxError() < 0.04 || gTrack->dEdxError() > 0.12 ) continue;
    const int index = gTrack->id();
    
    int nHftHitsInTrack = 0;
    if(gTrack->hasPxl1Hit()) nHftHitsInTrack++;
    if(gTrack->hasPxl2Hit()) nHftHitsInTrack++;
    if(gTrack->hasIstHit()) nHftHitsInTrack++;
//       if(gTrack->hasSstHit()) nHftHitsInTrack++;
    if(fCollectTrackHistograms) fTrackHistograms[0]->Fill(nHftHitsInTrack);
    
    if(fUseHFTTracksOnly && nHftHitsInTrack < 3) continue;
    
    StPicoTrackCovMatrix *cov = picoDst->trackCovMatrix(iTrack);
    const StDcaGeometry dcaG = cov->dcaGeometry();
    Int_t q = 1; if (gTrack->charge() < 0) q = -1;
    KFPTrack track;
    if( !GetTrack(dcaG, track, q, index) ) continue;
    
    if(fCollectTrackHistograms)
    {
      fTrackHistograms2D[0]->Fill(track.GetP(), gTrack->dEdx());
      if(q>0) fTrackHistograms2D[1]->Fill(track.GetP(), gTrack->dEdx());
      else    fTrackHistograms2D[2]->Fill(track.GetP(), gTrack->dEdx());  
    }
    
    double m2tof = -1.e6;
    bool isTofm2 = false;
    if(gTrack->bTofPidTraitsIndex() > 0)
    {
      const StPicoBTofPidTraits* btofPid = picoDst->btofPidTraits(gTrack->bTofPidTraitsIndex());
      double betaTof2 = btofPid->btofBeta() * btofPid->btofBeta();
      if(fabs(betaTof2) > 1.e-6)
      {
        m2tof = track.GetP()*track.GetP()*(1./betaTof2 - 1.);
        isTofm2 = true;
      }
      if(fCollectTrackHistograms)
      {
        fTrackHistograms2D[3]->Fill(track.GetP(), gTrack->dEdx());
        fTrackHistograms2D[4]->Fill(track.GetP(), m2tof);
      }
    }

    double dEdXPull[7] = { fabs(gTrack->dEdxPull(0.139570, fdEdXMode, 1)),   //0 - pi
                           fabs(gTrack->dEdxPull(0.493677, fdEdXMode, 1)),   //1 - K
                           fabs(gTrack->dEdxPull(0.938272, fdEdXMode, 1)),   //2 - p
                           fabs(gTrack->dEdxPull(1.876124, fdEdXMode, 1)),   //3 - d
                           fabs(gTrack->dEdxPull(2.809432, fdEdXMode, 1)),   //4 - t
                           fabs(gTrack->dEdxPull(2.809413, fdEdXMode, 2)),   //5 - He3
                           fabs(gTrack->dEdxPull(3.728400, fdEdXMode, 2))};  //6 - He4
    
    vector<int> totalPDG = GetPID(m2tof, track.GetP(), q, gTrack->dEdx(), dEdXPull, isTofm2, index);
    
    int nPartSaved0 = nPartSaved;
    AddTrackToParticleList(track, nHftHitsInTrack, index, totalPDG, primaryVertex, primaryTrackList, fNHftHits, fParticlesPdg, fParticles, nPartSaved); 
    
    if(nPartSaved > nPartSaved0) 
      triggeredTracks.push_back(iTrack);
    
    //fill PID histograms if they are created
    if(fCollectPIDHistograms) FillPIDHistograms(gTrack, totalPDG, isTofm2, m2tof);
    
    nUsedTracks++;
  }
  
  fParticles.resize(nPartSaved);
  fParticlesPdg.resize(nPartSaved);
  fNHftHits.resize(nPartSaved);
  
  if( fCleanLowPVTrackEvents && ( 10*primaryTrackList.size() < (nUsedTracks - primaryTrackList.size()) ) ) return 0;
  
  const Double_t field = picoEvent->bField();  
  SetField(field);

  CleanPV();
  InitParticles();

  //read PV
  AddPV(primaryVertex, primaryTrackList);
  if(fCollectTrackHistograms)
  {
    fTrackHistograms[1]->Fill(sqrt(dx*dx + dy*dy));
    fTrackHistograms2D[5]->Fill( nPartSaved, sqrt(dx*dx + dy*dy) );
    fTrackHistograms2D[6]->Fill( primaryTrackList.size(), sqrt(dx*dx + dy*dy) );
  }  
  //reconstruct short-lived particles
  ReconstructParticles();
  
  return 1;
}

bool StKFParticleInterface::ProcessEvent(StMuDst* muDst, vector<KFMCTrack>& mcTracks, vector<int>& mcIndices, bool processSignal)
{  
  mcTracks.resize(muDst->numberOfMcTracks());

//  cout<<"WHERE!!!!!!-------------"<<endl;
//  cout<<"numberOfMcTracks:"<<muDst->numberOfMcTracks()<<endl;

  for (unsigned int iMCTrack=0; iMCTrack<muDst->numberOfMcTracks(); iMCTrack++) 
  {
    StMuMcTrack *mcTrack = muDst->MCtrack(iMCTrack);
    if (! mcTrack) continue;    
    KFMCTrack &mcTrackKF = mcTracks[iMCTrack];
    mcTrack->FillKFMCTrack(mcTrackKF);
    mcTrackKF.SetNMCPixelPoints(mcTrack->No_ist_hit() + mcTrack->No_ssd_hit() + mcTrack->No_pix_hit());
  }
  
  //read PV
  KFVertex primaryVertex;
  vector<int> primaryTrackList;

  float bestRank=-1000000;
  int bestPV=-1;
  double dx = 0., dy = 0., dz = 0.;
  for(unsigned int iPV=0; iPV<muDst->numberOfPrimaryVertices(); iPV++) 
  {
    StMuPrimaryVertex *Vtx = muDst->primaryVertex(iPV);
    if(!Vtx) continue;
    if (bestRank < Vtx->ranking()) {
      bestRank = Vtx->ranking();
      bestPV = iPV;
    }
    else continue;
    
    //convert StMuPrimaryVertex to KFVertex
    KFPVertex kfVertex;
    kfVertex.SetXYZ(Vtx->position().x(), Vtx->position().y(), Vtx->position().z());
    dx = Vtx->posError().x();
    dy = Vtx->posError().y();
    dz = Vtx->posError().z();

//    cout<<"dx:"<<dx<<" "<<dy<<" "<<dz<<endl;

  //TOBEREVERTED
  /*
    TF1 *fa0x = new TF1("fa0x","gaus(0)",0.01,0.02);
    fa0x -> SetParameters(1,0.01476,0.00135);
    TF1 *fa0y = new TF1("fa0y","gaus(0)",0.01,0.02);
    fa0y -> SetParameters(1,0.01489,0.001415);
    TF1 *fa0z = new TF1("fa0z","gaus(0)",0.01,0.02);
    fa0z -> SetParameters(1,0.0119,0.0009366); 

    TF1 *fa1x = new TF1("fa1x","gaus(0)",0.01,0.04);
    fa1x -> SetParameters(1,0.02193,0.004218);
    TF1 *fa1y = new TF1("fa1y","gaus(0)",0.01,0.04);
    fa1y -> SetParameters(1,0.02215,0.004292);
    TF1 *fa1z = new TF1("fa1z","gaus(0)",0.01,0.04);
    fa1z -> SetParameters(1,0.01789,0.003246);

    TF1 *fa3 = new TF1("fa3","1/x/x/x",0.03,0.15); 
    


    int drefmult = muDst->event()->refMult();
    if(drefmult<=69){
//    dx = 0.05767;//0.03435
//    dy = 0.05852;//0.03624
//    dz = 0.04704;//0.02577
    dx = fa3->GetRandom(0.028,0.15);
    dy = fa3->GetRandom(0.028,0.15);
    dz = fa3->GetRandom(0.028,0.15);
    
    }else if(drefmult>69 && drefmult<=226){
//    dx = 0.02193;//0.004218
//    dy = 0.02215;//0.004292
//    dz = 0.01789;//0.003246
    dx = fa1x->GetRandom(0.01,0.04);
    dy = fa1y->GetRandom(0.01,0.04);
    dz = fa1z->GetRandom(0.01,0.04);
    
    }else{
//    dx = 0.01476;//0.00135
//    dy = 0.01489;//0.001415
//    dz = 0.0119;//0.0009366
    dx = fa0x->GetRandom(0.01,0.02);    
    dy = fa0y->GetRandom(0.01,0.02);
    dz = fa0z->GetRandom(0.01,0.02); 
    }
*/
//    cout<<"new dx:"<<dx<<" "<<dy<<" "<<dz<<" "<<drefmult<<endl;
	
     StRefMultCorr* refmultCorrUtil = CentralityMaker::instance()->getRefMultCorr() ;
     refmultCorrUtil -> init(muDst->event()->runId());
     int dcent9 = refmultCorrUtil -> getCentralityBin9();
     //cout<<"dcent9:"<<dcent9<<endl;
     
     TF1 *f0gaus = new TF1("f0gaus","gaus(0)",0,0.25);
     TF1 *f0land = new TF1("f0land","[2]*TMath::Landau(x,[0],[1],0)",0,0.25);

     if(dcent9==8){
     f0gaus->SetParameters(1025.05, 0.0138299, -0.000607724);
     dx = f0gaus->GetRandom(0.001,0.25);
     f0gaus->SetParameters(983.332, 0.0139099, -0.000635113);
     dy = f0gaus->GetRandom(0.001,0.25);
     f0gaus->SetParameters(1275.22, 0.0111875, -0.000494996);
     dz = f0gaus->GetRandom(0.001,0.25);
     }else if(dcent9==7){
f0gaus->SetParameters(884.35, 0.0150516, -0.000669005);
dx = f0gaus->GetRandom(0.001,0.25);
f0gaus->SetParameters(842.556, 0.0151476, -0.000700332);
dy = f0gaus->GetRandom(0.001,0.25);
f0gaus->SetParameters(1169.83, 0.0122328, 0.000506295);
dz = f0gaus->GetRandom(0.001,0.25);
     }else if(dcent9==6){
f0gaus->SetParameters(1102.43, 0.0172002, -0.00111562);
dx = f0gaus->GetRandom(0.001,0.25);
f0gaus->SetParameters(1081.19, 0.0173019, -0.00113609);
dy = f0gaus->GetRandom(0.001,0.25);
f0gaus->SetParameters(1382.73, 0.0140201, 0.000903965);
dz = f0gaus->GetRandom(0.001,0.25);
     }else if(dcent9==5){
f0gaus->SetParameters(796.088, 0.0207058, 0.0015559);
dx = f0gaus->GetRandom(0.001,0.25);
f0gaus->SetParameters(783.342, 0.0208394, 0.00157652);
dy = f0gaus->GetRandom(0.001,0.25);
f0gaus->SetParameters(1036.28, 0.0169582, 0.00119988);
dz = f0gaus->GetRandom(0.001,0.25);
}else if(dcent9==4){
f0gaus->SetParameters(557.038, 0.025574, 0.00229113);
dx = f0gaus->GetRandom(0.001,0.25);
f0gaus->SetParameters(552.145, 0.0257164, 0.0023139);
dy = f0gaus->GetRandom(0.001,0.25);
f0gaus->SetParameters(730.544, 0.0210279, 0.00175009);
dz = f0gaus->GetRandom(0.001,0.25);
}else if(dcent9==3){
f0gaus->SetParameters(358.873, 0.0325547, 0.00344104);
dx = f0gaus->GetRandom(0.001,0.25);
f0gaus->SetParameters(355.21, 0.0327452, -0.00346166);
dy = f0gaus->GetRandom(0.001,0.25);
f0gaus->SetParameters(469.756, 0.0268485, 0.00263286);
dz = f0gaus->GetRandom(0.001,0.25);
}else if(dcent9==2){
f0gaus->SetParameters(217.214, 0.0429084, 0.00542956);
dx = f0gaus->GetRandom(0.001,0.25);
f0gaus->SetParameters(213.714, 0.0432479, 0.00547868);
dy = f0gaus->GetRandom(0.001,0.25);
f0gaus->SetParameters(283.279, 0.0354599, 0.00424226);
dz = f0gaus->GetRandom(0.001,0.25);
}else if(dcent9==1){
f0gaus->SetParameters(117.541, 0.0593891, -0.00938783);
dx = f0gaus->GetRandom(0.001,0.25);
f0gaus->SetParameters(116.556, 0.0595198, 0.00956014);
dy = f0gaus->GetRandom(0.001,0.25);
f0gaus->SetParameters(157.167, 0.049147, -0.00720776);
dz = f0gaus->GetRandom(0.001,0.25);
}else if(dcent9==0){
//f0land->SetParameters(53.4994, 0.0857768, 0.0181559);
f0land->SetParameters(0.0743963, 0.00844436, 333.871);
dx = f0land->GetRandom(0.001,0.25);
//f0land->SetParameters(54.0003, 0.0863724, 0.0178055);
f0land->SetParameters(0.074785, 0.00844436, 331.799);
dy = f0land->GetRandom(0.001,0.25);
//f0land->SetParameters(73.9683, 0.0712166, -0.0133135);
f0land->SetParameters(0.0622631, 0.00606179, 451.85);
dz = f0land->GetRandom(0.001,0.25);
}else{
dx = 0;
dy = 0;
dz = 0;
}
//cout<<"dcent9:"<<dcent9<<endl;

    kfVertex.SetCovarianceMatrix( dx*dx, 0, dy*dy, 0, 0, dz*dz );

    UShort_t noTracks = Vtx->noTracks();
    kfVertex.SetNContributors(noTracks);
    kfVertex.SetChi2(Vtx->chiSquared());

//    cout<<"noTracks:"<<noTracks<<" "<<Vtx->chiSquared()<<endl;

    primaryVertex = KFVertex(kfVertex);
  }  
   if(!IsGoodPV(primaryVertex)) return 0;

  Int_t nGlobalTracks = muDst->numberOfGlobalTracks();
  
  fParticles.resize(nGlobalTracks*7);
  fNHftHits.resize(nGlobalTracks*7);
  fParticlesPdg.resize(nGlobalTracks*7);
  int nPartSaved = 0;
  int nUsedTracks = 0;
  
  for (Int_t iTrack = 0; iTrack < nGlobalTracks; iTrack++) 
  {
    StMuTrack *gTrack = muDst->globalTracks(iTrack);
    if (! gTrack)            continue;
    if (! gTrack->charge())  continue;
    if (  gTrack->flag() < 100 ||  gTrack->flag()%100 == 11) continue; // bad fit or short track pointing to EEMC
    if (  gTrack->flag() > 1000) continue;  // pile up track in TPC
    if (  gTrack->nHitsFit() < 15) continue;

//cout<<"gTrack->probPidTraits().dEdxErrorFit():"<<gTrack->probPidTraits().dEdxErrorFit()<<endl;
    if (  gTrack->probPidTraits().dEdxErrorFit() < 0.04 || gTrack->probPidTraits().dEdxErrorFit() > 0.12 ) continue;
//    if (  gTrack->probPidTraits().dEdxErrorFit() < 0.08 || gTrack->probPidTraits().dEdxErrorFit() > 0.24 ) continue;//the latter is for run 11 27 GeV only due to bug
    
    int nHftHitsInTrack = gTrack->nHitsFit(kIstId) + gTrack->nHitsFit(kSsdId) + gTrack->nHitsFit(kPxlId);
    if(fCollectTrackHistograms) fTrackHistograms[0]->Fill(nHftHitsInTrack);
    if(fUseHFTTracksOnly && nHftHitsInTrack < 3) continue;
    
    const int index = gTrack->id();
  //cout<<"gTrack->idTruth():"<<gTrack->idTruth()<<  " "<< processSignal<<endl   ;
    mcIndices[index] = gTrack->idTruth()-1;
    if(mcIndices[index] >= int(mcTracks.size()))
      mcIndices[index] = -1;
    if(mcIndices[index] > -1)
    {
      mcTracks[mcIndices[index]].SetReconstructed();
      if(!processSignal) continue;
    }
//    else if(processSignal) continue;
    
    Int_t dcaGeometryIndex = gTrack->index2Cov();
    if (dcaGeometryIndex < 0) continue;
    StDcaGeometry *dcaG = StMuDst::instance()->covGlobTracks(dcaGeometryIndex);
    if (! dcaG) continue;
      
    THelixTrack t = dcaG->thelix();
    StThreeVectorD V(muDst->primaryVertex()->position());
    Double_t dca3D = t.Dca(V.xyz());
//    if(dca3D > 50.) continue;
      
    Int_t q = 1; if (gTrack->charge() < 0) q = -1;
    KFPTrack track;
    if( !GetTrack(*dcaG, track, q, index) ) continue;
    
    if(fCollectTrackHistograms)
    {
      fTrackHistograms2D[0]->Fill(track.GetP(), gTrack->dEdx()*1.e6);
      if(q>0) fTrackHistograms2D[1]->Fill(track.GetP(), gTrack->dEdx()*1.e6);
      else    fTrackHistograms2D[2]->Fill(track.GetP(), gTrack->dEdx()*1.e6);  
    }
    
    const StMuBTofPidTraits &btofPid = gTrack->btofPidTraits();
    double timeTof = btofPid.timeOfFlight();
    double lengthTof = btofPid.pathLength();
    if(lengthTof < 0.)
    {
      const StThreeVectorF & tofPoint  = btofPid.position();
      const StThreeVectorF & dcaPoint  = gTrack->dca(bestPV);
      StPhysicalHelixD innerHelix = gTrack->helix();
      double dlDCA = fabs( innerHelix.pathLength( StThreeVector<double>(dcaPoint.x(), dcaPoint.y(), dcaPoint.z()) ) );
      StPhysicalHelixD outerHelix = gTrack->outerHelix();
      double dlTOF = fabs( outerHelix.pathLength( StThreeVector<double>(tofPoint.x(), tofPoint.y(), tofPoint.z()) ) );
      
      double l = gTrack->length();
      lengthTof = l + dlDCA + dlTOF;
    }
    double m2tof = -1.e6;
    bool isTofm2 = false;
    if(timeTof > 0. && lengthTof > 0.)
    {
      m2tof = track.GetP()*track.GetP()*(1./((lengthTof/timeTof/29.9792458)*(lengthTof/timeTof/29.9792458))-1.);
      isTofm2 = true;
      
      if(fCollectTrackHistograms)
      {
        fTrackHistograms2D[3]->Fill(track.GetP(), gTrack->dEdx()*1.e6);
        fTrackHistograms2D[4]->Fill(track.GetP(), m2tof);
      }
    }

     double dEdXPull[7] = { fabs(gTrack->dEdxPull(0.139570, fdEdXMode, 1)),   //0 - pi
                            fabs(gTrack->dEdxPull(0.493677, fdEdXMode, 1)),   //1 - K
                            fabs(gTrack->dEdxPull(0.938272, fdEdXMode, 1)),   //2 - p
                            fabs(gTrack->dEdxPull(1.876124, fdEdXMode, 1)),   //3 - d
                            fabs(gTrack->dEdxPull(2.809432, fdEdXMode, 1)),   //4 - t
                            fabs(gTrack->dEdxPull(2.809413, fdEdXMode, 2)),   //5 - He3
                            fabs(gTrack->dEdxPull(3.728400, fdEdXMode, 2))};  //6 - He4

//     cout<<"isTofm2:"<<isTofm2<<endl;

//the following is for bug in run 11
/*                            
    double dEdXPull[7] = { fabs(2*gTrack->nSigmaPion()),   //0 - pi
                           fabs(2*gTrack->nSigmaKaon()),   //1 - K
                           fabs(2*gTrack->nSigmaProton()),   //2 - p
                           fabs(gTrack->dEdxPull(1.876124, fdEdXMode, 1)),   //3 - d
                           fabs(gTrack->dEdxPull(2.809432, fdEdXMode, 1)),   //4 - t
                           fabs(gTrack->dEdxPull(2.809413, fdEdXMode, 2)),   //5 - He3
                           fabs(gTrack->dEdxPull(3.728400, fdEdXMode, 2))};  //6 - He4
*/   
    vector<int> totalPDG = GetPID(m2tof, track.GetP(), q, gTrack->dEdx()*1.e6, dEdXPull, isTofm2, index);
        
    AddTrackToParticleList(track, nHftHitsInTrack, index, totalPDG, primaryVertex, primaryTrackList, fNHftHits, fParticlesPdg, fParticles, nPartSaved);         
    nUsedTracks++;
  }

  fParticles.resize(nPartSaved);
  fParticlesPdg.resize(nPartSaved);
  fNHftHits.resize(nPartSaved);

  if( fCleanLowPVTrackEvents && ( 10*primaryTrackList.size() < (nUsedTracks - primaryTrackList.size()) ) ) return 0;

  const Double_t field = muDst->event()->magneticField();
  SetField(field);

  CleanPV();
  InitParticles();

  //read PV
  AddPV(primaryVertex, primaryTrackList);
  if(fCollectTrackHistograms)
  {
    fTrackHistograms[1]->Fill(sqrt(dx*dx + dy*dy));
    fTrackHistograms2D[5]->Fill( nPartSaved, sqrt(dx*dx + dy*dy) );
    fTrackHistograms2D[6]->Fill( primaryTrackList.size(), sqrt(dx*dx + dy*dy) );
  }  
  //reconstruct short-lived particles
  ReconstructParticles();
  
  return 1;
}

//b
bool StKFParticleInterface::ProcessEvent(StMuDst* muDst, std::vector<int>& triggeredTracks, bool processSignal)
{

  triggeredTracks.resize(0);
  
//  mcTracks.resize(muDst->numberOfMcTracks());
//  for (unsigned int iMCTrack=0; iMCTrack<muDst->numberOfMcTracks(); iMCTrack++) 
//  {
//    StMuMcTrack *mcTrack = muDst->MCtrack(iMCTrack);
//    if (! mcTrack) continue;    
//    KFMCTrack &mcTrackKF = mcTracks[iMCTrack];
//    mcTrack->FillKFMCTrack(mcTrackKF);
//    mcTrackKF.SetNMCPixelPoints(mcTrack->No_ist_hit() + mcTrack->No_ssd_hit() + mcTrack->No_pix_hit());
//  }
  
  //read PV
  KFVertex primaryVertex;
  vector<int> primaryTrackList;

  float bestRank=-1000000;
  int bestPV=-1;
  double dx = 0., dy = 0., dz = 0.;
  for(unsigned int iPV=0; iPV<muDst->numberOfPrimaryVertices(); iPV++) 
  {
    StMuPrimaryVertex *Vtx = muDst->primaryVertex(iPV);
    if(!Vtx) continue;
    if (bestRank < Vtx->ranking()) {
      bestRank = Vtx->ranking();
      bestPV = iPV;
    }
    else continue;
    
    //convert StMuPrimaryVertex to KFVertex
    KFPVertex kfVertex;
    kfVertex.SetXYZ(Vtx->position().x(), Vtx->position().y(), Vtx->position().z());
    dx = Vtx->posError().x();
    dy = Vtx->posError().y();
    dz = Vtx->posError().z();
    kfVertex.SetCovarianceMatrix( dx*dx, 0, dy*dy, 0, 0, dz*dz );
    UShort_t noTracks = Vtx->noTracks();
    kfVertex.SetNContributors(noTracks);
    kfVertex.SetChi2(Vtx->chiSquared());
    primaryVertex = KFVertex(kfVertex);
  }  
   if(!IsGoodPV(primaryVertex)) return 0;

  Int_t nGlobalTracks = muDst->numberOfGlobalTracks();
 
  //DEBUG
  //cout<<"No. of tracks:"<<nGlobalTracks<<endl;
 
  fParticles.resize(nGlobalTracks*7);
  fNHftHits.resize(nGlobalTracks*7);
  fParticlesPdg.resize(nGlobalTracks*7);
  int nPartSaved = 0;
  int nUsedTracks = 0;
  
  for (Int_t iTrack = 0; iTrack < nGlobalTracks; iTrack++) 
  {
    StMuTrack *gTrack = muDst->globalTracks(iTrack);
    if (! gTrack)            continue;
    if (! gTrack->charge())  continue;

    if (  gTrack->flag() < 100 ||  gTrack->flag()%100 == 11) continue; // bad fit or short track pointing to EEMC
    if (  gTrack->flag() > 1000) continue;  // pile up track in TPC

//        cout<<"nGlobalTracks:"<<nGlobalTracks<<" "<<iTrack<<" "<< gTrack->nHitsFit()<<endl;

    if (  gTrack->nHitsFit() < 15) continue;

 //           cout<<"nGlobalTracks:"<<nGlobalTracks<<" "<<iTrack<<" "<< gTrack->nHitsFit()<<endl; 

//    if (  gTrack->probPidTraits().dEdxErrorFit() < 0.04 || gTrack->probPidTraits().dEdxErrorFit() > 0.12 ) continue;
//DEBUG RUN11 #1of2    
    if (  gTrack->probPidTraits().dEdxErrorFit() < 0.08 || gTrack->probPidTraits().dEdxErrorFit() > 0.24 ) continue;
//    cout<<"gTrack->probPidTraits().dEdxErrorFit():"<<gTrack->probPidTraits().dEdxErrorFit()<<endl;
  //       cout<<"nGlobalTracks:"<<nGlobalTracks<<" "<<iTrack<<" "<< gTrack->probPidTraits().dEdxErrorFit() <<endl;

    int nHftHitsInTrack = gTrack->nHitsFit(kIstId) + gTrack->nHitsFit(kSsdId) + gTrack->nHitsFit(kPxlId);
    if(fCollectTrackHistograms) fTrackHistograms[0]->Fill(nHftHitsInTrack);
    if(fUseHFTTracksOnly && nHftHitsInTrack < 3) continue;
    const int index = gTrack->id();
//    mcIndices[index] = gTrack->idTruth()-1;
//    if(mcIndices[index] >= int(mcTracks.size()))
//      mcIndices[index] = -1;
//    if(mcIndices[index] > -1)
//    {
//      mcTracks[mcIndices[index]].SetReconstructed();
//      if(!processSignal) continue;
//    }
//    else if(processSignal) continue;
    
    Int_t dcaGeometryIndex = gTrack->index2Cov();
    if (dcaGeometryIndex < 0) continue;
    StDcaGeometry *dcaG = StMuDst::instance()->covGlobTracks(dcaGeometryIndex);
    if (! dcaG) {cout <<"CANNOT FIND THE COVARIANCE MATRIX"<<endl;}
    if (! dcaG) continue;

    THelixTrack t = dcaG->thelix();
    StThreeVectorD V(muDst->primaryVertex()->position());
    Double_t dca3D = t.Dca(V.xyz());
//    if(dca3D > 50.) continue;
      
    Int_t q = 1; if (gTrack->charge() < 0) q = -1;
    KFPTrack track;
    if( !GetTrack(*dcaG, track, q, index) ) continue;
    
    if(fCollectTrackHistograms)
    {
      fTrackHistograms2D[0]->Fill(track.GetP(), gTrack->dEdx()*1.e6);
      if(q>0) fTrackHistograms2D[1]->Fill(track.GetP(), gTrack->dEdx()*1.e6);
      else    fTrackHistograms2D[2]->Fill(track.GetP(), gTrack->dEdx()*1.e6);  
    }
    
    const StMuBTofPidTraits &btofPid = gTrack->btofPidTraits();
    double timeTof = btofPid.timeOfFlight();
    double lengthTof = btofPid.pathLength();
    if(lengthTof < 0.)
    {
      const StThreeVectorF & tofPoint  = btofPid.position();
      const StThreeVectorF & dcaPoint  = gTrack->dca(bestPV);
      StPhysicalHelixD innerHelix = gTrack->helix();
      double dlDCA = fabs( innerHelix.pathLength( StThreeVector<double>(dcaPoint.x(), dcaPoint.y(), dcaPoint.z()) ) );
      StPhysicalHelixD outerHelix = gTrack->outerHelix();
      double dlTOF = fabs( outerHelix.pathLength( StThreeVector<double>(tofPoint.x(), tofPoint.y(), tofPoint.z()) ) );
      
      double l = gTrack->length();
      lengthTof = l + dlDCA + dlTOF;
    }
    double m2tof = -1.e6;
    bool isTofm2 = false;
    if(timeTof > 0. && lengthTof > 0.)
    {
      m2tof = track.GetP()*track.GetP()*(1./((lengthTof/timeTof/29.9792458)*(lengthTof/timeTof/29.9792458))-1.);
      isTofm2 = true;
      
      if(fCollectTrackHistograms)
      {
        fTrackHistograms2D[3]->Fill(track.GetP(), gTrack->dEdx()*1.e6);
        fTrackHistograms2D[4]->Fill(track.GetP(), m2tof);
      }
    }
//b testing
/*
     double dEdXPull[7] = { fabs(gTrack->dEdxPull(0.139570, fdEdXMode, 1)),   //0 - pi
                            fabs(gTrack->dEdxPull(0.493677, fdEdXMode, 1)),   //1 - K
                            fabs(gTrack->dEdxPull(0.938272, fdEdXMode, 1)),   //2 - p
                            fabs(gTrack->dEdxPull(1.876124, fdEdXMode, 1)),   //3 - d
                            fabs(gTrack->dEdxPull(2.809432, fdEdXMode, 1)),   //4 - t
                            fabs(gTrack->dEdxPull(2.809413, fdEdXMode, 2)),   //5 - He3
                            fabs(gTrack->dEdxPull(3.728400, fdEdXMode, 2))};  //6 - He4

*/
//
////DEBUG RUN11 #2of2

    double dEdXPull[7] = { fabs(2*gTrack->nSigmaPion()),   //0 - pi
                           fabs(2*gTrack->nSigmaKaon()),   //1 - K
                           fabs(2*gTrack->nSigmaProton()),   //2 - p
                           fabs(gTrack->dEdxPull(1.876124, fdEdXMode, 1)),   //3 - d
                           fabs(gTrack->dEdxPull(2.809432, fdEdXMode, 1)),   //4 - t
                           fabs(gTrack->dEdxPull(2.809413, fdEdXMode, 2)),   //5 - He3
                           fabs(gTrack->dEdxPull(3.728400, fdEdXMode, 2))};  //6 - He4



//cout<<m2tof<<" "<<track.GetP()<<" "<<q<<" "<<gTrack->dEdx()*1.e6<<" "<<dEdXPull[0]  <<" "<<isTofm2<<" "<<index<<endl;
//cout<<gTrack->dEdx()*1.e6<<" "<<dEdXPull[0]<<" " <<gTrack->nSigmaPion()<<endl;
    
    vector<int> totalPDG = GetPID(m2tof, track.GetP(), q, gTrack->dEdx()*1.e6, dEdXPull, isTofm2, index);
    
    int nPartSaved0 = nPartSaved;    
    AddTrackToParticleList(track, nHftHitsInTrack, index, totalPDG, primaryVertex, primaryTrackList, fNHftHits, fParticlesPdg, fParticles, nPartSaved);         
 

if(nPartSaved > nPartSaved0)
      triggeredTracks.push_back(iTrack);

    nUsedTracks++;
  }

  fParticles.resize(nPartSaved);
  fParticlesPdg.resize(nPartSaved);
  fNHftHits.resize(nPartSaved);

  if( fCleanLowPVTrackEvents && ( 10*primaryTrackList.size() < (nUsedTracks - primaryTrackList.size()) ) ) return 0;

  const Double_t field = muDst->event()->magneticField();
  SetField(field);

  CleanPV();
  InitParticles();

  //read PV
  AddPV(primaryVertex, primaryTrackList);
  if(fCollectTrackHistograms)
  {
    fTrackHistograms[1]->Fill(sqrt(dx*dx + dy*dy));
    fTrackHistograms2D[5]->Fill( nPartSaved, sqrt(dx*dx + dy*dy) );
    fTrackHistograms2D[6]->Fill( primaryTrackList.size(), sqrt(dx*dx + dy*dy) );
  }  
  //reconstruct short-lived particles
  ReconstructParticles();
  
  return 1;

}
