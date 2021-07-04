//*-- Author : Yuri Fisyak 02/02/2016
#include "StKFParticleAnalysisMaker.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TSystem.h"
//--- KF particle classes ---
#include "KFVertex.h"
#include "KFParticle.h"
#include "KFParticleSIMD.h"
#include "KFPTrack.h"
#include "KFParticleTopoReconstructor.h"
#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"
//--- Pico classes ---
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
//--- Mu classes ---
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"
//--- TMVA classes ---
#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/Factory.h"
//--- StRefMult class ---
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"
ClassImp(StKFParticleAnalysisMaker);

//________________________________________________________________________________
StKFParticleAnalysisMaker::StKFParticleAnalysisMaker(const char *name) : StMaker(name), fNTrackTMVACuts(0), fIsPicoAnalysis(true), fsnn(-999), fdEdXMode(1), 
  fStoreTmvaNTuples(false), fProcessSignal(false), fCollectTrackHistograms(false), fCollectPIDHistograms(false),fTMVAselection(false), fStoremctree(false),
  fFlowAnalysis(false), fFlowChain(NULL), fFlowRunId(-1), fFlowEventId(-1), fCentrality(-1), fFlowFiles(), fFlowMap(), 
  fRunCentralityAnalysis(0), fRefmultCorrUtil(0), fCentralityFile(""), fAnalyseDsPhiPi(false)
{
  memset(mBeg,0,mEnd-mBeg+1);

/*  
  fNTuplePDG[0] = 421;
  fNTuplePDG[1] = 411;
  fNTuplePDG[2] = 431;
  fNTuplePDG[3] = 4122;
  fNTuplePDG[4] = 426;
  fNTuplePDG[5] = 429;
  fNTuplePDG[6] = 521;
  fNTuplePDG[7] = 511;
  fNTuplePDG[8] = 3122;
*/
  fNTuplePDG[0] = 3122;
  //fNTuplePDG[1] = 3334;
  //TODO
  fNTuplePDG[1] = 3334;
  fNTuplePDG[2] = 3312;
  
/*
  fNtupleNames[0] = "D0"; 
  fNtupleNames[1] = "DPlus"; 
  fNtupleNames[2] = "Ds"; 
  fNtupleNames[3] = "Lc";
  fNtupleNames[4] = "D0KK";
  fNtupleNames[5] = "D04";
  fNtupleNames[6] = "BPlus";
  fNtupleNames[7] = "B0";
  fNtupleNames[8] = "Ld";
*/
  fNtupleNames[0] = "Ld";
  fNtupleNames[1] = "Om";
  fNtupleNames[2] = "Xi";
  
  vector<TString> trackCutNames;
  trackCutNames.push_back("pt_");
  trackCutNames.push_back("chi2Primary_");
  trackCutNames.push_back("dEdXPi_");
  trackCutNames.push_back("dEdXK_");
  trackCutNames.push_back("dEdXP_");
  trackCutNames.push_back("ToFPi_");
  trackCutNames.push_back("ToFK_");
  trackCutNames.push_back("ToFP_");
  fNTrackTMVACuts = trackCutNames.size();

/*  
  fDaughterNames[0].push_back("K");     fDaughterNames[0].push_back("Pi");                                                                              //D0 -> Kpi
  fDaughterNames[1].push_back("K");     fDaughterNames[1].push_back("Pi1");    fDaughterNames[1].push_back("Pi2");                                      //D+ -> Kpipi
  fDaughterNames[2].push_back("KPlus"); fDaughterNames[2].push_back("KMinus"); fDaughterNames[2].push_back("Pi");                                       //Ds -> KKpi
  fDaughterNames[3].push_back("K");     fDaughterNames[3].push_back("Pi");     fDaughterNames[3].push_back("P");                                        //Lc -> pKpi
  fDaughterNames[4].push_back("KPlus"); fDaughterNames[4].push_back("KMinus");                                                                          //D0 -> KK
  fDaughterNames[5].push_back("K");     fDaughterNames[5].push_back("Pi1");    fDaughterNames[5].push_back("Pi2");  fDaughterNames[5].push_back("Pi3"); //D0 -> Kpipipi
  fDaughterNames[6].push_back("PiD");   fDaughterNames[6].push_back("KD");     fDaughterNames[6].push_back("Pi");                                       //B+ -> D0_bpi
  fDaughterNames[7].push_back("Pi1D");  fDaughterNames[7].push_back("KD");     fDaughterNames[7].push_back("Pi2D"); fDaughterNames[7].push_back("Pi");  //B0 -> D-pi+
//  fDaughterNames[8].push_back("Proton");fDaughterNames[8].push_back("Pi");
*/
  fDaughterNames[0].push_back("Proton");fDaughterNames[0].push_back("Pi");
  fDaughterNames[1].push_back("Kaon");fDaughterNames[1].push_back("Proton");fDaughterNames[1].push_back("Pion");
  fDaughterNames[2].push_back("Pion_bach");fDaughterNames[2].push_back("Proton");fDaughterNames[2].push_back("Pion");

  for(int iDecay=0; iDecay<fNNTuples; iDecay++)
  {
    for(unsigned int iDaughter=0; iDaughter<fDaughterNames[iDecay].size(); iDaughter++)
    {
      for(int iTrackTMVACut=0; iTrackTMVACut<fNTrackTMVACuts; iTrackTMVACut++)
      {
        if(iDaughter==0 && iTrackTMVACut==0)
          fNtupleCutNames[iDecay] = trackCutNames[iTrackTMVACut];  
        else
          fNtupleCutNames[iDecay] += trackCutNames[iTrackTMVACut];
        fNtupleCutNames[iDecay] += fDaughterNames[iDecay][iDaughter];
        fNtupleCutNames[iDecay] += ":";
      }
    }
    if(iDecay==0) 
      fNtupleCutNames[iDecay] += "Chi2NDF:LdL:L:Chi2Topo:refMult:refMultCor:reweight:cent9:mass:px:py:pz:vx:vy:vz:evx:evy:evz";	
    else if(iDecay==1 || iDecay==2)
      //fNtupleCutNames[iDecay] += "Chi2NDF:LdL:L:Chi2Topo:refMult:mass:pt";
      fNtupleCutNames[iDecay] += "Chi2NDF:LdL:L:Chi2Topo:refMult:refMultCor:reweight:cent9:mass:px:py:pz:vx:vy:vz:evx:evy:evz:Chi2NDF_ld:LdL_ld:L_ld:Chi2Topo_ld:pid";
    else if(iDecay>2 && iDecay<6 )
      fNtupleCutNames[iDecay] += "Chi2NDF:LdL:Chi2Topo:refMult";
    else if(iDecay>=6 && iDecay<8)
    {
      fNtupleCutNames[iDecay] += "Chi2NDF_D:LdL_D:Chi2Topo_D:Chi2NDF:LdL:Chi2Topo:refMult";
    } 
    
    SetTMVABins(iDecay);
  }
}
//________________________________________________________________________________
StKFParticleAnalysisMaker::~StKFParticleAnalysisMaker() 
{
  SafeDelete(fStKFParticleInterface);
  SafeDelete(fStKFParticlePerformanceInterface);
}

//_____________________________________________________________________________
Int_t StKFParticleAnalysisMaker::Init()
{
  TFile *f = GetTFile();
  if(f) 
  {
    f->cd();
    BookVertexPlots();
    if(fCollectTrackHistograms)
      fStKFParticleInterface->CollectTrackHistograms();
    if(fCollectPIDHistograms)
      fStKFParticleInterface->CollectPIDHistograms();
  }
  
  if(fTMVAselection || fStoreTmvaNTuples)
  {
    for(int iReader=0; iReader<fNNTuples; iReader++)
    {
      TString cutName;
      int firstSymbolOfCutName = 0;
      
      int nCuts = 0;
      while(fNtupleCutNames[iReader].Tokenize(cutName,firstSymbolOfCutName,":"))
        nCuts++;
      fTMVAParticleParameters[iReader].resize(nCuts);
    }
  }
  
  if(fTMVAselection)
  {
    for(int iReader=0; iReader<fNNTuples; iReader++)
    {
      const int nCentralityBins = fTMVACentralityBins[iReader].size() - 1;
      const int nPtBins = fTMVAPtBins[iReader].size() - 1;
      
      for(int iCentralityBin=0; iCentralityBin<nCentralityBins; iCentralityBin++)
      {
        for(int iPtBin=0; iPtBin<nPtBins; iPtBin++)
        {
          fTMVAReader[iReader][iCentralityBin][iPtBin] = new TMVA::Reader("Silent");

          TString cutName;
          int firstSymbolOfCutName = 0;      
          unsigned int iCut = 0;
          while(fNtupleCutNames[iReader].Tokenize(cutName,firstSymbolOfCutName,":"))
          {
            fTMVAReader[iReader][iCentralityBin][iPtBin] -> AddVariable( cutName.Data(), &fTMVAParticleParameters[iReader][iCut] );
            iCut++;
            if(iCut == (fTMVAParticleParameters[iReader].size()-1)) break;
          }
          
          fTMVAReader[iReader][iCentralityBin][iPtBin] -> BookMVA("BDT", fTMVACutFile[iReader][iCentralityBin][iPtBin].Data());
        }
      }
    }
  }
      
  //Create file with NTuples for cut optimization
  if(fStoreTmvaNTuples)
  {  
    TFile* curFile = gFile;
    TDirectory* curDirectory = gDirectory;
    for(int iNtuple=0; iNtuple<fNNTuples; iNtuple++)
    {
      TString SignalPrefix = "_Signal";
      if(!fProcessSignal) SignalPrefix = "_BG";
      TString currentNTupleFileName = fNtupleNames[iNtuple]+SignalPrefix+TString(".root");
      fNTupleFile[iNtuple] = new TFile(currentNTupleFileName.Data(),"RECREATE");
      fCutsNTuple[iNtuple] = new TNtuple(fNtupleNames[iNtuple].Data(), fNtupleNames[iNtuple].Data(), fNtupleCutNames[iNtuple].Data());
    }
    gFile = curFile;
    gDirectory = curDirectory;
  }
  
  //fRefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr_P16id();
  //fRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
  //fRefmultCorrUtil->readScaleForWeight("/gpfs01/star/pwg/pfederic/qVectors/StRoot/StRefMultCorr/macros/weight_grefmult_VpdnoVtx_Vpd5_Run16.txt"); //for new StRefMultCorr, Run16, SL16j
 fRefmultCorrUtil = CentralityMaker::instance()->getRefMultCorr() ;
 
  //Initialise the chain with files containing centrality and reaction plane
  if(fFlowAnalysis)
  {
    std::cout << "StKFParticleAnalysisMaker: run flow analysis. Flow file list:"<<std::endl;
    
    //fFlowChain = new TChain("mTree");
    fFlowChain = new TChain("psi_tree");
    for(unsigned int iFlowFile=0; iFlowFile<fFlowFiles.size(); iFlowFile++)
    {
      std::cout << "      " << fFlowFiles[iFlowFile] << std::endl;
      fFlowChain->Add(fFlowFiles[iFlowFile].Data());
    }
    
    fFlowChain->SetBranchStatus("*",0);
    //fFlowChain->SetBranchAddress("runid",   &fFlowRunId);   fFlowChain->SetBranchStatus("runid", 1);
    fFlowChain->SetBranchAddress("runnumber",   &fFlowRunId);   fFlowChain->SetBranchStatus("runnumber", 1);
    //fFlowChain->SetBranchAddress("eventid", &fFlowEventId); fFlowChain->SetBranchStatus("eventid", 1);
    fFlowChain->SetBranchAddress("eventid", &fFlowEventId); fFlowChain->SetBranchStatus("eventid", 1);
    //fFlowChain->SetBranchAddress("cent", &fCentrality);  fFlowChain->SetBranchStatus("cent", 1);
    fFlowChain->SetBranchAddress("centnumber", &fCentrality);  fFlowChain->SetBranchStatus("centnumber", 1);

    fFlowChain->SetBranchAddress("psi_1_EPD_0", &psi_1_EPD_0);
    fFlowChain->SetBranchStatus("psi_1_EPD_0", 1);
    fFlowChain->SetBranchAddress("psi_1_EPD_1", &psi_1_EPD_1);
    fFlowChain->SetBranchStatus("psi_1_EPD_1", 1);
    fFlowChain->SetBranchAddress("psi_1_EPD_2", &psi_1_EPD_2);
    fFlowChain->SetBranchStatus("psi_1_EPD_2", 1);
    fFlowChain->SetBranchAddress("psi_1_EPD_3", &psi_1_EPD_3);
    fFlowChain->SetBranchStatus("psi_1_EPD_3", 1);
    fFlowChain->SetBranchAddress("psi_1_EPD_4", &psi_1_EPD_4);
    fFlowChain->SetBranchStatus("psi_1_EPD_4", 1);
    fFlowChain->SetBranchAddress("psi_1_EPD_5", &psi_1_EPD_5);
    fFlowChain->SetBranchStatus("psi_1_EPD_5", 1);
    fFlowChain->SetBranchAddress("psi_1_EPD_6", &psi_1_EPD_6);
    fFlowChain->SetBranchStatus("psi_1_EPD_6", 1);
    fFlowChain->SetBranchAddress("psi_1_EPD", &psi_1_EPD);
    fFlowChain->SetBranchStatus("psi_1_EPD", 1);
    fFlowChain->SetBranchAddress("gweight", &gweight);
    fFlowChain->SetBranchStatus("gweight", 1);
    
    std::cout << "StKFParticleAnalysisMaker: number of entries in the flow chain" << fFlowChain->GetEntries() << std::endl;
    for(int iEntry=0; iEntry<fFlowChain->GetEntries(); iEntry++)
    {
      fFlowChain->GetEvent(iEntry);
      fFlowMap[GetUniqueEventId(fFlowRunId, fFlowEventId)] = iEntry;
    }
  }


   file_out = new TFile("ana_tree.root","RECREATE");
   lambda_tree = new TTree("lambda_tree","ana_tree");

//  lambda_tree->Branch("bVz",&bVz,"bVz/F");
//  lambda_tree->Branch("bVr",&bVr,"bVr/F");
//  lambda_tree->Branch("bVzerr",&bVzerr,"bVzerr/F");
//  lambda_tree->Branch("bVrerr",&bVrerr,"bVrerr/F");
 
   lambda_tree->Branch("brunId",&brunid,"brunid/I");
   lambda_tree->Branch("beventid",&beventid,"beventid/I");
   lambda_tree->Branch("bVz",&bVz,"bVz/F");
   lambda_tree->Branch("brefmult",&brefmult,"brefmult/I");
   lambda_tree->Branch("btofmult",&btofmult,"btofmult/I");
   lambda_tree->Branch("countrefmult",&countrefmult,"countrefmult/I");
 
   lambda_tree->Branch("ld_chi2topo",&ld_chi2topo,"ld_chi2topo/F");
   lambda_tree->Branch("ld_chi2ndf",&ld_chi2ndf,"ld_chi2ndf/F");
   lambda_tree->Branch("ld_ldl",&ld_ldl,"ld_ldl/F");
   lambda_tree->Branch("ld_l",&ld_l,"ld_l/F");
   lambda_tree->Branch("ld_dl",&ld_dl,"ld_dl/F");

   lambda_tree->Branch("dca_proton",&dca_proton,"dca_proton/F");
   lambda_tree->Branch("dca_pi",&dca_pi,"dca_pi/F");    

   lambda_tree->Branch("chi2primary_proton",&chi2primary_proton,"chi2primary_proton/F");
   lambda_tree->Branch("chi2primary_pi",&chi2primary_pi,"chi2primary_pi/F");

   lambda_tree->Branch("nhits_ld_proton",&nhits_ld_proton,"nhits_ld_proton/I");
   lambda_tree->Branch("nhits_ld_pi",&nhits_ld_pi,"nhits_ld_pi/I");

   lambda_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
   lambda_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");

   lambda_tree->Branch("bx",&bx,"bx/F");
   lambda_tree->Branch("by",&by,"by/F");
   lambda_tree->Branch("bz",&bz,"bz/F");
   lambda_tree->Branch("bpx",&bpx,"bpx/F");
   lambda_tree->Branch("bpy",&bpy,"bpy/F");
   lambda_tree->Branch("bpz",&bpz,"bpz/F");
   lambda_tree->Branch("bpl",&bpl,"bpl/F");
 
   lambda_tree->Branch("notbadrun",&notbadrun,"notbadrun/I");
 
   lambda_tree->Branch("bmcpx",&bmcpx,"bmcpx/F");
   lambda_tree->Branch("bmcpy",&bmcpy,"bmcpy/F");
   lambda_tree->Branch("bmcpz",&bmcpz,"bmcpz/F");
   lambda_tree->Branch("bmcx",&bmcx,"bmcx/F");
   lambda_tree->Branch("bmcy",&bmcy,"bmcy/F");
   lambda_tree->Branch("bmcz",&bmcz,"bmcz/F");
   lambda_tree->Branch("bmcl",&bmcl,"bmcl/F");
   lambda_tree->Branch("bmcpl",&bmcpl,"bmcpl/F");

   lambda_tree->Branch("ld_bdfvtx",&ld_bdfvtx,"ld_bdfvtx/F");
   lambda_tree->Branch("ld_bdfvtx2",&ld_bdfvtx2,"ld_bdfvtx2/F");
 
   lambda_tree->Branch("ld_bdfvtx_xy",&ld_bdfvtx_xy,"ld_bdfvtx_xy/F");
   lambda_tree->Branch("ld_bdfvtxdev_xy",&ld_bdfvtxdev_xy,"ld_bdfvtxdev_xy/F");
   lambda_tree->Branch("ld_lifetime",&ld_lifetime,"ld_lifetime/F");

   lambda_tree->Branch("bismc",&bismc,"bismc/I");
 
   lambda_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
   lambda_tree->Branch("reweight",&reweight,"reweight/D");
   lambda_tree->Branch("cent9",&cent9,"cent9/I");


  if(fFlowAnalysis){
   lambda_tree->Branch("fCentrality",&fCentrality,"fCentrality/I");
   lambda_tree->Branch("psi_1_EPD_0",&psi_1_EPD_0,"psi_1_EPD_0/D");
   lambda_tree->Branch("psi_1_EPD_1",&psi_1_EPD_1,"psi_1_EPD_1/D");
   lambda_tree->Branch("psi_1_EPD_2",&psi_1_EPD_2,"psi_1_EPD_2/D");
   lambda_tree->Branch("psi_1_EPD_3",&psi_1_EPD_3,"psi_1_EPD_3/D");
   lambda_tree->Branch("psi_1_EPD_4",&psi_1_EPD_4,"psi_1_EPD_4/D");
   lambda_tree->Branch("psi_1_EPD_5",&psi_1_EPD_5,"psi_1_EPD_5/D");
   lambda_tree->Branch("psi_1_EPD_6",&psi_1_EPD_6,"psi_1_EPD_6/D");
   lambda_tree->Branch("psi_1_EPD",&psi_1_EPD,"psi_1_EPD/D");
   lambda_tree->Branch("gweight",&gweight,"gweight/F");
   }
 
   lambda_mc_tree = new TTree("lambda_mc_tree","ana_tree");
   lambda_mc_tree->Branch("brunid",&brunid,"brunid/I");
   lambda_mc_tree->Branch("brefmult",&brefmult,"brefmult/I");
   lambda_mc_tree->Branch("bmcparticleid",&bmcparticleid,"bmcparticleid/I");
   lambda_mc_tree->Branch("bmcrawpx",&bmcrawpx,"bmcrawpx/F");
   lambda_mc_tree->Branch("bmcrawpy",&bmcrawpy,"bmcrawpy/F");
   lambda_mc_tree->Branch("bmcrawpz",&bmcrawpz,"bmcrawpz/F");
   lambda_mc_tree->Branch("bmcrawl",&bmcrawl,"bmcrawl/F");
   lambda_mc_tree->Branch("bmcrawpl",&bmcrawpl,"bmcrawpl/F");
 
   lambda_mc_tree->Branch("countrefmult",&countrefmult,"countrefmult/I");
   lambda_mc_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
   lambda_mc_tree->Branch("reweight",&reweight,"reweight/D");
   lambda_mc_tree->Branch("cent9",&cent9,"cent9/I");

   ks_mc_tree = new TTree("ks_mc_tree","ana_tree");
   ks_mc_tree->Branch("brunid",&brunid,"brunid/I");
   ks_mc_tree->Branch("brefmult",&brefmult,"brefmult/I");
   ks_mc_tree->Branch("bmcparticleid",&bmcparticleid,"bmcparticleid/I");
   ks_mc_tree->Branch("bmcrawpx",&bmcrawpx,"bmcrawpx/F");
   ks_mc_tree->Branch("bmcrawpy",&bmcrawpy,"bmcrawpy/F");
   ks_mc_tree->Branch("bmcrawpz",&bmcrawpz,"bmcrawpz/F");
   ks_mc_tree->Branch("bmcrawl",&bmcrawl,"bmcrawl/F");
   ks_mc_tree->Branch("bmcrawpl",&bmcrawpl,"bmcrawpl/F");
 
   ks_mc_tree->Branch("countrefmult",&countrefmult,"countrefmult/I");
   ks_mc_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
   ks_mc_tree->Branch("reweight",&reweight,"reweight/D");
   ks_mc_tree->Branch("cent9",&cent9,"cent9/I");
 
   htriton_mc_tree = new TTree("htriton_mc_tree","ana_tree");
   htriton_mc_tree->Branch("brunid",&brunid,"brunid/I");
   htriton_mc_tree->Branch("brefmult",&brefmult,"brefmult/I");
   htriton_mc_tree->Branch("bmcparticleid",&bmcparticleid,"bmcparticleid/I");
   htriton_mc_tree->Branch("bmcrawpx",&bmcrawpx,"bmcrawpx/F");
   htriton_mc_tree->Branch("bmcrawpy",&bmcrawpy,"bmcrawpy/F");
   htriton_mc_tree->Branch("bmcrawpz",&bmcrawpz,"bmcrawpz/F");
   htriton_mc_tree->Branch("bmcrawl",&bmcrawl,"bmcrawl/F");
   htriton_mc_tree->Branch("bmcrawpl",&bmcrawpl,"bmcrawpl/F");
 
   htriton_mc_tree->Branch("countrefmult",&countrefmult,"countrefmult/I");
   htriton_mc_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
   htriton_mc_tree->Branch("reweight",&reweight,"reweight/D");
   htriton_mc_tree->Branch("cent9",&cent9,"cent9/I");

   ks_tree = new TTree("ks_tree","ana_tree");
   ks_tree->Branch("brunId",&brunid,"brunid/I");
   ks_tree->Branch("beventid",&beventid,"beventid/I");
   ks_tree->Branch("bVz",&bVz,"bVz/F");
   ks_tree->Branch("brefmult",&brefmult,"brefmult/I");
   ks_tree->Branch("btofmult",&btofmult,"btofmult/I");
   ks_tree->Branch("ld_chi2topo",&ld_chi2topo,"ld_chi2topo/F");
   ks_tree->Branch("ld_chi2ndf",&ld_chi2ndf,"ld_chi2ndf/F");
   ks_tree->Branch("ld_ldl",&ld_ldl,"ld_ldl/F");
   ks_tree->Branch("ld_l",&ld_l,"ld_l/F");
   ks_tree->Branch("ld_dl",&ld_dl,"ld_dl/F");
   ks_tree->Branch("chi2primary_proton",&chi2primary_proton,"chi2primary_proton/F");
   ks_tree->Branch("chi2primary_pi",&chi2primary_pi,"chi2primary_pi/F");
   ks_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
   ks_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");
   ks_tree->Branch("bpx",&bpx,"bpx/F");
   ks_tree->Branch("bpy",&bpy,"bpy/F");
   ks_tree->Branch("bpz",&bpz,"bpz/F");
   ks_tree->Branch("bmcpx",&bmcpx,"bmcpx/F");
   ks_tree->Branch("bmcpy",&bmcpy,"bmcpy/F");
   ks_tree->Branch("bmcpz",&bmcpz,"bmcpz/F");
   ks_tree->Branch("bismc",&bismc,"bismc/I");
   ks_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
   ks_tree->Branch("reweight",&reweight,"reweight/D");
   ks_tree->Branch("cent9",&cent9,"cent9/I");

   if(fFlowAnalysis){
   ks_tree->Branch("fCentrality",&fCentrality,"fCentrality/I");
   ks_tree->Branch("psi_1_EPD_0",&psi_1_EPD_0,"psi_1_EPD_0/D");
   ks_tree->Branch("psi_1_EPD_1",&psi_1_EPD_1,"psi_1_EPD_1/D");
   ks_tree->Branch("psi_1_EPD_2",&psi_1_EPD_2,"psi_1_EPD_2/D");
   ks_tree->Branch("psi_1_EPD_3",&psi_1_EPD_3,"psi_1_EPD_3/D");
   ks_tree->Branch("psi_1_EPD_4",&psi_1_EPD_4,"psi_1_EPD_4/D");
   ks_tree->Branch("psi_1_EPD_5",&psi_1_EPD_5,"psi_1_EPD_5/D");
   ks_tree->Branch("psi_1_EPD_6",&psi_1_EPD_6,"psi_1_EPD_6/D");
   ks_tree->Branch("psi_1_EPD",&psi_1_EPD,"psi_1_EPD/D");
   ks_tree->Branch("gweight",&gweight,"gweight/F");
   }


   cascade_tree = new TTree("cascade_tree","ana_tree");
   cascade_tree->Branch("brunid",&brunid,"brunid/I");
   cascade_tree->Branch("beventid",&beventid,"beventid/I");

   cascade_tree->Branch("bVz",&bVz,"bVz/F");
//  cascade_tree->Branch("bVr",&bVr,"bVr/F");
//  cascade_tree->Branch("bVzerr",&bVzerr,"bVzerr/F");
//  cascade_tree->Branch("bVrerr",&bVrerr,"bVrerr/F");

   cascade_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
   cascade_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");
//  cascade_tree->Branch("bx",&bx,"bx/F");
//  cascade_tree->Branch("by",&by,"by/F");
//  cascade_tree->Branch("bz",&bz,"bz/F");
   cascade_tree->Branch("bpx",&bpx,"bpx/F");
   cascade_tree->Branch("bpy",&bpy,"bpy/F");
   cascade_tree->Branch("bpz",&bpz,"bpz/F");
//  cascade_tree->Branch("bdl",&bdl,"bdl/F");

   cascade_tree->Branch("brefmult",&brefmult,"brefmult/I");
   cascade_tree->Branch("btofmult",&btofmult,"btofmult/I");

   cascade_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
   cascade_tree->Branch("reweight",&reweight,"reweight/D");
   cascade_tree->Branch("cent9",&cent9,"cent9/I");

//  cascade_tree->Branch("bbachid",&bbachid,"bbachid/I");
//  cascade_tree->Branch("bbachpx",&bbachpx,"bbachpx/F");
//  cascade_tree->Branch("bbachpy",&bbachpy,"bbachpy/F");
//  cascade_tree->Branch("bbachpz",&bbachpz,"bbachpz/F");
//  cascade_tree->Branch("bbachmass",&bbachmass,"bbachmass/F");

//  cascade_tree->Branch("bpionid",&bpionid,"bpionid/I");
//  cascade_tree->Branch("bpionpx",&bpionpx,"bpionpx/F");
//  cascade_tree->Branch("bpionpy",&bpionpy,"bpionpy/F");
//  cascade_tree->Branch("bpionpz",&bpionpz,"bpionpz/F");
//  cascade_tree->Branch("bpionmass",&bpionmass,"bpionmass/F");

//  cascade_tree->Branch("bprotonid",&bprotonid,"bprotonid/I");
//  cascade_tree->Branch("bprotonpx",&bprotonpx,"bprotonpx/F");
//  cascade_tree->Branch("bprotonpy",&bprotonpy,"bprotonpy/F");
//  cascade_tree->Branch("bprotonpz",&bprotonpz,"bprotonpz/F");
//  cascade_tree->Branch("bprotonmass",&bprotonmass,"bprotonmass/F");

   cascade_tree->Branch("xi_chi2topo",&xi_chi2topo,"xi_chi2topo/F");
   cascade_tree->Branch("xi_chi2ndf",&xi_chi2ndf,"xi_chi2ndf/F");
   cascade_tree->Branch("xi_ldl",&xi_ldl,"xi_ldl/F");

   cascade_tree->Branch("xi_ld_chi2topo",&xi_ld_chi2topo,"xi_ld_chi2topo/F");
   cascade_tree->Branch("xi_ld_chi2ndf",&xi_ld_chi2ndf,"xi_ld_chi2ndf/F");
   cascade_tree->Branch("xi_ld_ldl",&xi_ld_ldl,"xi_ld_ldl/F");
   cascade_tree->Branch("xi_ld_l",&xi_ld_l,"xi_ld_l/F");

   cascade_tree->Branch("xi_l",&xi_l,"xi_l/F");
   cascade_tree->Branch("xi_dl",&xi_dl,"xi_dl/F");

   cascade_tree->Branch("chi2primary_xi_proton",&chi2primary_xi_proton,"chi2primary_xi_proton/F");
   cascade_tree->Branch("chi2primary_xi_pi",&chi2primary_xi_pi,"chi2primary_xi_pi/F");
   cascade_tree->Branch("chi2primary_xi_bach",&chi2primary_xi_bach,"chi2primary_xi_bach/F");
   cascade_tree->Branch("chi2primary_xi_ld",&chi2primary_xi_ld,"chi2primary_xi_ld/F");

   cascade_tree->Branch("notbadrun",&notbadrun,"notbadrun/I");

   cascade_tree->Branch("bmcpx",&bmcpx,"bmcpx/F");
   cascade_tree->Branch("bmcpy",&bmcpy,"bmcpy/F");
   cascade_tree->Branch("bmcpz",&bmcpz,"bmcpz/F");
   cascade_tree->Branch("bismc",&bismc,"bismc/I");


   cascade_mc_tree = new TTree("cascade_mc_tree","ana_tree");
   cascade_mc_tree->Branch("brunid",&brunid,"brunid/I");
   cascade_mc_tree->Branch("brefmult",&brefmult,"brefmult/I");
   cascade_mc_tree->Branch("bmcparticleid",&bmcparticleid,"bmcparticleid/I");
   cascade_mc_tree->Branch("bmcrawpx",&bmcrawpx,"bmcrawpx/F");
   cascade_mc_tree->Branch("bmcrawpy",&bmcrawpy,"bmcrawpy/F");
   cascade_mc_tree->Branch("bmcrawpz",&bmcrawpz,"bmcrawpz/F");

   cascade_mc_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
   cascade_mc_tree->Branch("reweight",&reweight,"reweight/D");
   cascade_mc_tree->Branch("cent9",&cent9,"cent9/I");

   omega_tree = new TTree("omega_tree","ana_tree");
   omega_tree->Branch("brunid",&brunid,"brunid/I");
   // omega_tree->Branch("beventid",&beventid,"beventid/I");

   omega_tree->Branch("bVz",&bVz,"bVz/F");
   // omega_tree->Branch("bVr",&bVr,"bVr/F");
   // omega_tree->Branch("bVzerr",&bVzerr,"bVzerr/F");
   // omega_tree->Branch("bVrerr",&bVrerr,"bVrerr/F");
   omega_tree->Branch("brefmult",&brefmult,"brefmult/I");
   omega_tree->Branch("btofmult",&btofmult,"btofmult/I");

   omega_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
   omega_tree->Branch("reweight",&reweight,"reweight/D");
   omega_tree->Branch("cent9",&cent9,"cent9/I");

   omega_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
   omega_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");
   // omega_tree->Branch("bx",&bx,"bx/F");
   // omega_tree->Branch("by",&by,"by/F");
   // omega_tree->Branch("bz",&bz,"bz/F");
   omega_tree->Branch("bpx",&bpx,"bpx/F");
   omega_tree->Branch("bpy",&bpy,"bpy/F");
   omega_tree->Branch("bpz",&bpz,"bpz/F");
   // omega_tree->Branch("bdl",&bdl,"bdl/F");

   omega_tree->Branch("om_chi2topo",&om_chi2topo,"om_chi2topo/F");
   omega_tree->Branch("om_chi2ndf",&om_chi2ndf,"om_chi2ndf/F");
   omega_tree->Branch("om_ldl",&om_ldl,"om_ldl/F");

   omega_tree->Branch("om_ld_chi2topo",&om_ld_chi2topo,"om_ld_chi2topo/F");
   omega_tree->Branch("om_ld_chi2ndf",&om_ld_chi2ndf,"om_ld_chi2ndf/F");
   omega_tree->Branch("om_ld_ldl",&om_ld_ldl,"om_ld_ldl/F");
   omega_tree->Branch("om_ld_l",&om_ld_l,"om_ld_l/F");

   omega_tree->Branch("om_l",&om_l,"om_l/F");
   omega_tree->Branch("om_dl",&om_dl,"om_dl/F");

   omega_tree->Branch("chi2primary_om_proton",&chi2primary_om_proton,"chi2primary_om_proton/F");
   omega_tree->Branch("chi2primary_om_pi",&chi2primary_om_pi,"chi2primary_om_pi/F");
   omega_tree->Branch("chi2primary_om_bach",&chi2primary_om_bach,"chi2primary_om_bach/F");
   omega_tree->Branch("chi2primary_om_ld",&chi2primary_om_ld,"chi2primary_om_ld/F");

   omega_tree->Branch("nhits_om_proton",&nhits_om_proton,"nhits_om_proton/I");
   omega_tree->Branch("nhits_om_pi",&nhits_om_pi,"nhits_om_pi/I");
   omega_tree->Branch("nhits_om_bach",&nhits_om_bach,"nhits_om_bach/I");
   omega_tree->Branch("dedx_om_proton",&dedx_om_proton,"dedx_om_proton/F");
   omega_tree->Branch("dedx_om_pi",&dedx_om_pi,"dedx_om_pi/F");
   omega_tree->Branch("dedx_om_bach",&dedx_om_bach,"dedx_om_bach/F");

   omega_tree->Branch("bbachpx",&bbachpx,"bbachpx/F");
   omega_tree->Branch("bbachpy",&bbachpy,"bbachpy/F");
   omega_tree->Branch("bbachpz",&bbachpz,"bbachpz/F");

   omega_tree->Branch("bprotonpx",&bprotonpx,"bprotonpx/F");
   omega_tree->Branch("bprotonpy",&bprotonpy,"bprotonpy/F");
   omega_tree->Branch("bprotonpz",&bprotonpz,"bprotonpz/F");

   omega_tree->Branch("bpionpx",&bpionpx,"bpionpx/F");
   omega_tree->Branch("bpionpy",&bpionpy,"bpionpy/F");
   omega_tree->Branch("bpionpz",&bpionpz,"bpionpz/F");

   omega_tree->Branch("notbadrun",&notbadrun,"notbadrun/I");

   omega_tree->Branch("bmcpx",&bmcpx,"bmcpx/F");
   omega_tree->Branch("bmcpy",&bmcpy,"bmcpy/F");
   omega_tree->Branch("bmcpz",&bmcpz,"bmcpz/F");
   omega_tree->Branch("bismc",&bismc,"bismc/I");


   omega_mc_tree = new TTree("omega_mc_tree","ana_tree");
   omega_mc_tree->Branch("brunid",&brunid,"brunid/I");
   omega_mc_tree->Branch("brefmult",&brefmult,"brefmult/I");
   omega_mc_tree->Branch("bmcparticleid",&bmcparticleid,"bmcparticleid/I");
   omega_mc_tree->Branch("bmcrawpx",&bmcrawpx,"bmcrawpx/F");
   omega_mc_tree->Branch("bmcrawpy",&bmcrawpy,"bmcrawpy/F");
   omega_mc_tree->Branch("bmcrawpz",&bmcrawpz,"bmcrawpz/F");

   omega_mc_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
   omega_mc_tree->Branch("reweight",&reweight,"reweight/D");
   omega_mc_tree->Branch("cent9",&cent9,"cent9/I");

/*
  omega_tree->Branch("bbachid",&bbachid,"bbachid/I");
  omega_tree->Branch("bbachpx",&bbachpx,"bbachpx/F");
  omega_tree->Branch("bbachpy",&bbachpy,"bbachpy/F");
  omega_tree->Branch("bbachpz",&bbachpz,"bbachpz/F");
  omega_tree->Branch("bbachmass",&bbachmass,"bbachmass/F");

  omega_tree->Branch("bpionid",&bpionid,"bpionid/I");
  omega_tree->Branch("bpionpx",&bpionpx,"bpionpx/F");
  omega_tree->Branch("bpionpy",&bpionpy,"bpionpy/F");
  omega_tree->Branch("bpionpz",&bpionpz,"bpionpz/F");
  omega_tree->Branch("bpionmass",&bpionmass,"bpionmass/F");

  omega_tree->Branch("bprotonid",&bprotonid,"bprotonid/I");
  omega_tree->Branch("bprotonpx",&bprotonpx,"bprotonpx/F");
  omega_tree->Branch("bprotonpy",&bprotonpy,"bprotonpy/F");
  omega_tree->Branch("bprotonpz",&bprotonpz,"bprotonpz/F");
  omega_tree->Branch("bprotonmass",&bprotonmass,"bprotonmass/F");
*/

   htriton_tree = new TTree("htriton_tree","ana_tree");
   htriton_tree->Branch("brunid",&brunid,"brunid/I");
   htriton_tree->Branch("beventid",&beventid,"beventid/I");
 
   htriton_tree->Branch("bVz",&bVz,"bVz/F");
     // omega_tree->Branch("bVr",&bVr,"bVr/F");
   htriton_tree->Branch("bVzerr",&bVzerr,"bVzerr/F");
     // omega_tree->Branch("bVrerr",&bVrerr,"bVrerr/F"); 
   htriton_tree->Branch("brefmult",&brefmult,"brefmult/I");
   htriton_tree->Branch("btofmult",&btofmult,"btofmult/I");
   htriton_tree->Branch("countrefmult",&countrefmult,"countrefmult/I");
 
   htriton_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
   htriton_tree->Branch("reweight",&reweight,"reweight/D");
   htriton_tree->Branch("cent9",&cent9,"cent9/I");
 
   htriton_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
   htriton_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");
                     // omega_tree->Branch("bx",&bx,"bx/F");
                     // omega_tree->Branch("by",&by,"by/F");
                     // omega_tree->Branch("bz",&bz,"bz/F");
   htriton_tree->Branch("bpx",&bpx,"bpx/F");
   htriton_tree->Branch("bpy",&bpy,"bpy/F");
   htriton_tree->Branch("bpz",&bpz,"bpz/F");
   htriton_tree->Branch("bpl",&bpl,"bpl/F");

   htriton_tree->Branch("chi2primary_he",&chi2primary_he,"chi2primary_he/F");
   htriton_tree->Branch("chi2primary_pi",&chi2primary_pi,"chi2primary_pi/F");
   htriton_tree->Branch("dca_he",&dca_he,"dca_he/F");
   htriton_tree->Branch("dca_pi",&dca_pi,"dca_pi/F");

   htriton_tree->Branch("ht_ldl",&ht_ldl,"ht_ldl/F");
   htriton_tree->Branch("ht_dl",&ht_dl,"ht_dl/F");
   htriton_tree->Branch("ht_l",&ht_l,"ht_l/F");
   htriton_tree->Branch("ht_chi2topo",&ht_chi2topo,"ht_chi2topo/F");
   htriton_tree->Branch("ht_chi2ndf",&ht_chi2ndf,"ht_chi2ndf/F");

   htriton_tree->Branch("nhits_he",&nhits_he,"nhits_he/I");
   htriton_tree->Branch("nhits_pi",&nhits_pi,"nhits_pi/I");

   htriton_tree->Branch("dedx_he",&dedx_he,"dedx_he/F");
   htriton_tree->Branch("dedx_pi",&dedx_pi,"dedx_pi/F");

   htriton_tree->Branch("ht_chi2",&ht_chi2,"ht_chi2/F");
   htriton_tree->Branch("ht_NDF",&ht_NDF,"ht_NDF/F");

   htriton_tree->Branch("ht_bdfvtx",&ht_bdfvtx,"ht_bdfvtx/F");
   htriton_tree->Branch("ht_bdfvtx2",&ht_bdfvtx2,"ht_bdfvtx2/F");
   htriton_tree->Branch("ht_lifetime",&ht_lifetime,"ht_lifetime/F");

   htriton_tree->Branch("px_pi",&px_pi,"px_pi/F");
   htriton_tree->Branch("py_pi",&py_pi,"py_pi/F");
   htriton_tree->Branch("pz_pi",&pz_pi,"pz_pi/F");
   htriton_tree->Branch("px_he",&px_he,"px_he/F");
   htriton_tree->Branch("py_he",&py_he,"py_he/F");
   htriton_tree->Branch("pz_he",&pz_he,"pz_he/F");

   htriton_tree->Branch("bismc",&bismc,"bismc/I");
   htriton_tree->Branch("bmcpx",&bmcpx,"bmcpx/F");
   htriton_tree->Branch("bmcpy",&bmcpy,"bmcpy/F");
   htriton_tree->Branch("bmcpz",&bmcpz,"bmcpz/F");
   htriton_tree->Branch("bmcl",&bmcl,"bmcl/F");
   htriton_tree->Branch("bmcpl",&bmcpl,"bmcpl/F");

   if(fFlowAnalysis){
   htriton_tree->Branch("fCentrality",&fCentrality,"fCentrality/I");
   htriton_tree->Branch("psi_1_EPD_0",&psi_1_EPD_0,"psi_1_EPD_0/D");
   htriton_tree->Branch("psi_1_EPD_1",&psi_1_EPD_1,"psi_1_EPD_1/D");
   htriton_tree->Branch("psi_1_EPD_2",&psi_1_EPD_2,"psi_1_EPD_2/D");
   htriton_tree->Branch("psi_1_EPD_3",&psi_1_EPD_3,"psi_1_EPD_3/D");
   htriton_tree->Branch("psi_1_EPD_4",&psi_1_EPD_4,"psi_1_EPD_4/D");
   htriton_tree->Branch("psi_1_EPD_5",&psi_1_EPD_5,"psi_1_EPD_5/D");
   htriton_tree->Branch("psi_1_EPD_6",&psi_1_EPD_6,"psi_1_EPD_6/D");
   htriton_tree->Branch("psi_1_EPD",&psi_1_EPD,"psi_1_EPD/D");
   htriton_tree->Branch("gweight",&gweight,"gweight/F");
   }

   he3_tree = new TTree("he3_tree","ana_tree");
   he3_tree->Branch("brunid",&brunid,"brunid/I");
   he3_tree->Branch("bVz",&bVz,"bVz/F");
 
   he3_tree->Branch("brefmult",&brefmult,"brefmult/I");
   he3_tree->Branch("btofmult",&btofmult,"btofmult/I");

   he3_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
   he3_tree->Branch("reweight",&reweight,"reweight/D");
   he3_tree->Branch("cent9",&cent9,"cent9/I");

   he3_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
   he3_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");
 
   he3_tree->Branch("bpx",&bpx,"bpx/F");
   he3_tree->Branch("bpy",&bpy,"bpy/F");
   he3_tree->Branch("bpz",&bpz,"bpz/F");
 
   he3_tree->Branch("chi2primary_he",&chi2primary_he,"chi2primary_he/F");
   he3_tree->Branch("bnhits",&bnhits,"bnhits/I");
   he3_tree->Branch("bdedx",&bdedx,"bdedx/F");
   he3_tree->Branch("bdca",&bdca,"bdca/F");

   he3_tree->Branch("bmcpx",&bmcpx,"bmcpx/F");
   he3_tree->Branch("bmcpy",&bmcpy,"bmcpy/F");
   he3_tree->Branch("bmcpz",&bmcpz,"bmcpz/F");
   he3_tree->Branch("bismc",&bismc,"bismc/I");
 
   h4lambda_tree = new TTree("h4lambda_tree","ana_tree");
   h4lambda_tree->Branch("brunid",&brunid,"brunid/I");
   h4lambda_tree->Branch("beventid",&beventid,"beventid/I");
   h4lambda_tree->Branch("bVz",&bVz,"bVz/F");
   h4lambda_tree->Branch("brefmult",&brefmult,"brefmult/I");
   h4lambda_tree->Branch("btofmult",&btofmult,"btofmult/I");
   h4lambda_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
   h4lambda_tree->Branch("reweight",&reweight,"reweight/D");
   h4lambda_tree->Branch("cent9",&cent9,"cent9/I");
   h4lambda_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
   h4lambda_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");
   h4lambda_tree->Branch("bpx",&bpx,"bpx/F");
   h4lambda_tree->Branch("bpy",&bpy,"bpy/F");
   h4lambda_tree->Branch("bpz",&bpz,"bpz/F");
   h4lambda_tree->Branch("chi2primary_h4",&chi2primary_h4,"chi2primary_h4/F");
   h4lambda_tree->Branch("chi2primary_pi",&chi2primary_pi,"chi2primary_pi/F");
   h4lambda_tree->Branch("hl_ldl",&hl_ldl,"hl_ldl/F");
   h4lambda_tree->Branch("hl_dl",&hl_dl,"hl_dl/F");
   h4lambda_tree->Branch("hl_l",&hl_l,"hl_l/F");
   h4lambda_tree->Branch("hl_chi2topo",&hl_chi2topo,"hl_chi2topo/F");
   h4lambda_tree->Branch("hl_chi2ndf",&hl_chi2ndf,"hl_chi2ndf/F");
   h4lambda_tree->Branch("hl_chi2",&hl_chi2,"hl_chi2/F");
   h4lambda_tree->Branch("hl_NDF",&hl_NDF,"hl_NDF/F");


   he3_mc_tree = new TTree("he3_mc_tree","ana_tree");
   he3_mc_tree->Branch("brunid",&brunid,"brunid/I");
   he3_mc_tree->Branch("brefmult",&brefmult,"brefmult/I");
   he3_mc_tree->Branch("bmcparticleid",&bmcparticleid,"bmcparticleid/I");
   he3_mc_tree->Branch("bmcrawpx",&bmcrawpx,"bmcrawpx/F");
   he3_mc_tree->Branch("bmcrawpy",&bmcrawpy,"bmcrawpy/F");
   he3_mc_tree->Branch("bmcrawpz",&bmcrawpz,"bmcrawpz/F");
   he3_mc_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
   he3_mc_tree->Branch("reweight",&reweight,"reweight/D");
   he3_mc_tree->Branch("cent9",&cent9,"cent9/I");

   h_tree = new TTree("h_tree","ana_tree");
   h_tree->Branch("brunid",&brunid,"brunid/I");
   h_tree->Branch("beventid",&beventid,"beventid/I");
   h_tree->Branch("bVz",&bVz,"bVz/F");
   h_tree->Branch("brefmult",&brefmult,"brefmult/I");
   h_tree->Branch("btofmult",&btofmult,"btofmult/I");
   h_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
   h_tree->Branch("reweight",&reweight,"reweight/D");
   h_tree->Branch("cent9",&cent9,"cent9/I");
   h_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
   h_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");
   h_tree->Branch("bpx",&bpx,"bpx/F");
   h_tree->Branch("bpy",&bpy,"bpy/F");
   h_tree->Branch("bpz",&bpz,"bpz/F");
   h_tree->Branch("hl_ldl",&hl_ldl,"hl_ldl/F");
   h_tree->Branch("hl_dl",&hl_dl,"hl_dl/F");
   h_tree->Branch("hl_l",&hl_l,"hl_l/F");
   h_tree->Branch("hl_chi2topo",&hl_chi2topo,"hl_chi2topo/F");
   h_tree->Branch("hl_chi2ndf",&hl_chi2ndf,"hl_chi2ndf/F");
   hvtx      = new TH1F("hvtx",    "Vz;Vz(cm);Counts",400,0.,400);
   hvtxgood  = new TH1F("hvtxgood","Vz;Vz(cm);Counts",400,0.,400);
   hrefmult  = new TH1F("hrefmult", "refmult; hrefmult; N_{evt}", 600,0,600);
   wrefmult  = new TH1F("wrefmult", "refmult; wrefmult; N_{evt}", 600,0,600);

   return kStOK;
}
//________________________________________________________________________________
Int_t StKFParticleAnalysisMaker::InitRun(Int_t runumber) 
{
//   assert(StPicoDstMaker::instance());
//   if (StPicoDstMaker::instance()->IOMode() == StPicoDstMaker::ioRead) {
    //TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO Ask Yuri
//     StPicoDstMaker::instance()->SetStatus("*",0);
//     const Char_t *ActiveBranches[] = {
//       "MuEvent"
//       ,"PrimaryVertices"
//       ,"PrimaryTracks"
//       ,"GlobalTracks"
//       ,"StStMuMcVertex"
//       ,"StStMuMcTrack"
//       ,"CovPrimTrack"
//       ,"CovGlobTrack"
//       ,"StStMuMcVertex"
//       ,"StStMuMcTrack"
//       ,"KFTracks"
//       ,"KFVertices"
//       ,"StBTofHit"
//       ,"StBTofHeader"
//     }; 
//     Int_t Nb = sizeof(ActiveBranches)/sizeof(Char_t *);
//     for (Int_t i = 0; i < Nb; i++) StPicoDstMaker::instance()->SetStatus(ActiveBranches[i],1); // Set Active braches
//   }
  return StMaker::InitRun(runumber);
}
//_____________________________________________________________________________
void StKFParticleAnalysisMaker::PrintMem(const Char_t *opt)
{
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  cout << opt 
       << "\tMemory : Total = " << info.fMemTotal 
       << "\tUsed = " << info.fMemUsed
       << "\tFree = " << info.fMemFree
       << "\tSwap Total = " << info.fSwapTotal
       << "\tUsed = " << info.fSwapUsed
       << "\tFree = " << info.fSwapFree << endl;
}
//_____________________________________________________________________________
void StKFParticleAnalysisMaker::BookVertexPlots()
{
  TDirectory *dirs[2] = {0};
  dirs[0] = TDirectory::CurrentDirectory(); assert(dirs[0]);
  dirs[0]->cd();
  if (! dirs[0]->GetDirectory("Particles")) {
    dirs[0]->mkdir("Particles");
  }
  dirs[1] = dirs[0]->GetDirectory("Particles"); assert(dirs[1]);
  dirs[1]->cd();
  PrintMem(dirs[1]->GetPath());
  
  fStKFParticleInterface = new StKFParticleInterface;
  bool storeMCHistograms = false;
  if(!fIsPicoAnalysis && fProcessSignal) storeMCHistograms = true;
  fStKFParticlePerformanceInterface = new StKFParticlePerformanceInterface(fStKFParticleInterface->GetTopoReconstructor(), storeMCHistograms);
  dirs[0]->cd();
  PrintMem(dirs[1]->GetPath());
}
//_____________________________________________________________________________
Int_t StKFParticleAnalysisMaker::Make()
{  
  if(fIsPicoAnalysis)
  {
    fPicoDst = StPicoDst::instance();
    if(!fPicoDst) return kStOK;
  }
  else
  {  
    fMuDst = StMuDst::instance();
    if(!fMuDst) return kStOK;
    else { if(StMuDst::instance()->numberOfPrimaryVertices() == 0 ) return kStOK; }
  }

  //_fill_lambda_tree = false;
  _fill_lambda_tree = true;
 
  //_fill_he3_tree = true;
  _fill_he3_tree = false;
 
  //find max global track index
  int maxGBTrackIndex = -1;
  if(fIsPicoAnalysis)
  {
    for(unsigned int iTrack = 0; iTrack < fPicoDst->numberOfTracks(); iTrack++) 
    {
      StPicoTrack *gTrack = fPicoDst->track(iTrack);
      if (! gTrack) continue;
      int index = gTrack->id();
      if(index > maxGBTrackIndex)
        maxGBTrackIndex = index;
    }
  }
  else
  {
    for(unsigned int iTrack = 0; iTrack < fMuDst->numberOfGlobalTracks(); iTrack++) 
    {
      StMuTrack *gTrack = fMuDst->globalTracks(iTrack);
      if (! gTrack) continue;
      int index = gTrack->id();
      if(index > maxGBTrackIndex)
        maxGBTrackIndex = index;
    }
  }
  vector<KFMCTrack> mcTracks(0);
  vector<int> mcIndices(maxGBTrackIndex+1);
  for(unsigned int iIndex=0; iIndex<mcIndices.size(); iIndex++)
    mcIndices[iIndex] = -1;
  
//   fStKFParticleInterface->SetTriggerMode();
//   fStKFParticleInterface->SetSoftKaonPIDMode();
//   fStKFParticleInterface->SetSoftTofPidMode();
//   fStKFParticleInterface->SetChiPrimaryCut(10);
//   
//   fStKFParticleInterface->SetPtCutCharm(0.5);
//   fStKFParticleInterface->SetChiPrimaryCutCharm(8);
//   fStKFParticleInterface->SetLdLCutCharmManybodyDecays(3);
//   fStKFParticleInterface->SetChi2TopoCutCharmManybodyDecays(10);
//   fStKFParticleInterface->SetChi2CutCharmManybodyDecays(3);
//   fStKFParticleInterface->SetLdLCutCharm2D(3);
//   fStKFParticleInterface->SetChi2TopoCutCharm2D(10);
//   fStKFParticleInterface->SetChi2CutCharm2D(3);
 
  //work out the primary track count (for FXT)
  if(fIsPicoAnalysis)
  {
    countrefmult = 0;
    for(unsigned int iTrack = 0; iTrack < fPicoDst->numberOfTracks(); iTrack++)
    {
      const StPicoTrack *ptrk = (StPicoTrack*)fPicoDst->track(iTrack);
      if(! ptrk) continue;
      if(!ptrk->isPrimary())  continue;  // now selecting primary tracks
      const float dca = ptrk->gDCA( fPicoDst->event()->primaryVertex() ).Mag();
      const int nHitsFit = ptrk->nHitsFit();
      const int nHitsPoss = ptrk->nHitsMax();
      const float quality = (float)nHitsFit/(float)nHitsPoss;
      if( fabs(dca)>3.0 ) continue;
      if( nHitsFit < 15 )  continue;;
      if( quality < 0.52 )  continue;
      countrefmult++;
    }
  }
  else
  {
    countrefmult = 0;
    float bestRank=-1000000;
    int bestPV=-1;
    for(unsigned int iPV=0; iPV<fMuDst->numberOfPrimaryVertices(); iPV++)
    {
      StMuPrimaryVertex *Vtx = fMuDst->primaryVertex(iPV);
      if(!Vtx) continue;
      if (bestRank < Vtx->ranking()) {
      bestRank = Vtx->ranking();
      bestPV = iPV;
      }
      else continue;
    }
    if(bestPV!=-1){
    for(unsigned int iTrack = 0; iTrack < fMuDst->numberOfPrimaryTracks(); iTrack++)
      {
      StMuTrack *gTrack = fMuDst->primaryTracks(iTrack);
      if (! gTrack) continue;
      const int bnHitsFit = gTrack->nHitsFit();
      const int bnHitsPoss = gTrack->nHitsPoss();
      const float bquality = (float)bnHitsFit/(float)bnHitsPoss;
      const float bdca = gTrack->dcaGlobal(bestPV).mag();
      if ( bnHitsFit < 15 ) continue;
      if ( bquality < 0.52 ) continue;
      if ( fabs(bdca) > 3.0 ) continue;
      countrefmult++;
      }
    }
  } 
  //end work out the primary track count (for FXT)
 
  vector<int> triggeredTracks;
  bool isGoodEvent = false;
  
  //Process the event
  if(maxGBTrackIndex > 0)
    fStKFParticleInterface->ResizeTrackPidVectors(maxGBTrackIndex+1);
  if(fIsPicoAnalysis)
    isGoodEvent = fStKFParticleInterface->ProcessEvent(fPicoDst, triggeredTracks);
  else//b
//no embedding
//    isGoodEvent = fStKFParticleInterface->ProcessEvent(fMuDst, triggeredTracks, fProcessSignal);
//embedding
    isGoodEvent = fStKFParticleInterface->ProcessEvent(fMuDst, mcTracks, mcIndices, fProcessSignal);

//   bool openCharmTrigger = false;
//   if(isGoodEvent) openCharmTrigger =  fStKFParticleInterface->OpenCharmTrigger();
//   fStKFParticleInterface->OpenCharmTriggerCompression(triggeredTracks.size(), fPicoDst->numberOfTracks(), openCharmTrigger);
  //collect histograms
 

//events removed at this level include no track events and 
//whatever is in Process event skip level

    bool bWriteTree = true;
//if(bWriteTree){ 
       if(fIsPicoAnalysis){
	beventid = fPicoDst->event()->eventId();
	brunid = fPicoDst->event()->runId();

  const TVector3 picoPV = fPicoDst->event()->primaryVertex();
  const TVector3 picoPVError = fPicoDst->event()->primaryVertexError();

	bVx = picoPV.x();
	bVy = picoPV.y();
        bVr = sqrt(picoPV.x()*picoPV.x()+picoPV.y()*picoPV.y());
        bVz = picoPV.z(); 
        bVrerr = sqrt(picoPVError.x()*picoPVError.x()+picoPVError.y()*picoPVError.y());
	bVxerr = picoPVError.x();
	bVyerr = picoPVError.y();
        bVzerr = picoPVError.z();

	}
       else{
	beventid = fMuDst->event()->eventId();
	brunid = fMuDst->event()->runId();

//cout<<"brunid:"<<brunid<<endl;
  float bestRank=-1000000;
  int bestPV=-1;
  for(unsigned int iPV=0; iPV<fMuDst->numberOfPrimaryVertices(); iPV++)
  {
    StMuPrimaryVertex *Vtx = fMuDst->primaryVertex(iPV);
    if(!Vtx) continue;
    if (bestRank < Vtx->ranking()) {
      bestRank = Vtx->ranking();
      bestPV = iPV;
    }
    else continue;

    bVr = sqrt(Vtx->position().x()*Vtx->position().x()+Vtx->position().y()*Vtx->position().y());
    bVx = Vtx->position().x();
    bVy = Vtx->position().y();
    bVz = Vtx->position().z();
    bVrerr = sqrt(Vtx->posError().x()*Vtx->posError().x()+Vtx->posError().y()*Vtx->posError().y());
    bVxerr = Vtx->posError().x();
    bVyerr = Vtx->posError().y();
    bVzerr = Vtx->posError().z();
	}

}

	hvtx->Fill(bVz);
//}//if writetree


//cut on refmult and tofmult, trigger, to cut on good event. this cut is only for real analysis!! 
//when fStoreTmvaNTuples is true, there is no need to cut on these variables
int crefmult;
int ctofmult;
   if(fIsPicoAnalysis)
 {
   //brunid   = fPicoDst->event()->runId();
   crefmult = fPicoDst->event()->refMult();
   ctofmult = fPicoDst->event()->btofTrayMultiplicity();
 }
 else
 {
   //brunid   = fMuDst->event()->runId();
   crefmult = fMuDst->event()->refMult();
   ctofmult = fMuDst->event()->btofTrayMultiplicity();
 }
  
 //if(crefmult>(0.42*ctofmult+16) || crefmult<(0.02*ctofmult-5)) isGoodEvent = false;

  int trigger=0;
  notbadrun=0;

   if(fIsPicoAnalysis){
	if(fsnn==27){
  		if(fPicoDst->event()->isTrigger(610001)) trigger += 1;
  		if(fPicoDst->event()->isTrigger(610011)) trigger += 10;
		if(fPicoDst->event()->isTrigger(610021)) trigger += 100;
  		if(fPicoDst->event()->isTrigger(610031)) trigger += 1000;
  		if(fPicoDst->event()->isTrigger(610041)) trigger += 10000;
  		if(fPicoDst->event()->isTrigger(610051)) trigger += 100000;
	}else if(fsnn==3){
  		if(fPicoDst->event()->isTrigger(620052)) trigger += 1;
	}else{
  		trigger += 1;
	}
    }
    else{
	trigger+=1; //auto pass trigger because we mudst already has trigger seltcion
    }

   if(fsnn==27){
   if(brunid==19130085) notbadrun+=1;
   if(brunid==19131009) notbadrun+=1;
   if(brunid==19131010) notbadrun+=1;
   if(brunid==19131012) notbadrun+=1;
   if(brunid==19132063) notbadrun+=1;
   if(brunid==19133009) notbadrun+=1;
   if(brunid==19133010) notbadrun+=1;
   if(brunid==19133012) notbadrun+=1;
   if(brunid==19133013) notbadrun+=1;
   if(brunid==19133014) notbadrun+=1;
   if(brunid==19133018) notbadrun+=1;
   if(brunid==19134010) notbadrun+=1;
   if(brunid==19134011) notbadrun+=1;
   if(brunid==19135011) notbadrun+=1;
   if(brunid==19135013) notbadrun+=1;
   if(brunid==19135014) notbadrun+=1;
   if(brunid==19136016) notbadrun+=1;
   if(brunid==19137003) notbadrun+=1;
   if(brunid==19137022) notbadrun+=1;
   if(brunid==19137047) notbadrun+=1;
   if(brunid==19137050) notbadrun+=1;
   if(brunid==19137051) notbadrun+=1;
   if(brunid==19137052) notbadrun+=1;
   if(brunid==19137053) notbadrun+=1;
   if(brunid==19137056) notbadrun+=1;
   if(brunid==19137057) notbadrun+=1;
   if(brunid==19138008) notbadrun+=1;
   if(brunid==19138009) notbadrun+=1;
   if(brunid==19138014) notbadrun+=1;
   if(brunid==19139022) notbadrun+=1;
   if(brunid==19139023) notbadrun+=1;
   if(brunid==19139024) notbadrun+=1;
   if(brunid==19139026) notbadrun+=1;
   if(brunid==19139027) notbadrun+=1;
   if(brunid==19139028) notbadrun+=1;
   if(brunid==19139032) notbadrun+=1;
   if(brunid==19139033) notbadrun+=1;
   if(brunid==19139034) notbadrun+=1;
   if(brunid==19139037) notbadrun+=1;
   if(brunid==19140009) notbadrun+=1;
   if(brunid==19140014) notbadrun+=1;
   if(brunid==19141008) notbadrun+=1;
   if(brunid==19142005) notbadrun+=1;
   if(brunid==19142048) notbadrun+=1;
   if(brunid==19143008) notbadrun+=1;
   if(brunid==19143009) notbadrun+=1;
   if(brunid==19143010) notbadrun+=1;
   if(brunid==19143011) notbadrun+=1;
   if(brunid==19143012) notbadrun+=1;
   if(brunid==19143013) notbadrun+=1;
   if(brunid==19143014) notbadrun+=1;
   if(brunid==19143015) notbadrun+=1;
   if(brunid==19143016) notbadrun+=1;
   if(brunid==19143017) notbadrun+=1;
   if(brunid==19146016) notbadrun+=1;
   if(brunid==19147007) notbadrun+=1;
   if(brunid==19147008) notbadrun+=1;
   if(brunid==19147009) notbadrun+=1;
   if(brunid==19147010) notbadrun+=1;
   if(brunid==19147014) notbadrun+=1;
   if(brunid==19147015) notbadrun+=1;
   if(brunid==19147016) notbadrun+=1;
   if(brunid==19156002) notbadrun+=1;
   if(brunid==19156032) notbadrun+=1;
   if(brunid==19156044) notbadrun+=1;
   if(brunid==19156045) notbadrun+=1;
   if(brunid==19156046) notbadrun+=1;
   if(brunid==19157013) notbadrun+=1;
   if(brunid==19157018) notbadrun+=1;
   if(brunid==19158003) notbadrun+=1;
   if(brunid==19158007) notbadrun+=1;
   if(brunid==19158009) notbadrun+=1;
   if(brunid==19158010) notbadrun+=1;
   if(brunid==19158011) notbadrun+=1;
   if(brunid==19158013) notbadrun+=1;
   if(brunid==19158014) notbadrun+=1;
   if(brunid==19158015) notbadrun+=1;
   if(brunid==19158017) notbadrun+=1;
   if(brunid==19158018) notbadrun+=1;
   if(brunid==19158019) notbadrun+=1;
   if(brunid==19160018) notbadrun+=1;
   if(brunid==19162002) notbadrun+=1;
   if(brunid==19162005) notbadrun+=1;
   if(brunid==19165015) notbadrun+=1;
   if(brunid==19165020) notbadrun+=1;
   if(brunid==19165021) notbadrun+=1;
   if(brunid==19167042) notbadrun+=1;
   }
   if(fsnn==3){
/*
if(brunid==19160032) notbadrun+=1;
if(brunid==19160033) notbadrun+=1;
if(brunid==19160034) notbadrun+=1;
if(brunid==19160035) notbadrun+=1;
if(brunid==19160036) notbadrun+=1;
if(brunid==19160037) notbadrun+=1;
if(brunid==19160038) notbadrun+=1;
if(brunid==19160039) notbadrun+=1;
if(brunid==19160040) notbadrun+=1;
if(brunid==19160041) notbadrun+=1;
if(brunid==19160042) notbadrun+=1;
if(brunid==19160043) notbadrun+=1;
if(brunid==19160044) notbadrun+=1;
if(brunid==19161001) notbadrun+=1;
if(brunid==19161020) notbadrun+=1;
if(brunid==19161021) notbadrun+=1;
if(brunid==19161022) notbadrun+=1;
if(brunid==19161023) notbadrun+=1;
if(brunid==19161024) notbadrun+=1;
if(brunid==19161025) notbadrun+=1;
if(brunid==19161026) notbadrun+=1;
if(brunid==19161027) notbadrun+=1;
if(brunid==19161028) notbadrun+=1;
if(brunid==19161029) notbadrun+=1;
if(brunid==19161030) notbadrun+=1;
if(brunid==19161034) notbadrun+=1;
if(brunid==19161035) notbadrun+=1;
if(brunid==19161036) notbadrun+=1;
if(brunid==19161037) notbadrun+=1;
if(brunid==19164001) notbadrun+=1;
if(brunid==19164022) notbadrun+=1;
if(brunid==19164024) notbadrun+=1;
if(brunid==19167053) notbadrun+=1;
*/
   if(brunid==19151029) notbadrun+=1;
   if(brunid==19151045) notbadrun+=1;
   if(brunid==19152001) notbadrun+=1;
   if(brunid==19152034) notbadrun+=1;
   if(brunid==19152038) notbadrun+=1;
   if(brunid==19152045) notbadrun+=1;
   if(brunid==19152078) notbadrun+=1;
   if(brunid==19153023) notbadrun+=1;
   if(brunid==19153032) notbadrun+=1;
   if(brunid==19154051) notbadrun+=1;
   if(brunid==19154012) notbadrun+=1;
   if(brunid==19154013) notbadrun+=1;
   if(brunid==19154014) notbadrun+=1;
   if(brunid==19154015) notbadrun+=1;
   if(brunid==19154016) notbadrun+=1;
   if(brunid==19154017) notbadrun+=1;
   if(brunid==19154018) notbadrun+=1;
   if(brunid==19154019) notbadrun+=1;
   if(brunid==19154020) notbadrun+=1;
   if(brunid==19154021) notbadrun+=1;
   if(brunid==19154022) notbadrun+=1;
   if(brunid==19154023) notbadrun+=1;
   if(brunid==19154024) notbadrun+=1;
   if(brunid==19154026) notbadrun+=1;
   }

  //cout<<"trigger:"<<trigger<<endl;
   if(trigger==0) isGoodEvent = false;
   if(notbadrun!=0) isGoodEvent = false;


//if(isGoodEvent){
//        hvtxgood->Fill(bVz);
//        hrefmult->Fill(crefmult);
//}
 
   reweight = 0;
   refmultcor = 0;

  if(isGoodEvent)
//  if(1)
  {
    int centralityBin = -1;
    float centralityWeight = 0.;
    
    if(fRunCentralityAnalysis)
    {
	if(fIsPicoAnalysis){

      fRefmultCorrUtil->init(fPicoDst->event()->runId());
      if(! (fRefmultCorrUtil->isBadRun(fPicoDst->event()->runId())) )
      {
        //fRefmultCorrUtil->initEvent(fPicoDst->event()->grefMult(), fPicoDst->event()->primaryVertex().z(), fPicoDst->event()->ZDCx()) ;
        fRefmultCorrUtil->initEvent(fPicoDst->event()->refMult(), fPicoDst->event()->primaryVertex().z(), fPicoDst->event()->ZDCx()) ;
        centralityBin = fRefmultCorrUtil->getCentralityBin9();
        centralityWeight = fRefmultCorrUtil->getWeight();
        
        cent9 = centralityBin;
        reweight = centralityWeight; 
      }else{
        isGoodEvent = false;
      }

      Bool_t isPileUpEvt = !fRefmultCorrUtil->passnTofMatchRefmultCut(1.*fPicoDst->event()->refMult(), 1.*fPicoDst->event()->nBTOFMatch()); //reject pileup events
      if(isPileUpEvt) isGoodEvent = false;

        //refmultCor = fRefmultCorrUtil->getRefMultCorr();
        refmultcor = fRefmultCorrUtil->getRefMultCorr();
        //refmultcor = refmultCor;
      }else{

	fRefmultCorrUtil->init(fMuDst->event()->runId());
	if(! (fRefmultCorrUtil->isBadRun(fMuDst->event()->runId())) )
        {
	fRefmultCorrUtil->initEvent(fMuDst->event()->refMult(), bVz, fMuDst->event()->runInfo().zdcCoincidenceRate()) ;
	//cout<<"zdcrate:"<<fMuDst->event()->runInfo().zdcCoincidenceRate()<<endl;
   
	centralityBin = fRefmultCorrUtil->getCentralityBin9();
        centralityWeight = fRefmultCorrUtil->getWeight();
	cent9 = centralityBin;
        reweight = centralityWeight;
	//skipt rejecting pileup for mudst
	refmultcor = fRefmultCorrUtil->getRefMultCorr();
	}else{
          isGoodEvent = false;
        }
      }
    }
    
    if(fTMVAselection){
	cout<<"TMVA SELECTION!!!!!!!!!!!!! SOMETHING IS WRONG-----------------"<<endl;
    }
    
    if(fTMVAselection)//check this bool
    {
      for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
      {
        KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];
              
        for(int iReader=0; iReader<fNNTuples; iReader++)
        {
          if( abs(particle.GetPDG()) == fNTuplePDG[iReader] )
          {
            GetParticleParameters(iReader, particle);
            
            const int iTMVACentralityBin = GetTMVACentralityBin(iReader, centralityBin);
            const int iTMVAPtBin = GetTMVAPtBin(iReader, particle.GetPt());
            
            if(iTMVACentralityBin<0 || iTMVAPtBin<0) 
            {
              fStKFParticleInterface->RemoveParticle(iParticle);
              continue;
            }
            
            if(fTMVAReader[iReader][iTMVACentralityBin][iTMVAPtBin]->EvaluateMVA("BDT") < fTMVACut[iReader][iTMVACentralityBin][iTMVAPtBin])
              fStKFParticleInterface->RemoveParticle(iParticle);
           /* 
            if(fAnalyseDsPhiPi && abs(fStKFParticleInterface->GetParticles()[iParticle].GetPDG()) == 431)
            {              
              KFParticle phi;
              if(particle.GetPDG() == 431)
                phi += fStKFParticleInterface->GetParticles()[particle.DaughterIds()[0]];
              else
                phi += fStKFParticleInterface->GetParticles()[particle.DaughterIds()[1]];
              phi += fStKFParticleInterface->GetParticles()[particle.DaughterIds()[2]];
              float mass = 0.f, dmass = 0.f;
              phi.GetMass(mass, dmass);
              if( fabs(mass - 1.01946) > 0.015)
                fStKFParticleInterface->RemoveParticle(iParticle);
            }
          */
          }
        }
      }      
    }
    
    //clean H3L, H4L, Ln, Lnn
    for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
    {
      KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];
      if( abs(particle.GetPDG())==3003 || abs(particle.GetPDG())==3103 || abs(particle.GetPDG())==3004 || abs(particle.GetPDG())==3005)
      {
//         if(particle.GetP() < 1.)
//         {
//           fStKFParticleInterface->RemoveParticle(iParticle);
//           continue;
//         }

//         if(particle.GetPhi() > -0.8 && particle.GetPhi() < -0.4)
//         {
//           fStKFParticleInterface->RemoveParticle(iParticle);
//           continue;
//         }
        
        for(int iD=0; iD<particle.NDaughters(); iD++)
        {
          const int daughterId = particle.DaughterIds()[iD];
          const KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterId];
//          if(abs(daughter.GetPDG())==211 && daughter.GetP() > 0.5)
//            fStKFParticleInterface->RemoveParticle(iParticle);
        }
      }
    }
  
	if(isGoodEvent){
        hvtxgood->Fill(bVz);
        hrefmult->Fill(crefmult);
        if(crefmult>60){
	wrefmult->Fill(refmultcor);
	}else{
        wrefmult->Fill(refmultcor, reweight);
	}
	}
  
    int eventId = -1;
    int runId = -1;
    
    if(fFlowAnalysis)
    {
      if(fIsPicoAnalysis) 
      {
        runId   = fPicoDst->event()->runId();
        eventId = fPicoDst->event()->eventId();
      }
      else
      {
        runId   = fMuDst->event()->runId();
        eventId = fMuDst->event()->eventId();
      }
    
      long entryId = GetUniqueEventId(runId, eventId);
      //cout<<"entryId:"<<entryId<<" "<<runId<<" "<<eventId<<endl;
      std::map<long,int>::iterator flowMapIterator = fFlowMap.find(entryId);
      if (flowMapIterator != fFlowMap.end())
      {
        fFlowChain->GetEvent(fFlowMap[GetUniqueEventId(runId, eventId)]);       
        centralityBin = fCentrality;
        //cout<<"fFlowRunId:"<<fFlowRunId<<" "<< fFlowEventId<< " "<< fCentrality<<" "<< endl;
        //cout<<"flow:"<<psi_1_EPD_0 <<" "<< psi_1_EPD_1<<" "<<psi_1_EPD_2<<" "<<psi_1_EPD_3<<endl;
      }
    }
    
    centralityWeight = 1;
    
    fStKFParticlePerformanceInterface->SetMCTracks(mcTracks);
    fStKFParticlePerformanceInterface->SetMCIndexes(mcIndices);    
    fStKFParticlePerformanceInterface->SetCentralityBin(centralityBin);
    fStKFParticlePerformanceInterface->SetCentralityWeight(centralityWeight);
    Int_t nevent = 100000;
    fStKFParticlePerformanceInterface->SetPrintEffFrequency(nevent);

 //   cout<<"eventIdcheck:"<<fPicoDst->event()->eventId()<<endl;
    fStKFParticlePerformanceInterface->PerformanceAnalysis();
//    cout<<"fStKFParticlePerformanceInterface->GetNReconstructedParticles()"<<
//    fStKFParticlePerformanceInterface->GetNReconstructedParticles()<<endl;
//
//
    
    //b
    bool bDEBUG = false;
    if(bDEBUG){/*
    //check track ID:
    for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
    {
       KFParticle particle;
       fStKFParticlePerformanceInterface->GetParticle(particle, iParticle);
       if(fabs(particle.GetPDG())==3312){ 
       cout<<"Cascade found"<<endl;
       cout<<"eventIdcheck:"<<fPicoDst->event()->eventId()<<endl;
       vector<int> cascadedaughterid;
       cascadedaughterid = particle.DaughterIds();
       //cout << "no. of daughers:" << particle.NDaughters() << " " << cascadedaughterid[0] << 
       //" " << cascadedaughterid[1] <<  endl;
       KFParticle daughter1;
       fStKFParticlePerformanceInterface->GetParticle(daughter1, cascadedaughterid[0]);
       KFParticle daughter2;
       fStKFParticlePerformanceInterface->GetParticle(daughter2, cascadedaughterid[1]);
       //cout << "daughter PDG:" << cascadedaughterid.size() << " " <<
       // daughter1.GetPDG() << " " << daughter2.GetPDG() << endl;
       //cout << "daughter no daughters:" << daughter1.DaughterIds().size() << 
       //" " << daughter2.DaughterIds().size() << endl;
       vector<int> lambdadaughterid;
       lambdadaughterid = daughter2.DaughterIds();
       
       KFParticle lambdadaughter1;
       KFParticle lambdadaughter2;

       fStKFParticlePerformanceInterface->GetParticle(lambdadaughter1, lambdadaughterid[0]);
       fStKFParticlePerformanceInterface->GetParticle(lambdadaughter2, lambdadaughterid[1]);

       cout << daughter1.GetPDG() << " " <<lambdadaughter1.GetPDG() << " " << lambdadaughter2.GetPDG() << endl;

       cout << daughter1.GetPx() << " " << daughter1.GetPy() << " " << daughter1.GetPz() << endl;
       cout << lambdadaughter1.GetPx() << " " << lambdadaughter1.GetPy() << " " << lambdadaughter1.GetPz() << endl;
       cout << lambdadaughter2.GetPx() << " " << lambdadaughter2.GetPy() << " " << lambdadaughter2.GetPz() << endl;
       cout << "------" << endl;
       int picopion1id = -999;
       int picopion2id = -999;
       int picoprotonid = -999;
       //Matching: to be improved
       for(unsigned int i=0; i<fPicoDst->numberOfTracks(); i++){
       int picotempid;
       StPicoTrack *t = (StPicoTrack*)fPicoDst->track(i);       
       picotempid = t->id();

       StPicoTrackCovMatrix *covc = fPicoDst->trackCovMatrix(i);
       const StDcaGeometry dcaGc = covc->dcaGeometry();
       Double_t xyzpc[6], CovXyzpc[21];
       dcaGc.GetXYZ(xyzpc,CovXyzpc);
       //cout<<i<<" " << picotempid<<" " <<xyzpc[3]<<" "<<xyzpc[4]<<" "<<xyzpc[5]<<endl;
       if(xyzpc[3]<daughter1.GetPx()+0.00001 && xyzpc[3]>daughter1.GetPx()-0.00001 && xyzpc[4]<daughter1.GetPy()+0.00001 && xyzpc[4]>daughter1.GetPy()-0.00001){ picopion1id = i;} 
//       StPicoTrackCovMatrix *covc = fPicoDst->trackCovMatrix(i);   
//       const StDcaGeometry dcaGc = covc->dcaGeometry();
//       Double_t xyzpc[6], CovXyzpc[21];
//       dcaGc.GetXYZ(xyzpc,CovXyzpc);
//       cout<<xyzpc[3]<<" "<<xyzpc[4]<<" "<<xyzpc[5]<<endl;

//       if(picotempid == lambdadaughterid[0]) { picopion2id = i; }
//       if(picotempid == lambdadaughterid[1]) { picoprotonid = i; }
       if(xyzpc[3]<lambdadaughter1.GetPx()+0.00001 && xyzpc[3]>lambdadaughter1.GetPx()-0.00001 && xyzpc[4]<lambdadaughter1.GetPy()+0.00001 && xyzpc[4]>lambdadaughter1.GetPy()-0.00001){ picopion2id = i;}
       if(xyzpc[3]<lambdadaughter2.GetPx()+0.00001 && xyzpc[3]>lambdadaughter2.GetPx()-0.00001 && xyzpc[4]<lambdadaughter2.GetPy()+0.00001 && xyzpc[4]>lambdadaughter2.GetPy()-0.00001){ picoprotonid = i;}
       }
       cout <<"------" <<endl;
       //cout << cascadedaughterid[0] << " " << lambdadaughterid[0] << " " << lambdadaughterid[1] << endl;
       cout << picopion1id << " " << picopion2id << " " << picoprotonid << endl;

       if(picopion1id>0&&picopion2id>0&&picoprotonid>0){
       StPicoTrack *t_pion1  = (StPicoTrack*)fPicoDst->track(picopion1id);
       StPicoTrack *t_pion2  = (StPicoTrack*)fPicoDst->track(picopion2id);
       StPicoTrack *t_proton = (StPicoTrack*)fPicoDst->track(picoprotonid);

       cout << "nsigma:" << t_pion1->nSigmaPion() << " " << t_pion1->nSigmaProton() <<endl;      
       cout << "nsigma:" << t_pion2->nSigmaPion() << " " << t_pion2->nSigmaProton() <<endl;
       cout << "nsigma:" << t_proton->nSigmaPion() << " " << t_proton->nSigmaProton() <<endl;       
       }else{cout<<"particle not found"<<endl;}
       //cout << "lambda d:" << lambdadaughterid.size() << endl;
       
         
//       double pionmass    = 0.139570;
//       double protonmass  = 0.938;
       
//       cout << "cascade pT:" << particle.GetPt() << endl;
//       cout << "sumpico pT:" << t_pion1->gPt() + t_pion2->gPt() + t_proton->gPt()<<endl;
//       cout << "pico pT:" << t_pion1->pMom().X() << " " << t_pion1->pMom().Y() << " " << t_pion1->pMom().Z()<<endl;
//       cout << "pico pT:" << t_pion2->pMom().X() << " " << t_pion2->pMom().Y() << " " << t_pion2->pMom().Z()<<endl;
//       cout << "pico pT:" << t_proton->pMom().X() << " " << t_proton->pMom().Y() << " " << t_proton->pMom().Z()<<endl;
//
       cout << "cascade mass:" << particle.GetMass() << endl;
       cout << endl;
       }
//GetListOfDaughterTracks
    }
*/
  }




    if(bWriteTree){

/* 
       if(fIsPicoAnalysis){
	beventid = fPicoDst->event()->eventId();
	brunid = fPicoDst->event()->runId();

  const TVector3 picoPV = fPicoDst->event()->primaryVertex();
  const TVector3 picoPVError = fPicoDst->event()->primaryVertexError();

        bVr = sqrt(picoPV.x()*picoPV.x()+picoPV.y()*picoPV.y());
        bVz = picoPV.z(); 
        bVrerr = sqrt(picoPVError.x()*picoPVError.x()+picoPVError.y()*picoPVError.y());
        bVzerr = picoPVError.z();

	}
       else{
	beventid = fMuDst->event()->eventId();
	brunid = fMuDst->event()->runId();

  float bestRank=-1000000;
  int bestPV=-1;
  for(unsigned int iPV=0; iPV<fMuDst->numberOfPrimaryVertices(); iPV++)
  {
    StMuPrimaryVertex *Vtx = fMuDst->primaryVertex(iPV);
    if(!Vtx) continue;
    if (bestRank < Vtx->ranking()) {
      bestRank = Vtx->ranking();
      bestPV = iPV;
    }
    else continue;

    bVr = sqrt(Vtx->position().x()*Vtx->position().x()+Vtx->position().y()*Vtx->position().y());
    bVz = Vtx->position().z();
    bVrerr = sqrt(Vtx->posError().x()*Vtx->posError().x()+Vtx->posError().y()*Vtx->posError().y());
    bVzerr = Vtx->posError().z();
	}

}

	hvtx->Fill(bVz);
if(fabs(bVz)<70 && bVr<2){
        hvtxgood->Fill(bVz);
}
*/
    //check track ID:
    for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
    {

       //KFParticle particle;
       //fStKFParticlePerformanceInterface->GetParticle(particle, iParticle);
        KFParticle particle;
//       const KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];
             bool isMCParticle = fStKFParticlePerformanceInterface->GetParticle(particle, iParticle);

       //int beventid;
       //int brunid;
       //int bparticleid;
       //float bparticlemass;
        bismc =0;
        if(isMCParticle) bismc=1;

        bparticleid   = particle.GetPDG();
        bparticlemass = particle.GetMass();

	if(fIsPicoAnalysis)
      {
        brunid   = fPicoDst->event()->runId();
        brefmult = fPicoDst->event()->refMult();
        btofmult = fPicoDst->event()->btofTrayMultiplicity();
      }
      else
      {
        brunid   = fMuDst->event()->runId();
        brefmult = fMuDst->event()->refMult();
        btofmult = fMuDst->event()->btofTrayMultiplicity();
      }

        bx = particle.GetX();
        by = particle.GetY();
        bz = particle.GetZ();
        bpx = particle.GetPx();
        bpy = particle.GetPy();
        bpz = particle.GetPz();
        bdl = particle.GetDecayLength();
         //if(particle.IdTruth()!=0){ 
        //cout<<"IdTruth():"<<particle.IdTruth()<<" "<< particle.IdParentMcVx()  <<endl;
//         }

/*
       if(fabs(bparticleid)==3312 || fabs(bparticleid)==3334){ 
       vector<int> cascadedaughterid;
       cascadedaughterid = particle.DaughterIds();
       KFParticle daughter1;
       fStKFParticlePerformanceInterface->GetParticle(daughter1, cascadedaughterid[0]);
       KFParticle daughter2;
       fStKFParticlePerformanceInterface->GetParticle(daughter2, cascadedaughterid[1]);


       vector<int> lambdadaughterid;
       lambdadaughterid = daughter2.DaughterIds();
        
       KFParticle lambdadaughter1;
       KFParticle lambdadaughter2;

       fStKFParticlePerformanceInterface->GetParticle(lambdadaughter1, lambdadaughterid[0]);
       fStKFParticlePerformanceInterface->GetParticle(lambdadaughter2, lambdadaughterid[1]);

       if(fabs(lambdadaughter2.GetPDG())==221 || fabs(lambdadaughter2.GetPDG())==321){//daughter 1 is pion
              fStKFParticlePerformanceInterface->GetParticle(lambdadaughter2, lambdadaughterid[0]);
              fStKFParticlePerformanceInterface->GetParticle(lambdadaughter1, lambdadaughterid[1]);//swap
       }

        bbachid = daughter1.GetPDG();       
        bbachpx = daughter1.GetPx();
        bbachpy = daughter1.GetPy();
        bbachpz = daughter1.GetPz();
        bbachmass = daughter1.GetMass();

	bpionid = lambdadaughter1.GetPDG();
        bpionpx = lambdadaughter1.GetPx();
        bpionpy = lambdadaughter1.GetPy();
        bpionpz = lambdadaughter1.GetPz();
        bpionmass = lambdadaughter1.GetMass();

        bprotonid = lambdadaughter2.GetPDG();
        bprotonpx = lambdadaughter2.GetPx();
	bprotonpy = lambdadaughter2.GetPy();
        bprotonpz = lambdadaughter2.GetPz();
        bprotonmass = lambdadaughter2.GetMass();
       
   }//cascade loop
  */ 
       //cout << daughter1.GetPx() << " " << daughter1.GetPy() << " " << daughter1.GetPz() << endl;
       //cout << lambdadaughter1.GetPx() << " " << lambdadaughter1.GetPy() << " " << lambdadaughter1.GetPz() << endl;
       //cout << lambdadaughter2.GetPx() << " " << lambdadaughter2.GetPy() << " " << lambdadaughter2.GetPz() << endl;
//       cout << "------" << endl;
//       int picopion1id = -999;
//       int picopion2id = -999;
//       int picoprotonid = -999;


       //Matching: to be improved
/*       
       for(unsigned int i=0; i<fPicoDst->numberOfTracks(); i++){
       int picotempid;
       StPicoTrack *t = (StPicoTrack*)fPicoDst->track(i);       
       picotempid = t->id();

       StPicoTrackCovMatrix *covc = fPicoDst->trackCovMatrix(i);
       const StDcaGeometry dcaGc = covc->dcaGeometry();
       Double_t xyzpc[6], CovXyzpc[21];
       dcaGc.GetXYZ(xyzpc,CovXyzpc);
       //cout<<i<<" " << picotempid<<" " <<xyzpc[3]<<" "<<xyzpc[4]<<" "<<xyzpc[5]<<endl;
       if(xyzpc[3]<daughter1.GetPx()+0.00001 && xyzpc[3]>daughter1.GetPx()-0.00001 && xyzpc[4]<daughter1.GetPy()+0.00001 && xyzpc[4]>daughter1.GetPy()-0.00001){ picopion1id = i;} 
//       StPicoTrackCovMatrix *covc = fPicoDst->trackCovMatrix(i);   
//       const StDcaGeometry dcaGc = covc->dcaGeometry();
//       Double_t xyzpc[6], CovXyzpc[21];
//       dcaGc.GetXYZ(xyzpc,CovXyzpc);
//       cout<<xyzpc[3]<<" "<<xyzpc[4]<<" "<<xyzpc[5]<<endl;

//       if(picotempid == lambdadaughterid[0]) { picopion2id = i; }
//       if(picotempid == lambdadaughterid[1]) { picoprotonid = i; }
       if(xyzpc[3]<lambdadaughter1.GetPx()+0.00001 && xyzpc[3]>lambdadaughter1.GetPx()-0.00001 && xyzpc[4]<lambdadaughter1.GetPy()+0.00001 && xyzpc[4]>lambdadaughter1.GetPy()-0.00001){ picopion2id = i;}
       if(xyzpc[3]<lambdadaughter2.GetPx()+0.00001 && xyzpc[3]>lambdadaughter2.GetPx()-0.00001 && xyzpc[4]<lambdadaughter2.GetPy()+0.00001 && xyzpc[4]>lambdadaughter2.GetPy()-0.00001){ picoprotonid = i;}
       }
       cout <<"------" <<endl;
       cout << cascadedaughterid[0] << " " << lambdadaughterid[0] << " " << lambdadaughterid[1] << endl;
       cout << picopion1id << " " << picopion2id << " " << picoprotonid << endl;
       cout << daughter1.GetId() << " " << lambdadaughter1.GetId() << " "<< lambdadaughter2.GetId() << endl;

       if(picopion1id>0&&picopion2id>0&&picoprotonid>0){
       StPicoTrack *t_pion1  = (StPicoTrack*)fPicoDst->track(picopion1id);
       StPicoTrack *t_pion2  = (StPicoTrack*)fPicoDst->track(picopion2id);
       StPicoTrack *t_proton = (StPicoTrack*)fPicoDst->track(picoprotonid);

       cout << "nsigma:" << t_pion1->nSigmaPion() << " " << t_pion1->nSigmaProton() <<endl;      
       cout << "nsigma:" << t_pion2->nSigmaPion() << " " << t_pion2->nSigmaProton() <<endl;
       cout << "nsigma:" << t_proton->nSigmaPion() << " " << t_proton->nSigmaProton() <<endl;       
       }else{cout<<"particle not found"<<endl;}
       //cout << "lambda d:" << lambdadaughterid.size() << endl;
       
         
//       double pionmass    = 0.139570;
//       double protonmass  = 0.938;
       
//       cout << "cascade pT:" << particle.GetPt() << endl;
//       cout << "sumpico pT:" << t_pion1->gPt() + t_pion2->gPt() + t_proton->gPt()<<endl;
//       cout << "pico pT:" << t_pion1->pMom().X() << " " << t_pion1->pMom().Y() << " " << t_pion1->pMom().Z()<<endl;
//       cout << "pico pT:" << t_pion2->pMom().X() << " " << t_pion2->pMom().Y() << " " << t_pion2->pMom().Z()<<endl;
//       cout << "pico pT:" << t_proton->pMom().X() << " " << t_proton->pMom().Y() << " " << t_proton->pMom().Z()<<endl;
//
       cout << "cascade mass:" << particle.GetMass() << endl;

       }//GetListOfDaughterTracks
*/

	//bool _fill_lambda_tree;
	//_fill_lambda_tree = false;
	//_fill_lambda_tree = true;
	if(_fill_lambda_tree && isGoodEvent && bparticlemass<1.2){
		if(fabs(particle.GetPDG())==3122){

		KFParticleSIMD tempSIMDParticle(particle);
	        float_v l,dl;
                KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
                tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
                ld_ldl = l[0]/dl[0];
                ld_l = l[0];
                ld_dl = dl[0];

                bpl = ld_l/(sqrt(bpx*bpx+bpy*bpy+bpz*bpz)/bparticlemass);

		//cout<<"GetDistanceFromVertex:"<<tempSIMDParticle.GetDistanceFromVertex(pv)<<endl;
                //cout<<"GetLifeTime:"<<tempSIMDParticle.GetLifeTime()<<endl;
                ld_bdfvtx = tempSIMDParticle.GetDistanceFromVertex(pv)[0];
		ld_bdfvtx_xy = tempSIMDParticle.GetDistanceFromVertexXY(pv)[0];
		ld_bdfvtxdev_xy = tempSIMDParticle.GetDeviationFromVertexXY(pv)[0];

		//cout<<"ld_chi2topo:"<<tempSIMDParticle.Chi2()[0]<<"/"<<tempSIMDParticle.NDF()[0]<<endl;
                tempSIMDParticle.SetProductionVertex(pv);
		//cout<<"pv:"<<pv.GetX()<<" "<<pv.GetErrX()<<" "<<pv.GetY()<<" "<<pv.GetErrY()<<pv.GetZ()<<" "<<pv.GetErrZ()<<   endl;
		//cout<<"ld_chi2topo:"<<tempSIMDParticle.GetChi2()[0]<<"/"<<tempSIMDParticle.NDF()[0]<<endl;
		//cout<<"ld_dca:"<<tempSIMDParticle.GetDistanceFromVertex(pv)<<" "<<tempSIMDParticle.GetDeviationFromVertex(pv)<<endl;

		
                //cout<<"GetLifeTime:"<<tempSIMDParticle.GetLifeTime()<<endl;
                //cout<<"GetDistanceFromVertex2:"<<tempSIMDParticle.GetDistanceFromVertex(pv)<<endl;
	        ld_bdfvtx2 = tempSIMDParticle.GetDistanceFromVertex(pv)[0];
		ld_lifetime = tempSIMDParticle.GetLifeTime()[0];
		

                ld_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
                ld_chi2ndf = particle.Chi2()/particle.NDF();

			for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
				{
                                int order[4] = {0, 1, 2, 3};
                                const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
			        KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];
  					if(iDaughter==0){
//						cout<<"pdg:"<<daughter.GetPDG()<<endl;
						chi2primary_pi = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
                                                nhits_ld_pi = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
                                                dca_pi = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
					}
				        if(iDaughter==1){
//						cout<<"pdg:"<<daughter.GetPDG()<<endl;
						chi2primary_proton = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
                                                nhits_ld_proton = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
                                                dca_proton = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
					}

				}
  
		 bmcpx=-999;
                 bmcpy=-999;
                 bmcpz=-999;
                 bmcpl=-999;

			//if(isMCParticle && fProcessSignal){
			if(isMCParticle){
				int iMCPart = fStKFParticlePerformanceInterface->GetParticleMCId(iParticle);
				StMuMcTrack *mcTrack = fMuDst->MCtrack(iMCPart);
				bmcpx = mcTrack->Pxyz().x();
				bmcpy = mcTrack->Pxyz().y();
				bmcpz = mcTrack->Pxyz().z();
				bmcidvx = mcTrack->IdVx();
				bmcidvxend = mcTrack->IdVxEnd();
				StMuMcVertex *mcvx = fMuDst->MCvertex(bmcidvx-1);
			        StMuMcVertex *mcvxend = fMuDst->MCvertex(bmcidvxend-1);

				//cout<<"particle:"<<particle.X()<<" "<<particle.Y()<<" "<<particle.Z()<<" "<<endl;
	                        //cout<<"mcpvend:"<<mcvxend->XyzV().x()<<" "<<mcvxend->XyzV().y()<<" "<<mcvxend->XyzV().z()<<endl;
				//cout<<endl;
				//cout<<"particle p:"<<particle.Px()<<" "<<particle.Py()<<" "<<particle.Pz()<<" "<<endl;
			        //cout<<"mcpvend:"<<bmcpx<<" "<<bmcpy<<" "<<bmcpz<<endl;
				//cout<<endl;
                                //cout<<"pv:"<<pv.X()<<" "<<pv.Y()<<" "<<pv.Z()<<endl;
	                        //cout<<"mcpv:"<<mcvx->XyzV().x()<<" "<<mcvx->XyzV().y()<<" "<<mcvx->XyzV().z()<<endl;
				
				bmcx = mcvxend->XyzV().x();
				bmcy = mcvxend->XyzV().y();
				bmcz = mcvxend->XyzV().z(); 
				bx = particle.X();
				by = particle.Y();
				bz = particle.Z();
				
				bmcl = sqrt( (mcvx->XyzV().x()-mcvxend->XyzV().x())*(mcvx->XyzV().x()-mcvxend->XyzV().x()) + (mcvx->XyzV().y()-mcvxend->XyzV().y())*(mcvx->XyzV().y()-mcvxend->XyzV().y()) + (mcvx->XyzV().z()-mcvxend->XyzV().z())*(mcvx->XyzV().z()-mcvxend->XyzV().z()) ) ;
				bmcpl = bmcl * 1.115683 / sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz) ;
			}
		 lambda_tree->Fill();
	   	 }
         }



	if(isGoodEvent){
	if(particle.GetPDG()==310){
  	KFParticleSIMD tempSIMDParticle(particle);
  	float_v l,dl;
  	KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  	tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
  	ld_ldl = l[0]/dl[0];
  	ld_l = l[0];
  	ld_dl = dl[0];
  	tempSIMDParticle.SetProductionVertex(pv);
  	ld_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
  	ld_chi2ndf = particle.Chi2()/particle.NDF();
		for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
        	{
      		int order[4] = {0, 1, 2, 3};
  		const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
  		KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];  
  			if(iDaughter==0){
			chi2primary_proton = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
        		dca_proton = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
			}
  			if(iDaughter==1){
			chi2primary_pi = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
			dca_pi = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
			}
  	      	}
		bmcpx=-999;
		bmcpy=-999;
		bmcpz=-999;
		if(isMCParticle && fProcessSignal){
		int iMCPart = fStKFParticlePerformanceInterface->GetParticleMCId(iParticle);
		StMuMcTrack *mcTrack = fMuDst->MCtrack(iMCPart);
		bmcpx = mcTrack->Pxyz().x();
		bmcpy = mcTrack->Pxyz().y();
		bmcpz = mcTrack->Pxyz().z();
		}
      		ks_tree->Fill();
	}
	}


	if(fabs(particle.GetPDG())==3312 && isGoodEvent){

  	KFParticleSIMD tempSIMDParticle(particle);
  	float_v l,dl;
  	KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  	tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
  	xi_ldl = l[0]/dl[0];
  	xi_l = l[0];
  	xi_dl = dl[0];

  	tempSIMDParticle.SetProductionVertex(pv);
  	xi_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
  	xi_chi2ndf = particle.Chi2()/particle.NDF();

		for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
        	{
  		int order[4] = {0, 1, 2, 3};
  		const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
  		KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];

			if(iDaughter==0){
        		chi2primary_xi_bach = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
        }
			if(iDaughter==1){
        		chi2primary_xi_ld = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
		  	KFParticleSIMD ttempSIMDParticle(daughter);
  			float_v tl,tdl;

  			KFParticleSIMD tpv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  			ttempSIMDParticle.GetDistanceToVertexLine(tpv, tl, tdl);
  			xi_ld_ldl = tl[0]/tdl[0];
  			xi_ld_l = tl[0];
 	 		ttempSIMDParticle.SetProductionVertex(tpv);
  			xi_ld_chi2topo = double(ttempSIMDParticle.Chi2()[0])/double(ttempSIMDParticle.NDF()[0]);
  			xi_ld_chi2ndf = daughter.Chi2()/daughter.NDF();

        			for(int jDaughter=0; jDaughter<daughter.NDaughters(); jDaughter++){
                		int jorder[4] = {0, 1, 2, 3};
                		const int jdaughterParticleIndex = daughter.DaughterIds()[jorder[jDaughter]];
                  		KFParticle granddaughter = fStKFParticleInterface->GetParticles()[jdaughterParticleIndex];
                			if(jDaughter==0){
                        		chi2primary_xi_proton = granddaughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
                                	}
                			if(jDaughter==1){
                        		chi2primary_xi_pi = granddaughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
                                	}
                		}
        		}
		}
		bmcpx=-999;
		bmcpy=-999;
		bmcpz=-999;
      		if(isMCParticle && fProcessSignal){
        	int iMCPart = fStKFParticlePerformanceInterface->GetParticleMCId(iParticle);
        	StMuMcTrack *mcTrack = fMuDst->MCtrack(iMCPart);
         	bmcpx = mcTrack->Pxyz().x();
         	bmcpy = mcTrack->Pxyz().y();
         	bmcpz = mcTrack->Pxyz().z();
       		}
      		cascade_tree->Fill();
	}


	if(fabs(particle.GetPDG())==3334 && isGoodEvent){

  	KFParticleSIMD tempSIMDParticle(particle);
  	float_v l,dl;
  	KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  	tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
  	om_ldl = l[0]/dl[0];
  	om_l = l[0];
  	om_dl = dl[0];


  	tempSIMDParticle.SetProductionVertex(pv);
  	om_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
  	om_chi2ndf = particle.Chi2()/particle.NDF();

	for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
	{
  	int order[4] = {0, 1, 2, 3};
  	const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
  	KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];

		if(iDaughter==0){
		chi2primary_om_bach = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

        	dedx_om_bach = fStKFParticleInterface->GetdEdX(daughter.DaughterIds()[0]);
        	nhits_om_bach = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);

		bbachpx = daughter.GetPx();
		bbachpy = daughter.GetPy();
		bbachpz = daughter.GetPz();

		}
  		if(iDaughter==1){
		chi2primary_om_ld = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  		KFParticleSIMD ttempSIMDParticle(daughter);
  		float_v tl,tdl;

  		KFParticleSIMD tpv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  		ttempSIMDParticle.GetDistanceToVertexLine(tpv, tl, tdl);
  		om_ld_ldl = tl[0]/tdl[0];
  		om_ld_l = tl[0];
  		ttempSIMDParticle.SetProductionVertex(tpv);
  		om_ld_chi2topo = double(ttempSIMDParticle.Chi2()[0])/double(ttempSIMDParticle.NDF()[0]);
  		om_ld_chi2ndf = daughter.Chi2()/daughter.NDF();

       			for(int jDaughter=0; jDaughter<daughter.NDaughters(); jDaughter++){
			int jorder[4] = {0, 1, 2, 3};
			const int jdaughterParticleIndex = daughter.DaughterIds()[jorder[jDaughter]];
		  	KFParticle granddaughter = fStKFParticleInterface->GetParticles()[jdaughterParticleIndex];
				if(jDaughter==0){
        			chi2primary_om_proton = granddaughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

        			dedx_om_proton = fStKFParticleInterface->GetdEdX(granddaughter.DaughterIds()[0]);
        			nhits_om_proton = fStKFParticleInterface->GetNHits(granddaughter.DaughterIds()[0]);

       				bprotonpx = granddaughter.GetPx();
        			bprotonpy = granddaughter.GetPy();
        			bprotonpz = granddaughter.GetPz();

		        	}
                		if(jDaughter==1){
                        	chi2primary_om_pi = granddaughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
				dedx_om_pi = fStKFParticleInterface->GetdEdX(granddaughter.DaughterIds()[0]);
		        	nhits_om_pi = fStKFParticleInterface->GetNHits(granddaughter.DaughterIds()[0]);

        			bpionpx = granddaughter.GetPx();
        			bpionpy = granddaughter.GetPy();
        			bpionpz = granddaughter.GetPz();
			
		                }
		        }	
	  	}
	}
		bmcpx=-999;
		bmcpy=-999;
		bmcpz=-999;
   	if(isMCParticle && fProcessSignal){
        	int iMCPart = fStKFParticlePerformanceInterface->GetParticleMCId(iParticle);
		StMuMcTrack *mcTrack = fMuDst->MCtrack(iMCPart);
         	bmcpx = mcTrack->Pxyz().x();
         	bmcpy = mcTrack->Pxyz().y();
         	bmcpz = mcTrack->Pxyz().z();
                 }
      		omega_tree->Fill();
	}

//quick sanity check
/*
if(fabs(particle.GetPDG())==211 && isGoodEvent){
int partId  = particle.Id();//get the particle id
int trackId = particle.DaughterIds()[0];//get the track id
//if(fStKFParticleInterface->GetdEdXNSigmaPion(trackId)>3){
if(1){
cout<<"nsigmapion:"<<partId<<" "<<trackId<<" "<<fStKFParticleInterface->GetdEdXNSigmaPion(trackId)<<" "<<fStKFParticleInterface->GetdEdXNSigmaKaon(trackId)<<" "<<fStKFParticleInterface->GetdEdXNSigmaProton(trackId) <<endl;
cout<<"dEdX:"<<partId<<" "<<trackId<<" "<<fStKFParticleInterface->GetdEdX(trackId)<<" "<<fStKFParticleInterface->GetNHits(trackId)<<endl;
}
}
*/

//he3 tree
	if(fabs(particle.GetPDG())==1000020030 && isGoodEvent && _fill_he3_tree){
	int trackId = particle.DaughterIds()[0];//get the track id
	bdedx = fStKFParticleInterface->GetdEdX(trackId);
	bnhits = fStKFParticleInterface->GetNHits(trackId);
	bdca = fStKFParticleInterface->Getdca(trackId);
	chi2primary_he = particle.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

	if(isMCParticle){
		int iMCPart = fStKFParticlePerformanceInterface->GetParticleMCId(iParticle);
		StMuMcTrack *mcTrack = fMuDst->MCtrack(iMCPart);
		bmcpx = mcTrack->Pxyz().x();
		bmcpy = mcTrack->Pxyz().y();
		bmcpz = mcTrack->Pxyz().z();
	}else{
		bmcpx=-999;
		bmcpy=-999;
		bmcpz=-999;
	}

		he3_tree->Fill();
	}

	if(  ( fabs(particle.GetPDG())==3004 || fabs(particle.GetPDG())==3005  )  && isGoodEvent){

	KFParticleSIMD tempSIMDParticle(particle);
	float_v l,dl;
 	KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  	tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
  	ht_ldl = l[0]/dl[0];
  	ht_l = l[0];
  	ht_dl = dl[0];
        bpl = ht_l/(sqrt(bpx*bpx+bpy*bpy+bpz*bpz)/bparticlemass);
	ht_bdfvtx = tempSIMDParticle.GetDistanceFromVertex(pv)[0];

  	tempSIMDParticle.SetProductionVertex(pv);
        ht_bdfvtx2 = tempSIMDParticle.GetDistanceFromVertex(pv)[0];
        ht_lifetime = tempSIMDParticle.GetLifeTime()[0];
  	ht_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
  	ht_chi2 = double(tempSIMDParticle.Chi2()[0]);
  	ht_NDF = double(tempSIMDParticle.NDF()[0]);
  	ht_chi2ndf = particle.Chi2()/particle.NDF();

	for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
        {
      	int order[4] = {0, 1, 2, 3};
  	const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
  	KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];
            if(iDaughter==0)
		{
			chi2primary_pi = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
                        nhits_pi = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
                        dca_pi = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
                        dedx_pi = fStKFParticleInterface->GetdEdX(daughter.DaughterIds()[0]);
			px_pi = daughter.GetPx();
			py_pi = daughter.GetPy();
			pz_pi = daughter.GetPz();
		}
            if(iDaughter==1)
		{
//			cout<<"daughter:"<<sqrt(daughter.GetPx()*daughter.GetPx()+daughter.GetPy()*daughter.GetPy())<<endl;
			chi2primary_he = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
                        nhits_he = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
                        dca_he = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
                        dedx_he = fStKFParticleInterface->GetdEdX(daughter.DaughterIds()[0]);
			px_he = daughter.GetPx();
                        py_he = daughter.GetPy();
                        pz_he = daughter.GetPz();
		}
        }
	           bmcpx=-999;
                   bmcpy=-999;
                   bmcpz=-999;
		   bmcl=-999;
                   bmcpl=-999;

                          
                          if(isMCParticle){
                                  int iMCPart = fStKFParticlePerformanceInterface->GetParticleMCId(iParticle);
                                  StMuMcTrack *mcTrack = fMuDst->MCtrack(iMCPart);
                                  bmcpx = mcTrack->Pxyz().x();
                                  bmcpy = mcTrack->Pxyz().y();
                                  bmcpz = mcTrack->Pxyz().z();
                                  bmcidvx = mcTrack->IdVx();
                                  bmcidvxend = mcTrack->IdVxEnd();
                                  StMuMcVertex *mcvx = fMuDst->MCvertex(bmcidvx-1);
                                  StMuMcVertex *mcvxend = fMuDst->MCvertex(bmcidvxend-1);
//				cout<<"pv:"<<pv.X()<<" "<<pv.Y()<<" "<<pv.Z()<<endl;
//				cout<<"mcpv:"<<mcvx->XyzV().x()<<" "<<mcvx->XyzV().y()<<" "<<mcvx->XyzV().z()<<endl;
	
                                  bmcl = sqrt( (mcvx->XyzV().x()-mcvxend->XyzV().x())*(mcvx->XyzV().x()-mcvxend->XyzV().x()) + (mcvx->XyzV().y()-mcvxend->XyzV().y())*(mcvx->XyzV().y()-mcvxend->XyzV().y()) + (mcvx->XyzV().z()-mcvxend->XyzV().z())*(mcvx->XyzV().z()-mcvxend->XyzV().z()) ) ;
				if(fabs(particle.GetPDG())==3004){
                                  bmcpl = bmcl * 2.99131 / sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz) ;
				}
				if(fabs(particle.GetPDG())==3005){
                                  bmcpl = bmcl * 3.9239 / sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz) ;
                                }
                          }
      		htriton_tree->Fill();
	}
/*
if(fabs(particle.GetPDG())==3005 && isGoodEvent){

       KFParticleSIMD tempSIMDParticle(particle);
         float_v l,dl;
         KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
         tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
         hl_ldl = l[0]/dl[0];
         hl_l = l[0];
         hl_dl = dl[0];

         tempSIMDParticle.SetProductionVertex(pv);
         hl_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
         hl_chi2 = double(tempSIMDParticle.Chi2()[0]);
         hl_NDF = double(tempSIMDParticle.NDF()[0]);
         hl_chi2ndf = particle.Chi2()/particle.NDF();

       for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
               {
             int order[4] = {0, 1, 2, 3};
         const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
         KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];
                   if(iDaughter==0){chi2primary_h4 = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());}
                   if(iDaughter==1){chi2primary_pi = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());}

                               }

             h4lambda_tree->Fill();
       }
*/
	if((fabs(particle.GetPDG())==3009  || fabs(particle.GetPDG())==3007) && isGoodEvent && bparticlemass<5.5){

           KFParticleSIMD tempSIMDParticle(particle);
           float_v l,dl;
           KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
           tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
           hl_ldl = l[0]/dl[0];

           tempSIMDParticle.SetProductionVertex(pv);
           hl_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
           hl_chi2ndf = particle.Chi2()/particle.NDF();

           h_tree->Fill();
        	 }


    	}//loop through paritcles

    }//write tree bool

//cout<<"fStoreTmvaNTuples:"<<fStoreTmvaNTuples<<" "<< fStKFParticlePerformanceInterface->GetNReconstructedParticles()<<endl;



    if(fStoremctree && isGoodEvent){

    Int_t NoMuMcTracks = fMuDst->numberOfMcTracks();
         for (Int_t k = 0; k < NoMuMcTracks; k++) {

	      StMuMcTrack *mcTrack = fMuDst->MCtrack(k);

        bmcparticleid = mcTrack->GePid();
	if(bmcparticleid==40002){
	bmcparticleid = 3334;
	}
        if(bmcparticleid==99999){
        bmcparticleid = 3312;
        }
	if(bmcparticleid==10018 || bmcparticleid==18){
        bmcparticleid = 3122;
	}
        if(bmcparticleid==49){
        bmcparticleid = 1000020030;
        }
	if(bmcparticleid==61053){	
	bmcparticleid = 3004;
	}
        if(bmcparticleid==61055){
        bmcparticleid = 3005;
        }
        if(bmcparticleid==707){
        bmcparticleid = 310;
        }

	 bmcrawpx = mcTrack->Pxyz().x();
         bmcrawpy = mcTrack->Pxyz().y();
         bmcrawpz = mcTrack->Pxyz().z();		
         bmcrefmult = brefmult;

	if(abs(bmcparticleid)==3334){
	omega_mc_tree->Fill();
	}
        if(abs(bmcparticleid)==3312){
        cascade_mc_tree->Fill();
        }
	if(abs(bmcparticleid)==3122){
	bmcidvx = mcTrack->IdVx();
        bmcidvxend = mcTrack->IdVxEnd();
	if(bmcidvx!=1) continue;
	if(bmcidvxend==0) continue;
        StMuMcVertex *mcvx = fMuDst->MCvertex(bmcidvx-1);
        StMuMcVertex *mcvxend = fMuDst->MCvertex(bmcidvxend-1);
        bmcrawl = sqrt( (mcvx->XyzV().x()-mcvxend->XyzV().x())*(mcvx->XyzV().x()-mcvxend->XyzV().x()) + (mcvx->XyzV().y()-mcvxend->XyzV().y())*(mcvx->XyzV().y()-mcvxend->XyzV().y()) + (mcvx->XyzV().z()-mcvxend->XyzV().z())*(mcvx->XyzV().z()-mcvxend->XyzV().z()) );
	bmcrawpl = bmcrawl* 1.115683 / sqrt(bmcrawpx*bmcrawpx+bmcrawpy*bmcrawpy+bmcrawpz*bmcrawpz) ;
	lambda_mc_tree->Fill();
	}
        if(bmcparticleid==310){
        bmcidvx = mcTrack->IdVx();
        bmcidvxend = mcTrack->IdVxEnd();
        if(bmcidvx!=1) continue;
        if(bmcidvxend==0) continue;
        StMuMcVertex *mcvx = fMuDst->MCvertex(bmcidvx-1);
        StMuMcVertex *mcvxend = fMuDst->MCvertex(bmcidvxend-1);
        bmcrawl = sqrt( (mcvx->XyzV().x()-mcvxend->XyzV().x())*(mcvx->XyzV().x()-mcvxend->XyzV().x()) + (mcvx->XyzV().y()-mcvxend->XyzV().y())*(mcvx->XyzV().y()-mcvxend->XyzV().y()) + (mcvx->XyzV().z()-mcvxend->XyzV().z())*(mcvx->XyzV().z()-mcvxend->XyzV().z()) );
        bmcrawpl = bmcrawl* 0.497611 / sqrt(bmcrawpx*bmcrawpx+bmcrawpy*bmcrawpy+bmcrawpz*bmcrawpz) ;
        ks_mc_tree->Fill();
        }
        if(abs(bmcparticleid)==1000020030){
        he3_mc_tree->Fill();
        }
        if(abs(bmcparticleid)==3004 || abs(bmcparticleid)==3005){
        bmcidvx = mcTrack->IdVx();
        bmcidvxend = mcTrack->IdVxEnd();
        if(bmcidvx!=1) continue;
        if(bmcidvxend==0) continue;
        StMuMcVertex *mcvx = fMuDst->MCvertex(bmcidvx-1);
        StMuMcVertex *mcvxend = fMuDst->MCvertex(bmcidvxend-1);
        bmcrawl = sqrt( (mcvx->XyzV().x()-mcvxend->XyzV().x())*(mcvx->XyzV().x()-mcvxend->XyzV().x()) + (mcvx->XyzV().y()-mcvxend->XyzV().y())*(mcvx->XyzV().y()-mcvxend->XyzV().y()) + (mcvx->XyzV().z()-mcvxend->XyzV().z())*(mcvx->XyzV().z()-mcvxend->XyzV().z()) );
        if(abs(bmcparticleid)==3004){
        bmcrawpl = bmcrawl* 2.99131/ sqrt(bmcrawpx*bmcrawpx+bmcrawpy*bmcrawpy+bmcrawpz*bmcrawpz) ;
        }
        if(abs(bmcparticleid)==3005){
        bmcrawpl = bmcrawl* 3.924/ sqrt(bmcrawpx*bmcrawpx+bmcrawpy*bmcrawpy+bmcrawpz*bmcrawpz) ;
        }
 
	htriton_mc_tree->Fill();
	}


	}//mctrack loop

	}



    if(fStoreTmvaNTuples)
    {
      for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
      {
        KFParticle particle;
        bool isMCParticle = fStKFParticlePerformanceInterface->GetParticle(particle, iParticle);

 //       cout<<"isMCParticle:"<<isMCParticle<< " "<< particle.GetPDG()<<endl;      

//TOBEREVERTED //REVERTED
        if( !( (fProcessSignal && isMCParticle) || (!fProcessSignal && !isMCParticle) ) ) continue;
//if(isMCParticle) continue;          
       
        for(int iNTuple=0; iNTuple<fNNTuples; iNTuple++)
        {
          if( particle.GetPDG() == fNTuplePDG[iNTuple] )
          {
            GetParticleParameters(iNTuple, particle);
            double side_mass = particle.GetMass();

            //bool sideband = (side_mass>1.06 && side_mass<1.10) || (side_mass>1.13 && side_mass<1.17);

            bool sideband;
	if(iNTuple==0){ //lambda
            sideband = (side_mass>1.085 && side_mass<1.10) || (side_mass>1.125 && side_mass<1.14);//for consistency with long code
sideband = true;//for testing
	}
        if(iNTuple==1){ //omega

//	if(isMCParticle){

//	  int iMCPart = fStKFParticlePerformanceInterface->GetParticleMCId(iParticle);

//cout<<"mcPartid:"<<iMCPart<<endl;
//KFMCParticle &mcPart = vMCParticles[iMCPart];
//cout<<"mcPart:"<<mcPart.GetPDG()<<endl;
// mcTracks 
//StMuMcTrack *mcTrack = fMuDst->MCtrack(iMCPart);
//double mcpt = mcTrack->Pxyz().x();
//double mcpt = mcTrack->pT;
//cout<<"mcPt:"<<mcpt<<" "<< particle.GetPx()<<endl;

//	}

            //sideband = (side_mass>1.61 && side_mass<1.66) || (side_mass>1.685 && side_mass<1.85);//for testing
            sideband = (side_mass>1.61 && side_mass<1.66) || (side_mass>1.685 && side_mass<1.73);//for consistency with long code

//TOBEREVERTED
	//  sideband = true;//for testing
  	    sideband = side_mass<1.78;
//            cout<<"fProcessSignal:"<<fProcessSignal<<" isMCParticle:"<<isMCParticle<<endl;
	}
	if(iNTuple==2){
	    sideband = (side_mass>1.25 && side_mass<1.30) || (side_mass>1.35 && side_mass<1.40);//for consistency with long code
            //TOBEREVERTED
            sideband = true;
	}
            if( fProcessSignal || (!fProcessSignal && sideband) )
               fCutsNTuple[iNTuple]->Fill(fTMVAParticleParameters[iNTuple].data());
          }
        }
      }
    }
  }
  
  return kStOK;


}

void StKFParticleAnalysisMaker::GetDaughterParameters(const int iReader, int& iDaughterTrack, int& iDaughterParticle, KFParticle& particle)
{
  if(particle.NDaughters() == 1)
  {
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts]   = particle.GetPt();
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+1] = particle.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
    int trackId = particle.DaughterIds()[0];
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+2]   = fStKFParticleInterface->GetdEdXNSigmaPion(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+3]   = fStKFParticleInterface->GetdEdXNSigmaKaon(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+4]   = fStKFParticleInterface->GetdEdXNSigmaProton(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+5]   = fStKFParticleInterface->GetTofNSigmaPion(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+6]   = fStKFParticleInterface->GetTofNSigmaKaon(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+7]   = fStKFParticleInterface->GetTofNSigmaProton(trackId);
    
    iDaughterTrack++;
  }
  else if(particle.NDaughters() > 1)
  {
    int order[4] = {0, 1, 2, 3};
    if( particle.GetPDG() == -421 || particle.GetPDG() == -411 || particle.GetPDG() == -431 ||   
        particle.GetPDG() == -429 || particle.GetPDG() == -4122) 
    { 
      order[0] = 1; 
      order[1] = 0; 
    }
    
    for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
    {
      const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
      KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];
      //set pdg for correct order of cuts
      if(particle.GetPDG() == 521 && daughter.GetPDG() == -1) daughter.SetPDG(-421);
      if(particle.GetPDG() ==-521 && daughter.GetPDG() == -1) daughter.SetPDG( 421);
      if(particle.GetPDG() == 511 && daughter.GetPDG() == -1) daughter.SetPDG(-411);
      if(particle.GetPDG() ==-511 && daughter.GetPDG() == -1) daughter.SetPDG( 411);
        
      GetDaughterParameters(iReader, iDaughterTrack, iDaughterParticle, daughter);
    }
    
    fTMVAParticleParameters[iReader][fDaughterNames[iReader].size()*fNTrackTMVACuts + iDaughterParticle*3] = particle.Chi2()/particle.NDF();  
    
    KFParticleSIMD tempSIMDParticle(particle);
    float_v l,dl;
    KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
    tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
    fTMVAParticleParameters[iReader][fDaughterNames[iReader].size()*fNTrackTMVACuts + iDaughterParticle*3 + 1] = l[0]/dl[0];
    
    tempSIMDParticle.SetProductionVertex(pv);
    fTMVAParticleParameters[iReader][fDaughterNames[iReader].size()*fNTrackTMVACuts + iDaughterParticle*3 + 2] = 
      double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
    
    iDaughterParticle++;
  }
}

void StKFParticleAnalysisMaker::GetParticleParameters(const int iReader, KFParticle& particle)
{
  bool isBMeson = abs(particle.GetPDG()) == 511 || abs(particle.GetPDG()) == 521;
//   if( !isBMeson ) return;
  
  int iDaughterTrack = 0;
  int iDaughterParticle = 0;
  GetDaughterParameters(iReader, iDaughterTrack, iDaughterParticle, particle);

  int nDaughterParticleCut = 0;
  if(isBMeson) nDaughterParticleCut += 3;
  nDaughterParticleCut += fDaughterNames[iReader].size()*fNTrackTMVACuts;
  
  fTMVAParticleParameters[iReader][nDaughterParticleCut]   = particle.Chi2()/particle.NDF();  
  //if(_fill_lambda_tree){ld_chi2ndf = particle.Chi2()/particle.NDF();}
  
  KFParticleSIMD tempSIMDParticle(particle);
  float_v l,dl;
  KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
  fTMVAParticleParameters[iReader][nDaughterParticleCut + 1] = l[0]/dl[0];

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 2] = l[0];
  //b
  //if(_fill_lambda_tree){ld_ldl = l[0]/dl[0];}
  
  tempSIMDParticle.SetProductionVertex(pv);
  fTMVAParticleParameters[iReader][nDaughterParticleCut + 3] = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
  //if(_fill_lambda_tree){ld_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);}

  if(fIsPicoAnalysis)
    fTMVAParticleParameters[iReader][nDaughterParticleCut + 4] = fPicoDst->event()->refMult();
  else
    fTMVAParticleParameters[iReader][nDaughterParticleCut + 4] = fMuDst->event()->refMult();

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 5] = refmultcor;

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 6] = reweight;

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 7] = cent9;

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 8] = particle.GetMass();

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 9] = particle.GetPx();

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 10] = particle.GetPy();

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 11] = particle.GetPz();

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 12] = bVx;

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 13] = bVy;

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 14] = bVz;

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 15] = bVxerr;

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 16] = bVyerr;

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 17] = bVzerr;

	if(iReader>0){//cascade or omega
  int order[4] = {0, 1, 2, 3};
  const int daughterParticleIndex = particle.DaughterIds()[order[1]];
  KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 18]   = daughter.Chi2()/daughter.NDF();

  KFParticleSIMD ttempSIMDParticle(daughter);
  float_v tl,tdl;
  KFParticleSIMD tpv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  ttempSIMDParticle.GetDistanceToVertexLine(tpv, tl, tdl);

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 19] = tl[0]/tdl[0];

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 20] = tl[0];  

  ttempSIMDParticle.SetProductionVertex(tpv); 
  fTMVAParticleParameters[iReader][nDaughterParticleCut + 21] = double(ttempSIMDParticle.Chi2()[0])/double(ttempSIMDParticle.NDF()[0]);

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 22] = particle.GetPDG();
	
	}
}

Int_t StKFParticleAnalysisMaker::Finish() 
{
  if(fStoreTmvaNTuples)
  {
    TFile* curFile = gFile;
    TDirectory* curDirectory = gDirectory;
    for(int iNtuple=0; iNtuple<fNNTuples; iNtuple++)
    {
      fNTupleFile[iNtuple]->cd();
      fCutsNTuple[iNtuple]->Write();
    }
    gFile = curFile;
    gDirectory = curDirectory;
  }

//  file_out = new TFile("lambda_tree.root","RECREATE");
  file_out->cd();
//  lambda_tree->Write();
  file_out->Write();
   delete lambda_tree;
   delete file_out;

  
  return kStOK;
}

long StKFParticleAnalysisMaker::GetUniqueEventId(const int iRun, const int iEvent) const
{
  long id = 1000000000;
  return id*(iRun%1000) + iEvent;
}

int StKFParticleAnalysisMaker::GetTMVACentralityBin(int iReader, int centrality)
{
  for(unsigned int iBin=0; iBin<fTMVACentralityBins[iReader].size()-1; iBin++)
    if(centrality >= fTMVACentralityBins[iReader][iBin] && centrality < fTMVACentralityBins[iReader][iBin+1])
      return iBin;
  return -1;
}

int StKFParticleAnalysisMaker::GetTMVAPtBin(int iReader, double pt)
{
  for(unsigned int iBin=0; iBin<fTMVAPtBins[iReader].size()-1; iBin++)
    if(pt >= fTMVAPtBins[iReader][iBin] && pt < fTMVAPtBins[iReader][iBin+1])
      return iBin;
  return -1;
}

void StKFParticleAnalysisMaker::SetTMVACentralityBins(int iReader, TString bins)
{
  fTMVACentralityBins[iReader].clear();
  TString value; int firstSymbol = 0;      
  while(bins.Tokenize(value,firstSymbol,":"))
    fTMVACentralityBins[iReader].push_back(value.Atoi());
}

void StKFParticleAnalysisMaker::SetTMVAPtBins(int iReader, TString bins)
{
  fTMVAPtBins[iReader].clear();
  TString value; int firstSymbol = 0;      
  while(bins.Tokenize(value,firstSymbol,":"))
    fTMVAPtBins[iReader].push_back(value.Atof());
}

void StKFParticleAnalysisMaker::SetTMVABins(int iReader, TString centralityBins, TString ptBins)
{
  SetTMVACentralityBins(iReader, centralityBins);
  SetTMVAPtBins(iReader, ptBins);
  
  const int nCentralityBins = fTMVACentralityBins[iReader].size() - 1;
  const int nPtBins = fTMVAPtBins[iReader].size() - 1;
  
  fTMVACutFile[iReader].resize(nCentralityBins);
  fTMVACut[iReader].resize(nCentralityBins);
  fTMVAReader[iReader].resize(nCentralityBins);
  
  for(int iCentralityBin=0; iCentralityBin<nCentralityBins; iCentralityBin++)
  {
    fTMVACutFile[iReader][iCentralityBin].resize(nPtBins);
    fTMVACut[iReader][iCentralityBin].resize(nPtBins);
    fTMVAReader[iReader][iCentralityBin].resize(nPtBins);
  }
}
