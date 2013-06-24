// -*- C++ -*-
//
// Package:    CompareAmplitudes
// Class:      CompareAmplitudes
// 
/**\class CompareAmplitudes CompareAmplitudes.cc Tools/CompareAmplitudes/src/CompareAmplitudes.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Konstantinos Theofilatos
//         Created:  Fri Jan  8 18:38:05 CET 2010
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
//
// class declaration
//
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRawData/interface/EcalRawDataCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"

#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"


#include "SimCalorimetry/EcalSimAlgos/interface/EBShape.h"
#include "SimCalorimetry/EcalSimAlgos/interface/EEShape.h"


#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "TVector3.h"


#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalCleaningAlgo.h"



#include <sstream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TProfile2D.h"
#include "TStyle.h"


class CompareAmplitudes : public edm::EDAnalyzer {
   public:
      explicit CompareAmplitudes(const edm::ParameterSet&);
      ~CompareAmplitudes();
      void calculateR9(const DetId & det, const edm::Handle<EcalRecHitCollection> &hits, const edm::ESHandle<CaloTopology>  &caloTopo, 
                       float & R9var, int &gnc, int & nc);

      void calculateR9PD(const DetId & det, const edm::Handle<EcalRecHitCollection> &hits, const edm::ESHandle<CaloTopology>  &caloTopo, 
                       float & R9var, int &gnc, int & nc);


      void scan3x3(const DetId & det, const edm::Handle<EcalRecHitCollection> &hits, const edm::ESHandle<CaloTopology>  &caloTopo, 
                       float *E3x3, float *T3x3, float *J3x3, const edm::Handle<EcalUncalibratedRecHitCollection> &uhits);
      float getRatioPulsePredictionEB(float jitter, float iSample);
      float getSimPulsePredictionEB(float jitter, float iSample);

      float pulseShapeForChi2(float jitter,float iSample);
      void reset();
      enum sampleIndex{all,cleaned,good};


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      template<class T> std::string any2string(T i);


      edm::Service<TFileService> fTFileService_; 

      edm::InputTag EBDigiCollection_;
      edm::InputTag EEDigiCollection_;
      edm::Handle<EBDigiCollection> EBdigis_;
      edm::Handle<EEDigiCollection> EEdigis_;

      edm::InputTag EBURecHitCollection_;
      edm::InputTag EEURecHitCollection_;
      edm::Handle<EcalUncalibratedRecHitCollection> EBURecHits_;
      edm::Handle<EcalUncalibratedRecHitCollection> EEURecHits_;

      edm::InputTag EBRecHitCollection_;
      edm::InputTag EERecHitCollection_;
      edm::Handle<EcalRecHitCollection> EBRecHits_;
      edm::Handle<EcalRecHitCollection> EERecHits_;


      edm::InputTag EBFitRecHitCollection_;
      edm::InputTag EEFitRecHitCollection_;
      edm::Handle<EcalRecHitCollection> EBFitRecHits_;
      edm::Handle<EcalRecHitCollection> EEFitRecHits_;

      TH1F *energy_EB[3];
      TH1F *topo_EB[3];
      TH1F *chrono_EB[3];
      TH1F *shape_it_EB[3];
      TH1F *shape_oot_EB[3];
      TH2F *shape_itVSenergy_EB[3];
      TH2F *shape_itVStopo_EB[3];
      TH2F *shape_itVSchrono_EB[3];
      TH2F *shape_ootVSenergy_EB[3];
      TH2F *shape_ootVStopo_EB[3];
      TH2F *shape_ootVSchrono_EB[3];
      TH2F *topoVSenergy_EB[3];
      TH2F *chronoVSenergy_EB[3];

      TH1F *energy_EE[3];
      TH1F *topo_EE[3];
      TH1F *chrono_EE[3];
      TH1F *shape_it_EE[3];
      TH1F *shape_oot_EE[3];
      TH2F *shape_itVSenergy_EE[3];
      TH2F *shape_itVStopo_EE[3];
      TH2F *shape_itVSchrono_EE[3];
      TH2F *shape_ootVSenergy_EE[3];
      TH2F *shape_ootVStopo_EE[3];
      TH2F *shape_ootVSchrono_EE[3];
      TH2F *topoVSenergy_EE[3];
      TH2F *chronoVSenergy_EE[3];


      int topoBins ;
      float topoMin  ;
      float topoMax  ;

      int chronoBins ;
      float chronoMin  ;
      float chronoMax  ;

      int shapeBins ;
      float shapeMin  ;
      float shapeMax  ;

      int energyBins ;
      float energyMin  ;
      float energyMax  ;


  
     // TFile *fp_;
      TTree *myTree;
      float ebEnergyMin_;
      float eeEnergyMin_;
      float gncMinEnergy_;

      // ==== global rechit       

      float amp_;
      float ped_;
      float chi2_;

      float calibE_;
      float jitter_;
      float chi2Fixed_;
      float chi2Sliding_;
      float time_;
      float E3x3array_[9];
      float T3x3array_[9];
      float J3x3array_[9];
      float lambda3_;

      float sampleTime_[10];
      float mgpaSample_[10];
      float chi2Digi_[10];

      float R9var_;

      float rechitPt_;
      float rechitEta_;
      float rechitPhi_;

      unsigned int rechitFlag_;

      float calibEFit_;
      float timeError_;
      int isRecovered_;

      int gnc_;
      int nc_;
      int ncRelaxed_;

      int nhits_;
      int kPoorReco_;
      int kPoorReco_Uncalib_;
      int kGood_;
      int ecalSL_;
      int kOutOfTime_;


      int sampleGainId_[10];
      float sampleGainRatio_[10];
      float amplitudeAsynchFit_[10];


      int kDiWeird_;
      int kWeird_;

      //  ===============

      float ped1_;
      float ped2_;
      float ped3_;
      float pedDB_;
      float pedRMSDB_;
    

      float timeIC_;
      int gain_switch_;
      int iphi_;
      int ieta_;
      int iz_;
      int subDet_;

      int runNumber_;
      int lumi_;
      int eventNumber_;

      // =================


      TTree* global_data_;
      float mgpaSimPulseEB_[250];
      float mgpaSimPulseEE_[250];
      float mgpaRatioPulseEB_[250];
      float mgpaRatioPulseEE_[250];
      float mgpaTime_[250];

};

CompareAmplitudes::CompareAmplitudes(const edm::ParameterSet& iConfig)
{
    EBDigiCollection_ = iConfig.getParameter<edm::InputTag>("EBDigiCollection");
    EEDigiCollection_ = iConfig.getParameter<edm::InputTag>("EEDigiCollection");

    EBURecHitCollection_ = iConfig.getParameter<edm::InputTag>("EBURecHitCollection");
    EEURecHitCollection_ = iConfig.getParameter<edm::InputTag>("EEURecHitCollection");

    EBRecHitCollection_ = iConfig.getParameter<edm::InputTag>("EBRecHitCollection");
    EERecHitCollection_ = iConfig.getParameter<edm::InputTag>("EERecHitCollection");


    ebEnergyMin_ = iConfig.getUntrackedParameter<double>("ebEnergyMin",0.3);
    eeEnergyMin_ = iConfig.getUntrackedParameter<double>("eeEnergyMin",0.3);
    gncMinEnergy_ = iConfig.getUntrackedParameter<double>("gncMinEnergy",0.8);

}


CompareAmplitudes::~CompareAmplitudes()
{
}


void CompareAmplitudes::calculateR9(const DetId & det, const edm::Handle<EcalRecHitCollection> &hits, const edm::ESHandle<CaloTopology>  &caloTopo,
                                    float &R9var, int &gnc, int & nc)
{

  if(hits->find(det) !=hits->end())
  {

    CaloNavigator<DetId> cursor = CaloNavigator<DetId>(det,caloTopo->getSubdetectorTopology(det));
    int side_=3;

    EcalRecHit centralHit = *hits->find(det);
  
    float E1 =  0;
    float E3x3 = 0;
    int gnc_tmp=0;
    int nc_tmp=0;

    E1  =  centralHit.energy(); //consider only good rechits

    for(int j=side_/2; j>=-side_/2; --j)
    {
      for(int i=-side_/2; i<=side_/2; ++i)
      {
	cursor.home();
	cursor.offsetBy(i,j);
	
	if(hits->find(*cursor)!=hits->end())
	{
	  EcalRecHit tmpHit = *hits->find(*cursor);
	  {
	      E3x3 += tmpHit.energy();
	      nc_tmp++;
	      if(tmpHit.energy()>gncMinEnergy_)gnc_tmp++;
	  }
	}
      }
    }
    R9var = -1;
    if(E3x3!=0)R9var = E1/E3x3;
    gnc = gnc_tmp;
    nc = nc_tmp;

  }

}


void CompareAmplitudes::reset()
{

    amp_=0;
    ped_=0;
    chi2_=0;

    calibE_=0;
    jitter_=0;
    chi2Fixed_=0;
    chi2Sliding_=0;
    time_=0;
    for(int i=0;i<9;i++)E3x3array_[i]=0;
    for(int i=0;i<9;i++)T3x3array_[i]=0;
    for(int i=0;i<9;i++)J3x3array_[i]=0;
    lambda3_=0;

    for(int i=0;i<10;i++)sampleTime_[i]=0;
    for(int i=0;i<10;i++)mgpaSample_[i]=0;
    for(int i=0;i<10;i++)chi2Digi_[i]=0;

    R9var_=0;

    rechitPt_=0;
    rechitEta_=0;
    rechitPhi_=0;

    rechitFlag_=0;

    calibEFit_=0;
    timeError_=0;
    isRecovered_=0;

    gnc_=0;
    nc_=0;
    ncRelaxed_=0;

    nhits_=0;
    kPoorReco_=22;
    kPoorReco_Uncalib_=22;
    kGood_=22;
    kOutOfTime_=22;
    kDiWeird_=22;
    kWeird_=22;
    ecalSL_=0;


    for(int i=0;i<10;i++)sampleGainId_[i]=0;
    for(int i=0;i<10;i++)sampleGainRatio_[i]=0;
    for(int i=0;i<10;i++)amplitudeAsynchFit_[i]=0;




}

void CompareAmplitudes::calculateR9PD(const DetId & det, const edm::Handle<EcalRecHitCollection> &hits, const edm::ESHandle<CaloTopology>  &caloTopo,
                                    float &R9var, int &gnc, int & nc)
{

  if(hits->find(det) !=hits->end())
  {

    CaloNavigator<DetId> cursor = CaloNavigator<DetId>(det,caloTopo->getSubdetectorTopology(det));
    int side_=3;

    EcalRecHit centralHit = *hits->find(det);
  
    float E1 =  0;
    float E3x3 = 0;
    int gnc_tmp=0;
    int nc_tmp=0;

    E1  =  centralHit.energy(); //consider only good rechits

    for(int j=side_/2; j>=-side_/2; --j)
    {
      for(int i=-side_/2; i<=side_/2; ++i)
      {
	cursor.home();
	cursor.offsetBy(i,j);
	
	if(hits->find(*cursor)!=hits->end())
	{
	  EcalRecHit tmpHit = *hits->find(*cursor);
//	  if(tmpHit.recoFlag()==EcalRecHit::kGood) // consider only good energy
	  {
	      if(tmpHit.energy()>0.0)E3x3 += tmpHit.energy();
	      nc_tmp++;
	      if(tmpHit.energy()>gncMinEnergy_)gnc_tmp++;
	  }
	}
      }
    }
    R9var = -1;
    if(E3x3!=0)R9var = E1/E3x3;
    gnc = gnc_tmp;
    nc = nc_tmp;

  }

}

float CompareAmplitudes::getRatioPulsePredictionEB(float jitter, float iSample)
{
    double amplitudeFitParameters[] = {1.138,1.652};
    double alpha = amplitudeFitParameters[0];
    double beta = amplitudeFitParameters[1];

    double Tmax = 5.0 + jitter;

    double fi = 0;
    double termOne =  1 + (double(iSample) - Tmax) / (alpha * beta);
    if (termOne > 1.e-5) fi  = exp(alpha * log(termOne)) * exp(-(iSample - Tmax) / beta);

    return fi;

}

float CompareAmplitudes::pulseShapeForChi2(float jitter,float iSample)
{
    float fi=0;
    float t = (iSample - jitter)*25.0;

    TF1 *landau_thanks = new TF1("landau_thanks","landau",62.5,98);
    landau_thanks->SetParameter(0,4.99002);
    landau_thanks->SetParError(0,0.0239832);
    landau_thanks->SetParameter(1,113.497);
    landau_thanks->SetParError(1,0.209147);
    landau_thanks->SetParameter(2,14.5021);
    landau_thanks->SetParError(2,0.0768403);

    TF1 *sasha_fit = new TF1("sasha_fit","exp(1.138 * log(1+ (x-125)/(1.138*1.652*25.0))) * exp(-(x - 125.0) /(1.652*25.0) )",98,250);

    if(t<=98)fi=landau_thanks->Eval(t);
    if(t>98)fi=sasha_fit->Eval(t);
    if(t<62.5)fi=0;

    delete landau_thanks;
    delete sasha_fit;

    return fi;
}

float CompareAmplitudes::getSimPulsePredictionEB(float jitter, float iSample)
{
    double fi =0;
    EBShape mgpaSimPulseShapeEB;
    
    double ADC_clock = 25.0; // 25 ns
    double risingTimeEB = mgpaSimPulseShapeEB.timeToRise();
    double tzeroEB = risingTimeEB  - 5*ADC_clock;  // 5 samples before the peak
    fi = (mgpaSimPulseShapeEB)(tzeroEB  + float(iSample));


    return fi;

}

void CompareAmplitudes::scan3x3(const DetId & det, const edm::Handle<EcalRecHitCollection> &hits, const edm::ESHandle<CaloTopology>  &caloTopo,
                                    float *E3x3, float *T3x3, float *J3x3, const edm::Handle<EcalUncalibratedRecHitCollection> &uhits)
{
  for(int i=0;i<9;i++)E3x3[i]=0;
  for(int i=0;i<9;i++)T3x3[i]=0;
  for(int i=0;i<9;i++)J3x3[i]=0;

  if(hits->find(det) !=hits->end())
  {

    CaloNavigator<DetId> cursor = CaloNavigator<DetId>(det,caloTopo->getSubdetectorTopology(det));
    int side_=3;

  

    int counter=0;
    for(int j=side_/2; j>=-side_/2; --j)
    {
      for(int i=-side_/2; i<=side_/2; ++i)
      {
	cursor.home();
	cursor.offsetBy(i,j);
	
	if(hits->find(*cursor)!=hits->end())
	{
	  EcalRecHit tmpHit = *hits->find(*cursor);
	  E3x3[counter]=tmpHit.energy();
	  T3x3[counter]=tmpHit.time();
	}
	if(uhits->find(*cursor)!=uhits->end())
	{
	  EcalUncalibratedRecHit tmpUHit = *uhits->find(*cursor);
	  J3x3[counter]=tmpUHit.jitter();
	}

	counter++;
      }
    }
  }
}



// ------------ method called to for each event  ------------
void
CompareAmplitudes::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;


//   std::cout << "::: Trying to fetch the collections " << std::endl;
  iEvent.getByLabel(EBDigiCollection_, EBdigis_);
  iEvent.getByLabel(EEDigiCollection_, EEdigis_);

  iEvent.getByLabel(EBURecHitCollection_, EBURecHits_);
  iEvent.getByLabel(EEURecHitCollection_, EEURecHits_);

  iEvent.getByLabel(EBRecHitCollection_, EBRecHits_);
  iEvent.getByLabel(EERecHitCollection_, EERecHits_);

  ESHandle<CaloTopology> caloTopo;
  iSetup.get<CaloTopologyRecord>().get(caloTopo);
  
  edm::ESHandle<EcalPedestals> peds;
  iSetup.get<EcalPedestalsRcd>().get(peds);

  edm::ESHandle<EcalTimeCalibConstants> itime;
  iSetup.get<EcalTimeCalibConstantsRcd>().get(itime);

  EcalUncalibratedRecHitCollection::const_iterator urechitIt;
  EcalRecHitCollection::const_iterator rechitIt;

// Get Geometry
  edm::ESHandle<CaloGeometry> geometry;
  iSetup.get<CaloGeometryRecord>().get(geometry);

  edm::ESHandle<EcalSeverityLevelAlgo> ecalSevLvlAlgoHndl;
  iSetup.get<EcalSeverityLevelAlgoRcd>().get(ecalSevLvlAlgoHndl);
  const EcalSeverityLevelAlgo* ecalSevLvlAlgo = ecalSevLvlAlgoHndl.product();

   // common setup
  edm::ESHandle<EcalGainRatios>  gains;
  iSetup.get<EcalGainRatiosRcd>().get(gains);
  double gainRatios[3];


  TVector3 metVector[20];
  float sumEt[20];
  for(int i =0;i<20;i++)
  {
      metVector[i].SetXYZ(0.0,0.0,0.0);
      sumEt[i] = 0.0;
  }


  for(EBDigiCollection::const_iterator It = EBdigis_->begin(); It != EBdigis_->end(); ++It)
  {
    reset();

    EBDetId ebDetId  = (*It).id();
    int iphi = ebDetId.iphi();
    int ieta = ebDetId.ieta();
    int iz = ebDetId.zside();
    int subDet = 1;

    urechitIt = EBURecHits_->find(ebDetId);
    rechitIt = EBRecHits_->find(ebDetId);


    bool urechitFound = false;
    if(urechitIt!=EBURecHits_->end())urechitFound=true;

    bool rechitFound = false;
    if(rechitIt!=EBRecHits_->end())rechitFound=true;


    // compute the right bin of the pulse shape using time calibration constants
    EcalTimeCalibConstantMap::const_iterator itTime = itime->find( ebDetId );
    EcalTimeCalibConstant itimeconst = 0;
    if( itTime != itime->end() ) {
	      itimeconst = (*itTime);
    } else {  
	      edm::LogError("EcalRecHitError") << "No time intercalib const found for xtal "
	      << ebDetId.rawId()
	      << "! something wrong with EcalTimeCalibConstants in your DB? ";
    }



    const EcalMGPAGainRatio * aGain = 0;
    const EcalPedestals::Item * aped = 0;
    unsigned int hashedIndex = ebDetId.hashedIndex();
    aped = &peds->barrel(hashedIndex);
    aGain = &gains->barrel(hashedIndex);

    double pedVec[3];
    double pedRMSVec[3];

    pedVec[0] = aped->mean_x12;
    pedVec[1] = aped->mean_x6;
    pedVec[2] = aped->mean_x1;

    pedRMSVec[0] = aped->rms_x12;
    pedRMSVec[1] = aped->rms_x6;
    pedRMSVec[2] = aped->rms_x1;

    gainRatios[0] = 1.;
    gainRatios[1] = aGain->gain12Over6();
    gainRatios[2] = aGain->gain6Over1()*aGain->gain12Over6();


    bool gainSwitch = false;

    EcalDataFrame df;
    df = (*It);

    double minDigi=1.e+10;
    double maxDigi=1.e-10;
    for(int iSample = 0; iSample < 10; iSample++) 
    {
      int gainId = df.sample(iSample).gainId();
      if (gainId != 1) {gainSwitch=true;}
      double currentDigi = df.sample(iSample).adc();
      
      if(currentDigi>maxDigi)maxDigi=currentDigi;
      if(currentDigi<minDigi)minDigi=currentDigi;
    }

  

    if(urechitFound && rechitFound)
    {

      float R9var= 0;
      int gnc =0;
      int nc =0;
      int ncRelaxed =0;
      calculateR9(ebDetId, EBRecHits_, caloTopo, R9var , gnc, nc);



      int nhits = int(EBURecHits_->size()); 


      // ======= Fill the tree

      // ==== global rechit
      amp_ = urechitIt->amplitude();
      ped_ = urechitIt->pedestal();
      jitter_ = urechitIt->jitter();
      chi2_ = urechitIt->chi2();

      calibE_ = rechitIt->energy();
      calibEFit_ = rechitIt->outOfTimeEnergy();
      chi2Fixed_ = rechitIt->chi2();
      chi2Sliding_ = rechitIt->outOfTimeChi2();
      time_ = rechitIt->time();
      timeError_ = rechitIt->timeError();
      R9var_ = R9var; 
      isRecovered_ = rechitIt->isRecovered();

      kGood_ = rechitIt->checkFlag(EcalRecHit::kGood);
      kOutOfTime_ = rechitIt->checkFlag(EcalRecHit::kOutOfTime);
      kPoorReco_ = rechitIt->checkFlag(EcalRecHit::kPoorReco);
      kPoorReco_Uncalib_ = urechitIt->checkFlag(EcalUncalibratedRecHit::kPoorReco) ? 1:0;
      kWeird_ = rechitIt->checkFlag(EcalRecHit::kWeird);
      kDiWeird_ = rechitIt->checkFlag(EcalRecHit::kDiWeird);
      ecalSL_ =  ecalSevLvlAlgo->severityLevel(*rechitIt);// 

      gnc_ = gnc;
      nc_ = nc;
      ncRelaxed_ = ncRelaxed;
   
      rechitFlag_ = rechitIt->recoFlag();

      nhits_ = nhits;


      scan3x3(ebDetId, EBRecHits_, caloTopo, E3x3array_, T3x3array_ , J3x3array_, EBURecHits_);




      // ===
      iphi_ = iphi;
      ieta_ = ieta;
      iz_ = iz;
      timeIC_ = (float) itimeconst;	    

      gain_switch_ = gainSwitch ? 1:0;
      subDet_ = subDet;
      runNumber_ = iEvent.id().run();
      lumi_ = iEvent.id().luminosityBlock();
      eventNumber_ = iEvent.id().event();

      pedDB_ = pedVec[0];
      pedRMSDB_ = pedRMSVec[0];
      ped3_ = (1/3.0)*(df.sample(0).adc()+df.sample(1).adc()+df.sample(2).adc());
      ped2_ = (1/2.0)*(df.sample(0).adc()+df.sample(1).adc());
      ped1_ = df.sample(0).adc();

      for(int iSample = 0; iSample < 10; iSample++)
      {
        mgpaSample_[iSample] = float(df.sample(iSample).adc());
	sampleGainId_[iSample] = df.sample(iSample).gainId();
     	int effGain = df.sample(iSample).gainId()-1;
        if( df.sample(iSample).gainId()==0)effGain=3; // gainId=0 for saturation
        sampleGainRatio_[iSample] = gainRatios[effGain];
	sampleTime_[iSample] = iSample*25.0;
      }

      float chi2GS = (((mgpaSample_[3]-ped1_)/amp_-(sampleTime_[3] + timeIC_)*0.00366827 + 0.257046)/0.0045);
      float chi2SP = ((mgpaSample_[3]-ped1_)/(amp_*0.00694));
//      lambda3_ = (((mgpaSample_[3]-ped1_)/amp_-(sampleTime_[3] + timeIC_)*0.00366827 + 0.257046)/0.0045)**2 - ((mgpaSample_[3]-ped1)/(amp_*0.00694))**2;
      lambda3_ = chi2GS*chi2GS - chi2SP*chi2SP;

     // --- store geometry related
      const GlobalPoint p ( geometry->getPosition( ebDetId ) ) ;
      TVector3 hitPos(p.x(),p.y(),p.z());
      hitPos *= 1.0/hitPos.Mag();
      hitPos *= rechitIt->energy();

      rechitPt_ =  hitPos.Pt();
      if(hitPos.Pt()>0.01)rechitEta_ =  hitPos.Eta();
      rechitPhi_ = hitPos.Phi() ;


      if(calibE_>ebEnergyMin_ || calibE_<-1)myTree->Fill();

      // fill histograms

   float epsilon = 0.0001;
   float topoValue = R9var_; if(topoValue<topoMin)topoValue=topoMin+epsilon;if(topoValue>topoMax)topoValue=topoMax-epsilon;
   float chronoValue = time_; if(chronoValue<chronoMin)chronoValue=chronoMin+epsilon;if(chronoValue>chronoMax)chronoValue=chronoMax-epsilon;
   float shapeValue_it = chi2_; if(shapeValue_it<shapeMin)shapeValue_it=shapeMin+epsilon;if(shapeValue_it>shapeMax)shapeValue_it=shapeMax-epsilon;
   float shapeValue_oot = chi2Sliding_;if(shapeValue_oot<shapeMin)shapeValue_oot=shapeMin+epsilon;if(shapeValue_oot>shapeMax)shapeValue_oot=shapeMax-epsilon;
   float energyValue = calibE_; if(energyValue<energyMin)energyValue=energyMin+epsilon;if(energyValue>energyMax)energyValue=energyMax-epsilon;


      if(!isRecovered_ && !gain_switch_)
      { 
	energy_EB[all]->Fill(energyValue);
	if(energyValue>15)topo_EB[all]->Fill(topoValue); 
	if(energyValue>15)chrono_EB[all]->Fill(chronoValue);
	if(energyValue>15)shape_it_EB[all]->Fill(shapeValue_it);
	shape_itVSenergy_EB[all]->Fill(energyValue,shapeValue_it);
	if(energyValue>15)shape_itVStopo_EB[all]->Fill(topoValue,shapeValue_it);
	if(energyValue>15)shape_itVSchrono_EB[all]->Fill(chronoValue,shapeValue_it);
	if(energyValue>15)shape_oot_EB[all]->Fill(shapeValue_oot);
	shape_ootVSenergy_EB[all]->Fill(energyValue,shapeValue_oot);
	if(energyValue>15)shape_ootVStopo_EB[all]->Fill(topoValue,shapeValue_oot);
	if(energyValue>15)shape_ootVSchrono_EB[all]->Fill(chronoValue,shapeValue_oot);
	topoVSenergy_EB[all]->Fill(energyValue,topoValue);
	chronoVSenergy_EB[all]->Fill(energyValue,chronoValue);
      }
      if( !isRecovered_ && !gain_switch_ && ecalSL_==0)
      { 
        energy_EB[cleaned]->Fill(energyValue);
        if(energyValue>15)topo_EB[cleaned]->Fill(topoValue);
        if(energyValue>15)chrono_EB[cleaned]->Fill(chronoValue);
        if(energyValue>15)shape_it_EB[cleaned]->Fill(shapeValue_it);
        shape_itVSenergy_EB[cleaned]->Fill(energyValue,shapeValue_it);
        if(energyValue>15)shape_itVStopo_EB[cleaned]->Fill(topoValue,shapeValue_it);
        if(energyValue>15)shape_itVSchrono_EB[cleaned]->Fill(chronoValue,shapeValue_it);
        if(energyValue>15)shape_oot_EB[cleaned]->Fill(shapeValue_oot);
        shape_ootVSenergy_EB[cleaned]->Fill(energyValue,shapeValue_oot);
        if(energyValue>15)shape_ootVStopo_EB[cleaned]->Fill(topoValue,shapeValue_oot);
        if(energyValue>15)shape_ootVSchrono_EB[cleaned]->Fill(chronoValue,shapeValue_oot);
        topoVSenergy_EB[cleaned]->Fill(energyValue,topoValue);
        chronoVSenergy_EB[cleaned]->Fill(energyValue,chronoValue);
      }
      if(!isRecovered_ && !gain_switch_ && ecalSL_==0 && !kPoorReco_)
      { 
        energy_EB[good]->Fill(energyValue);
        if(energyValue>15)topo_EB[good]->Fill(topoValue);
        if(energyValue>15)chrono_EB[good]->Fill(chronoValue);
        if(energyValue>15)shape_it_EB[good]->Fill(shapeValue_it);
        shape_itVSenergy_EB[good]->Fill(energyValue,shapeValue_it);
        if(energyValue>15)shape_itVStopo_EB[good]->Fill(topoValue,shapeValue_it);
        if(energyValue>15)shape_itVSchrono_EB[good]->Fill(chronoValue,shapeValue_it);
        if(energyValue>15)shape_oot_EB[good]->Fill(shapeValue_oot);
        shape_ootVSenergy_EB[good]->Fill(energyValue,shapeValue_oot);
        if(energyValue>15)shape_ootVStopo_EB[good]->Fill(topoValue,shapeValue_oot);
        if(energyValue>15)shape_ootVSchrono_EB[good]->Fill(chronoValue,shapeValue_oot);
        topoVSenergy_EB[good]->Fill(energyValue,topoValue);
        chronoVSenergy_EB[good]->Fill(energyValue,chronoValue);
      }


    } // end of IF urehit && rechit true
  } // end of rechit loop



  for(EEDigiCollection::const_iterator It = EEdigis_->begin(); It != EEdigis_->end(); ++It)
  {
    reset();

    EEDetId eeDetId  = (*It).id();
    int iphi = eeDetId.ix();
    int ieta = eeDetId.iy();
    int iz = eeDetId.zside();
    int subDet = 2;

    urechitIt = EEURecHits_->find(eeDetId);
    rechitIt = EERecHits_->find(eeDetId);


    bool urechitFound = false;
    if(urechitIt!=EEURecHits_->end())urechitFound=true;

    bool rechitFound = false;
    if(rechitIt!=EERecHits_->end())rechitFound=true;


    // compute the right bin of the pulse shape using time calibration constants
    EcalTimeCalibConstantMap::const_iterator itTime = itime->find( eeDetId );
    EcalTimeCalibConstant itimeconst = 0;
    if( itTime != itime->end() ) {
	      itimeconst = (*itTime);
    } else {  
	      edm::LogError("EcalRecHitError") << "No time intercalib const found for xtal "
	      << eeDetId.rawId()
	      << "! something wrong with EcalTimeCalibConstants in your DB? ";
    }



    const EcalMGPAGainRatio * aGain = 0;
    const EcalPedestals::Item * aped = 0;
    unsigned int hashedIndex = eeDetId.hashedIndex();
    aped = &peds->barrel(hashedIndex);
    aGain = &gains->barrel(hashedIndex);

    double pedVec[3];
    double pedRMSVec[3];

    pedVec[0] = aped->mean_x12;
    pedVec[1] = aped->mean_x6;
    pedVec[2] = aped->mean_x1;

    pedRMSVec[0] = aped->rms_x12;
    pedRMSVec[1] = aped->rms_x6;
    pedRMSVec[2] = aped->rms_x1;

    gainRatios[0] = 1.;
    gainRatios[1] = aGain->gain12Over6();
    gainRatios[2] = aGain->gain6Over1()*aGain->gain12Over6();


    bool gainSwitch = false;

    EcalDataFrame df;
    df = (*It);

    double minDigi=1.e+10;
    double maxDigi=1.e-10;
    for(int iSample = 0; iSample < 10; iSample++) 
    {
      int gainId = df.sample(iSample).gainId();
      if (gainId != 1) {gainSwitch=true;}
      double currentDigi = df.sample(iSample).adc();
      
      if(currentDigi>maxDigi)maxDigi=currentDigi;
      if(currentDigi<minDigi)minDigi=currentDigi;
    }

  

    if(urechitFound && rechitFound)
    {

      float R9var= 0;
      int gnc =0;
      int nc =0;
      int ncRelaxed =0;
      calculateR9(eeDetId, EERecHits_, caloTopo, R9var , gnc, nc);



      int nhits = int(EEURecHits_->size()); 


      // ======= Fill the tree

      // ==== global rechit
      amp_ = urechitIt->amplitude();
      ped_ = urechitIt->pedestal();
      jitter_ = urechitIt->jitter();
      chi2_ = urechitIt->chi2();

      calibE_ = rechitIt->energy();
      calibEFit_ = rechitIt->outOfTimeEnergy();
      chi2Fixed_ = rechitIt->chi2();
      chi2Sliding_ = rechitIt->outOfTimeChi2();
      time_ = rechitIt->time();
      timeError_ = rechitIt->timeError();
      R9var_ = R9var; 
      isRecovered_ = rechitIt->isRecovered();

      kGood_ = rechitIt->checkFlag(EcalRecHit::kGood);
      kOutOfTime_ = rechitIt->checkFlag(EcalRecHit::kOutOfTime);
      kPoorReco_ = rechitIt->checkFlag(EcalRecHit::kPoorReco);
      kPoorReco_Uncalib_ = urechitIt->checkFlag(EcalUncalibratedRecHit::kPoorReco) ? 1:0;
      kWeird_ = rechitIt->checkFlag(EcalRecHit::kWeird);
      kDiWeird_ = rechitIt->checkFlag(EcalRecHit::kDiWeird);
      ecalSL_ =  ecalSevLvlAlgo->severityLevel(*rechitIt);// 

      gnc_ = gnc;
      nc_ = nc;
      ncRelaxed_ = ncRelaxed;
   
      rechitFlag_ = rechitIt->recoFlag();

      nhits_ = nhits;


      scan3x3(eeDetId, EERecHits_, caloTopo, E3x3array_, T3x3array_ , J3x3array_, EEURecHits_);




      // ===
      iphi_ = iphi;
      ieta_ = ieta;
      iz_ = iz;
      timeIC_ = (float) itimeconst;	    

      gain_switch_ = gainSwitch ? 1:0;
      subDet_ = subDet;
      runNumber_ = iEvent.id().run();
      lumi_ = iEvent.id().luminosityBlock();
      eventNumber_ = iEvent.id().event();

      pedDB_ = pedVec[0];
      pedRMSDB_ = pedRMSVec[0];
      ped3_ = (1/3.0)*(df.sample(0).adc()+df.sample(1).adc()+df.sample(2).adc());
      ped2_ = (1/2.0)*(df.sample(0).adc()+df.sample(1).adc());
      ped1_ = df.sample(0).adc();

      for(int iSample = 0; iSample < 10; iSample++)
      {
        mgpaSample_[iSample] = float(df.sample(iSample).adc());
	sampleGainId_[iSample] = df.sample(iSample).gainId();
	int effGain = df.sample(iSample).gainId()-1;
	if( df.sample(iSample).gainId()==0)effGain=3; // gainId=0 for saturation
	sampleGainRatio_[iSample] = gainRatios[effGain];
	sampleTime_[iSample] = iSample*25.0;
      }

      float chi2GS = (((mgpaSample_[3]-ped1_)/amp_-(sampleTime_[3] + timeIC_)*0.00366827 + 0.257046)/0.0045);
      float chi2SP = ((mgpaSample_[3]-ped1_)/(amp_*0.00694));
//      lambda3_ = (((mgpaSample_[3]-ped1_)/amp_-(sampleTime_[3] + timeIC_)*0.00366827 + 0.257046)/0.0045)**2 - ((mgpaSample_[3]-ped1)/(amp_*0.00694))**2;
      lambda3_ = chi2GS*chi2GS - chi2SP*chi2SP;

     // --- store geometry related
      const GlobalPoint p ( geometry->getPosition( eeDetId ) ) ;
      TVector3 hitPos(p.x(),p.y(),p.z());
      hitPos *= 1.0/hitPos.Mag();
      hitPos *= rechitIt->energy();

      rechitPt_ =  hitPos.Pt();
      if(hitPos.Pt()>0.01)rechitEta_ =  hitPos.Eta();
      rechitPhi_ = hitPos.Phi() ;


      if(calibE_>eeEnergyMin_ || calibE_<-1)myTree->Fill();

      // fill histograms

   float epsilon = 0.0001;
   float topoValue = R9var_; if(topoValue<topoMin)topoValue=topoMin+epsilon;if(topoValue>topoMax)topoValue=topoMax-epsilon;
   float chronoValue = time_; if(chronoValue<chronoMin)chronoValue=chronoMin+epsilon;if(chronoValue>chronoMax)chronoValue=chronoMax-epsilon;
   float shapeValue_it = chi2_; if(shapeValue_it<shapeMin)shapeValue_it=shapeMin+epsilon;if(shapeValue_it>shapeMax)shapeValue_it=shapeMax-epsilon;
   float shapeValue_oot = chi2Sliding_;if(shapeValue_oot<shapeMin)shapeValue_oot=shapeMin+epsilon;if(shapeValue_oot>shapeMax)shapeValue_oot=shapeMax-epsilon;
   float energyValue = calibE_; if(energyValue<energyMin)energyValue=energyMin+epsilon;if(energyValue>energyMax)energyValue=energyMax-epsilon;


      if(!isRecovered_ && !gain_switch_)
      { 
	energy_EE[all]->Fill(energyValue);
	if(energyValue>15)topo_EE[all]->Fill(topoValue); 
	if(energyValue>15)chrono_EE[all]->Fill(chronoValue);
	if(energyValue>15)shape_it_EE[all]->Fill(shapeValue_it);
	shape_itVSenergy_EE[all]->Fill(energyValue,shapeValue_it);
	if(energyValue>15)shape_itVStopo_EE[all]->Fill(topoValue,shapeValue_it);
	if(energyValue>15)shape_itVSchrono_EE[all]->Fill(chronoValue,shapeValue_it);
	if(energyValue>15)shape_oot_EE[all]->Fill(shapeValue_oot);
	shape_ootVSenergy_EE[all]->Fill(energyValue,shapeValue_oot);
	if(energyValue>15)shape_ootVStopo_EE[all]->Fill(topoValue,shapeValue_oot);
	if(energyValue>15)shape_ootVSchrono_EE[all]->Fill(chronoValue,shapeValue_oot);
	topoVSenergy_EE[all]->Fill(energyValue,topoValue);
	chronoVSenergy_EE[all]->Fill(energyValue,chronoValue);
      }
      if( !isRecovered_ && !gain_switch_ && ecalSL_==0)
      { 
        energy_EE[cleaned]->Fill(energyValue);
        if(energyValue>15)topo_EE[cleaned]->Fill(topoValue);
        if(energyValue>15)chrono_EE[cleaned]->Fill(chronoValue);
        if(energyValue>15)shape_it_EE[cleaned]->Fill(shapeValue_it);
        shape_itVSenergy_EE[cleaned]->Fill(energyValue,shapeValue_it);
        if(energyValue>15)shape_itVStopo_EE[cleaned]->Fill(topoValue,shapeValue_it);
        if(energyValue>15)shape_itVSchrono_EE[cleaned]->Fill(chronoValue,shapeValue_it);
        if(energyValue>15)shape_oot_EE[cleaned]->Fill(shapeValue_oot);
        shape_ootVSenergy_EE[cleaned]->Fill(energyValue,shapeValue_oot);
        if(energyValue>15)shape_ootVStopo_EE[cleaned]->Fill(topoValue,shapeValue_oot);
        if(energyValue>15)shape_ootVSchrono_EE[cleaned]->Fill(chronoValue,shapeValue_oot);
        topoVSenergy_EE[cleaned]->Fill(energyValue,topoValue);
        chronoVSenergy_EE[cleaned]->Fill(energyValue,chronoValue);
      }
      if(!isRecovered_ && !gain_switch_ && ecalSL_==0 && !kPoorReco_)
      { 
        energy_EE[good]->Fill(energyValue);
        if(energyValue>15)topo_EE[good]->Fill(topoValue);
        if(energyValue>15)chrono_EE[good]->Fill(chronoValue);
        if(energyValue>15)shape_it_EE[good]->Fill(shapeValue_it);
        shape_itVSenergy_EE[good]->Fill(energyValue,shapeValue_it);
        if(energyValue>15)shape_itVStopo_EE[good]->Fill(topoValue,shapeValue_it);
        if(energyValue>15)shape_itVSchrono_EE[good]->Fill(chronoValue,shapeValue_it);
        if(energyValue>15)shape_oot_EE[good]->Fill(shapeValue_oot);
        shape_ootVSenergy_EE[good]->Fill(energyValue,shapeValue_oot);
        if(energyValue>15)shape_ootVStopo_EE[good]->Fill(topoValue,shapeValue_oot);
        if(energyValue>15)shape_ootVSchrono_EE[good]->Fill(chronoValue,shapeValue_oot);
        topoVSenergy_EE[good]->Fill(energyValue,topoValue);
        chronoVSenergy_EE[good]->Fill(energyValue,chronoValue);
      }


    } // end of IF urehit && rechit true
  } // end of rechit loop


}

template<class T>
std::string CompareAmplitudes::any2string(T i)
{
        std::ostringstream buffer;
        buffer << i;
        return buffer.str();
}



// ------------ method called once each job just before starting event loop  ------------
void CompareAmplitudes::beginJob()
{



  myTree = fTFileService_->make<TTree>("events","events");
  // === global rechit
  myTree->Branch("amp",&amp_,"amp/F");
  myTree->Branch("chi2",&chi2_,"chi2/F");
  myTree->Branch("jitter",&jitter_,"jitter/F");
  myTree->Branch("ped",&ped_,"ped/F");
  myTree->Branch("calibE",&calibE_,"calibE/F");
  myTree->Branch("E3x3array",E3x3array_,"E3x3array[9]/F");
  myTree->Branch("T3x3array",T3x3array_,"T3x3array[9]/F");
  myTree->Branch("J3x3array",J3x3array_,"J3x3array[9]/F");
  myTree->Branch("lambda3",&lambda3_,"lambda3/F");

  myTree->Branch("chi2Fixed",&chi2Fixed_,"chi2Fixed/F");
  myTree->Branch("chi2Sliding",&chi2Sliding_,"chi2Sliding/F");
  myTree->Branch("time",&time_,"time/F");
  
  myTree->Branch("sampleTime",sampleTime_,"sampleTime[10]/F");
  myTree->Branch("mgpaSample",mgpaSample_,"mgpaSample[10]/F");

  myTree->Branch("R9var",&R9var_,"R9var/F");
  myTree->Branch("gnc",&gnc_,"gnc/I");
  myTree->Branch("nc",&nc_,"nc/I");

  myTree->Branch("ncRelaxed",&ncRelaxed_,"ncRelaxed/I");

  myTree->Branch("kGood",&kGood_,"kGood/I");
  myTree->Branch("kOutOfTime",&kOutOfTime_,"kOutOfTime/I");
  myTree->Branch("kPoorReco",&kPoorReco_,"kPoorReco/I");
  myTree->Branch("kPoorReco_Uncalib",&kPoorReco_Uncalib_,"kPoorReco_Uncalib/I");
  myTree->Branch("kWeird",&kWeird_,"kWeird/I");
  myTree->Branch("kDiWeird",&kDiWeird_,"kDiWeird/I");
  myTree->Branch("ecalSL",&ecalSL_,"ecalSL/I");

  myTree->Branch("calibEFit",&calibEFit_,"calibEFit/F");
  myTree->Branch("rechitPt",&rechitPt_,"rechitPt/F");
  myTree->Branch("rechitEta",&rechitEta_,"rechitEta/F");
  myTree->Branch("rechitPhi",&rechitPhi_,"rechitPhi/F");

  myTree->Branch("rechitFlag",&rechitFlag_,"rechitFlag/i");
  myTree->Branch("isRecovered",&isRecovered_,"isRecovered/I");



  myTree->Branch("timeError",&timeError_,"timeError/F");
  myTree->Branch("sampleGainId",sampleGainId_,"sampleGainId[10]/I");
  myTree->Branch("sampleGainRatio",sampleGainRatio_,"sampleGainRatio[10]/F");

 // ===

  myTree->Branch("ped1",&ped1_,"ped1/F");
  myTree->Branch("ped2",&ped2_,"ped2/F");
  myTree->Branch("ped3",&ped3_,"ped3/F");
  myTree->Branch("pedDB",&pedDB_,"pedDB/F");
  myTree->Branch("pedRMSDB",&pedRMSDB_,"pedRMSDB/F");
  myTree->Branch("nhits",&nhits_,"nhits/I");

  myTree->Branch("timeIC",&timeIC_,"timeIC/F");
  myTree->Branch("gain_switch",&gain_switch_,"gain_switch/I");
  myTree->Branch("iphi",&iphi_,"iphi/I");
  myTree->Branch("ieta",&ieta_,"ieta/I");
  myTree->Branch("iz",&iz_,"iz/I");
  myTree->Branch("runNumber",&runNumber_,"runNumber/I");
  myTree->Branch("lumi",&lumi_,"lumi/I");
  myTree->Branch("eventNumber",&eventNumber_,"eventNumber/I");
  myTree->Branch("subDet",&subDet_,"subDet/I");


  // === store global_data 
  global_data_ = fTFileService_->make<TTree>("global_data","global_data");
  global_data_->Branch("mgpaSimPulseEB",mgpaSimPulseEB_,"mgpaSimPulseEB[250]/F");
  global_data_->Branch("mgpaSimPulseEE",mgpaSimPulseEE_,"mgpaSimPulseEE[250]/F");
  global_data_->Branch("mgpaRatioPulseEB",mgpaRatioPulseEB_,"mgpaRatioPulseEB[250]/F");
  global_data_->Branch("mgpaRatioPulseEE",mgpaRatioPulseEE_,"mgpaRatioPulseEE[250]/F");
  global_data_->Branch("mgpaTime",mgpaTime_,"mgpaTime[250]/F");

/*
  // === store physics 
  physics_ = fTFileService_->make<TTree>("physics","physics");
  physics_->Branch("met",met_,"met[20]/F");
  physics_->Branch("metPhi",metPhi_,"metPhi[20]/F");
  physics_->Branch("metEta",metEta_,"metEta[20]/F");
  physics_->Branch("sumEt",sumEt_,"sumEt[20]/F");
  physics_->Branch("eve",&eve_,"eve/I");
  physics_->Branch("run",&run_,"run/I");
*/
  double ADC_clock = 25.0; // 25 ns

  EBShape mgpaSimPulseShapeEB;
  double risingTimeEB = mgpaSimPulseShapeEB.timeToRise();
  double tzeroEB = risingTimeEB  - 5*ADC_clock;  // 5 samples before the peak

  EEShape mgpaSimPulseShapeEE;
  double risingTimeEE = mgpaSimPulseShapeEE.timeToRise();
  double tzeroEE = risingTimeEE  - 5*ADC_clock;  // 5 samples before the peak

  double amplitudeFitParametersEB[] = {1.138,1.652};
  double amplitudeFitParametersEE[] = {1.890,1.400};

  double alphaEB = amplitudeFitParametersEB[0];
  double betaEB = amplitudeFitParametersEB[1]*ADC_clock;

  double alphaEE = amplitudeFitParametersEE[0];
  double betaEE = amplitudeFitParametersEE[1]*ADC_clock;

  for(int iSample=0;iSample<250;iSample++)
  {
    mgpaSimPulseEB_[iSample] = (mgpaSimPulseShapeEB)(tzeroEB  + float(iSample));
    mgpaSimPulseEE_[iSample] = (mgpaSimPulseShapeEE)(tzeroEE  + float(iSample));
    mgpaTime_[iSample] = float(iSample);

    double ratioPulseEB = 0;
    double termOneEB = 1 + (double(iSample) - 125.0) / (alphaEB * betaEB);
    if (termOneEB > 1.e-5) ratioPulseEB  = exp(alphaEB * log(termOneEB)) * exp(-(iSample - 125.0) / betaEB);

    double ratioPulseEE = 0;
    double termOneEE = 1 + (double(iSample) - 125.0) / (alphaEE * betaEE);
    if (termOneEE > 1.e-5) ratioPulseEE  = exp(alphaEE * log(termOneEE)) * exp(-(iSample - 125.0) / betaEE);

    mgpaRatioPulseEB_[iSample] = ratioPulseEB;
    mgpaRatioPulseEE_[iSample] = ratioPulseEE;
  }


  global_data_->Fill(); // this TTree is filled only once


  // --- define histograms' segmentation
  topoBins = 400;
  topoMin  = 0.0;
  topoMax  = 1.2;

  chronoBins = 400;
  chronoMin  = -15;
  chronoMax  = +15;

  shapeBins = 320;
  shapeMin  = 0;
  shapeMax  = 64;

  energyBins = 8*420;
  energyMin  = -20;
  energyMax  = 400;




 // --- create monitoring histograms for EB
  energy_EB[all] = fTFileService_->make<TH1F>("energy_EB_all","energy_EB_all;energy[GeV];events",energyBins,energyMin,energyMax);
  topo_EB[all] = fTFileService_->make<TH1F>("topo_EB_all","topo_EB_all;E1/E3x3;events",topoBins,topoMin,topoMax);
  chrono_EB[all] = fTFileService_->make<TH1F>("chrono_EB_all","chrono_EB_all;time[ns];events",chronoBins,chronoMin,chronoMax);
  shape_it_EB[all] = fTFileService_->make<TH1F>("shape_it_EB_all","shape_it_EB_all;log[3(#chi2/ndf) + 7];events",shapeBins,shapeMin,shapeMax);
  shape_oot_EB[all] = fTFileService_->make<TH1F>("shape_oot_EB_all","shape_oot_EB_all;log[3(#chi2/ndf) + 7];events",shapeBins,shapeMin,shapeMax);
  shape_itVSenergy_EB[all] = fTFileService_->make<TH2F>("shape_itVSenergy_EB_all","shape_itVSenergy_EB_all",energyBins,energyMin,energyMax,shapeBins,shapeMin,shapeMax);
  shape_itVStopo_EB[all] = fTFileService_->make<TH2F>("shape_itVStopo_EB_all","shape_itVStopo_EB_all",topoBins,topoMin,topoMax,shapeBins,shapeMin,shapeMax);
  shape_itVSchrono_EB[all] = fTFileService_->make<TH2F>("shape_itVSchrono_EB_all","shape_itVSchrono_EB_all",chronoBins,chronoMin,chronoMax,shapeBins,shapeMin,shapeMax);
  shape_ootVSenergy_EB[all] = fTFileService_->make<TH2F>("shape_ootVSenergy_EB_all","shape_ootVSenergy_EB_all",energyBins,energyMin,energyMax,shapeBins,shapeMin,shapeMax);
  shape_ootVStopo_EB[all] = fTFileService_->make<TH2F>("shape_ootVStopo_EB_all","shape_ootVStopo_EB_all",topoBins,topoMin,topoMax,shapeBins,shapeMin,shapeMax);
  shape_ootVSchrono_EB[all] = fTFileService_->make<TH2F>("shape_ootVSchrono_EB_all","shape_ootVSchrono_EB_all",chronoBins,chronoMin,chronoMax,shapeBins,shapeMin,shapeMax);
  topoVSenergy_EB[all] = fTFileService_->make<TH2F>("topoVSenergy_EB_all","topoVSenergy_EB_all",energyBins,energyMin,energyMax,topoBins,topoMin,topoMax);
  chronoVSenergy_EB[all] = fTFileService_->make<TH2F>("chronoVSenergy_EB_all","chronoVSenergy_EB_all",energyBins,energyMin,energyMax,chronoBins,chronoMin,chronoMax);

  energy_EB[cleaned] = fTFileService_->make<TH1F>("energy_EB_cleaned","energy_EB_cleaned;energy[GeV];events",energyBins,energyMin,energyMax);
  topo_EB[cleaned] = fTFileService_->make<TH1F>("topo_EB_cleaned","topo_EB_cleaned;E1/E3x3;events",topoBins,topoMin,topoMax);
  chrono_EB[cleaned] = fTFileService_->make<TH1F>("chrono_EB_cleaned","chrono_EB_cleaned;time[ns];events",chronoBins,chronoMin,chronoMax);
  shape_it_EB[cleaned] = fTFileService_->make<TH1F>("shape_it_EB_cleaned","shape_it_EB_cleaned;log[3(#chi2/ndf) + 7];events",shapeBins,shapeMin,shapeMax);
  shape_oot_EB[cleaned] = fTFileService_->make<TH1F>("shape_oot_EB_cleaned","shape_oot_EB_cleaned;log[3(#chi2/ndf) + 7];events",shapeBins,shapeMin,shapeMax);
  shape_itVSenergy_EB[cleaned] = fTFileService_->make<TH2F>("shape_itVSenergy_EB_cleaned","shape_itVSenergy_EB_cleaned",energyBins,energyMin,energyMax,shapeBins,shapeMin,shapeMax);
  shape_itVStopo_EB[cleaned] = fTFileService_->make<TH2F>("shape_itVStopo_EB_cleaned","shape_itVStopo_EB_cleaned",topoBins,topoMin,topoMax,shapeBins,shapeMin,shapeMax);
  shape_itVSchrono_EB[cleaned] = fTFileService_->make<TH2F>("shape_itVSchrono_EB_cleaned","shape_itVSchrono_EB_cleaned",chronoBins,chronoMin,chronoMax,shapeBins,shapeMin,shapeMax);
  shape_ootVSenergy_EB[cleaned] = fTFileService_->make<TH2F>("shape_ootVSenergy_EB_cleaned","shape_ootVSenergy_EB_cleaned",energyBins,energyMin,energyMax,shapeBins,shapeMin,shapeMax);
  shape_ootVStopo_EB[cleaned] = fTFileService_->make<TH2F>("shape_ootVStopo_EB_cleaned","shape_ootVStopo_EB_cleaned",topoBins,topoMin,topoMax,shapeBins,shapeMin,shapeMax);
  shape_ootVSchrono_EB[cleaned] = fTFileService_->make<TH2F>("shape_ootVSchrono_EB_cleaned","shape_ootVSchrono_EB_cleaned",chronoBins,chronoMin,chronoMax,shapeBins,shapeMin,shapeMax);
  topoVSenergy_EB[cleaned] = fTFileService_->make<TH2F>("topoVSenergy_EB_cleaned","topoVSenergy_EB_cleaned",energyBins,energyMin,energyMax,topoBins,topoMin,topoMax);
  chronoVSenergy_EB[cleaned] = fTFileService_->make<TH2F>("chronoVSenergy_EB_cleaned","chronoVSenergy_EB_cleaned",energyBins,energyMin,energyMax,chronoBins,chronoMin,chronoMax);

  energy_EB[good] = fTFileService_->make<TH1F>("energy_EB_good","energy_EB_good;energy[GeV];events",energyBins,energyMin,energyMax);
  topo_EB[good] = fTFileService_->make<TH1F>("topo_EB_good","topo_EB_good;E1/E3x3;events",topoBins,topoMin,topoMax);
  chrono_EB[good] = fTFileService_->make<TH1F>("chrono_EB_good","chrono_EB_good;time[ns];events",chronoBins,chronoMin,chronoMax);
  shape_it_EB[good] = fTFileService_->make<TH1F>("shape_it_EB_good","shape_it_EB_good;log[3(#chi2/ndf) + 7];events",shapeBins,shapeMin,shapeMax);
  shape_oot_EB[good] = fTFileService_->make<TH1F>("shape_oot_EB_good","shape_oot_EB_good;log[3(#chi2/ndf) + 7];events",shapeBins,shapeMin,shapeMax);
  shape_itVSenergy_EB[good] = fTFileService_->make<TH2F>("shape_itVSenergy_EB_good","shape_itVSenergy_EB_good",energyBins,energyMin,energyMax,shapeBins,shapeMin,shapeMax);
  shape_itVStopo_EB[good] = fTFileService_->make<TH2F>("shape_itVStopo_EB_good","shape_itVStopo_EB_good",topoBins,topoMin,topoMax,shapeBins,shapeMin,shapeMax);
  shape_itVSchrono_EB[good] = fTFileService_->make<TH2F>("shape_itVSchrono_EB_good","shape_itVSchrono_EB_good",chronoBins,chronoMin,chronoMax,shapeBins,shapeMin,shapeMax);
  shape_ootVSenergy_EB[good] = fTFileService_->make<TH2F>("shape_ootVSenergy_EB_good","shape_ootVSenergy_EB_good",energyBins,energyMin,energyMax,shapeBins,shapeMin,shapeMax);
  shape_ootVStopo_EB[good] = fTFileService_->make<TH2F>("shape_ootVStopo_EB_good","shape_ootVStopo_EB_good",topoBins,topoMin,topoMax,shapeBins,shapeMin,shapeMax);
  shape_ootVSchrono_EB[good] = fTFileService_->make<TH2F>("shape_ootVSchrono_EB_good","shape_ootVSchrono_EB_good",chronoBins,chronoMin,chronoMax,shapeBins,shapeMin,shapeMax);
  topoVSenergy_EB[good] = fTFileService_->make<TH2F>("topoVSenergy_EB_good","topoVSenergy_EB_good",energyBins,energyMin,energyMax,topoBins,topoMin,topoMax);
  chronoVSenergy_EB[good] = fTFileService_->make<TH2F>("chronoVSenergy_EB_good","chronoVSenergy_EB_good",energyBins,energyMin,energyMax,chronoBins,chronoMin,chronoMax);


 // --- create monitoring histograms for EE


  energy_EE[all] = fTFileService_->make<TH1F>("energy_EE_all","energy_EE_all;energy[GeV];events",energyBins,energyMin,energyMax);
  topo_EE[all] = fTFileService_->make<TH1F>("topo_EE_all","topo_EE_all;E1/E3x3;events",topoBins,topoMin,topoMax);
  chrono_EE[all] = fTFileService_->make<TH1F>("chrono_EE_all","chrono_EE_all;time[ns];events",chronoBins,chronoMin,chronoMax);
  shape_it_EE[all] = fTFileService_->make<TH1F>("shape_it_EE_all","shape_it_EE_all;log[3(#chi2/ndf) + 7];events",shapeBins,shapeMin,shapeMax);
  shape_oot_EE[all] = fTFileService_->make<TH1F>("shape_oot_EE_all","shape_oot_EE_all;log[3(#chi2/ndf) + 7];events",shapeBins,shapeMin,shapeMax);
  shape_itVSenergy_EE[all] = fTFileService_->make<TH2F>("shape_itVSenergy_EE_all","shape_itVSenergy_EE_all",energyBins,energyMin,energyMax,shapeBins,shapeMin,shapeMax);
  shape_itVStopo_EE[all] = fTFileService_->make<TH2F>("shape_itVStopo_EE_all","shape_itVStopo_EE_all",topoBins,topoMin,topoMax,shapeBins,shapeMin,shapeMax);
  shape_itVSchrono_EE[all] = fTFileService_->make<TH2F>("shape_itVSchrono_EE_all","shape_itVSchrono_EE_all",chronoBins,chronoMin,chronoMax,shapeBins,shapeMin,shapeMax);
  shape_ootVSenergy_EE[all] = fTFileService_->make<TH2F>("shape_ootVSenergy_EE_all","shape_ootVSenergy_EE_all",energyBins,energyMin,energyMax,shapeBins,shapeMin,shapeMax);
  shape_ootVStopo_EE[all] = fTFileService_->make<TH2F>("shape_ootVStopo_EE_all","shape_ootVStopo_EE_all",topoBins,topoMin,topoMax,shapeBins,shapeMin,shapeMax);
  shape_ootVSchrono_EE[all] = fTFileService_->make<TH2F>("shape_ootVSchrono_EE_all","shape_ootVSchrono_EE_all",chronoBins,chronoMin,chronoMax,shapeBins,shapeMin,shapeMax);
  topoVSenergy_EE[all] = fTFileService_->make<TH2F>("topoVSenergy_EE_all","topoVSenergy_EE_all",energyBins,energyMin,energyMax,topoBins,topoMin,topoMax);
  chronoVSenergy_EE[all] = fTFileService_->make<TH2F>("chronoVSenergy_EE_all","chronoVSenergy_EE_all",energyBins,energyMin,energyMax,chronoBins,chronoMin,chronoMax);

  energy_EE[cleaned] = fTFileService_->make<TH1F>("energy_EE_cleaned","energy_EE_cleaned;energy[GeV];events",energyBins,energyMin,energyMax);
  topo_EE[cleaned] = fTFileService_->make<TH1F>("topo_EE_cleaned","topo_EE_cleaned;E1/E3x3;events",topoBins,topoMin,topoMax);
  chrono_EE[cleaned] = fTFileService_->make<TH1F>("chrono_EE_cleaned","chrono_EE_cleaned;time[ns];events",chronoBins,chronoMin,chronoMax);
  shape_it_EE[cleaned] = fTFileService_->make<TH1F>("shape_it_EE_cleaned","shape_it_EE_cleaned;log[3(#chi2/ndf) + 7];events",shapeBins,shapeMin,shapeMax);
  shape_oot_EE[cleaned] = fTFileService_->make<TH1F>("shape_oot_EE_cleaned","shape_oot_EE_cleaned;log[3(#chi2/ndf) + 7];events",shapeBins,shapeMin,shapeMax);
  shape_itVSenergy_EE[cleaned] = fTFileService_->make<TH2F>("shape_itVSenergy_EE_cleaned","shape_itVSenergy_EE_cleaned",energyBins,energyMin,energyMax,shapeBins,shapeMin,shapeMax);
  shape_itVStopo_EE[cleaned] = fTFileService_->make<TH2F>("shape_itVStopo_EE_cleaned","shape_itVStopo_EE_cleaned",topoBins,topoMin,topoMax,shapeBins,shapeMin,shapeMax);
  shape_itVSchrono_EE[cleaned] = fTFileService_->make<TH2F>("shape_itVSchrono_EE_cleaned","shape_itVSchrono_EE_cleaned",chronoBins,chronoMin,chronoMax,shapeBins,shapeMin,shapeMax);
  shape_ootVSenergy_EE[cleaned] = fTFileService_->make<TH2F>("shape_ootVSenergy_EE_cleaned","shape_ootVSenergy_EE_cleaned",energyBins,energyMin,energyMax,shapeBins,shapeMin,shapeMax);
  shape_ootVStopo_EE[cleaned] = fTFileService_->make<TH2F>("shape_ootVStopo_EE_cleaned","shape_ootVStopo_EE_cleaned",topoBins,topoMin,topoMax,shapeBins,shapeMin,shapeMax);
  shape_ootVSchrono_EE[cleaned] = fTFileService_->make<TH2F>("shape_ootVSchrono_EE_cleaned","shape_ootVSchrono_EE_cleaned",chronoBins,chronoMin,chronoMax,shapeBins,shapeMin,shapeMax);
  topoVSenergy_EE[cleaned] = fTFileService_->make<TH2F>("topoVSenergy_EE_cleaned","topoVSenergy_EE_cleaned",energyBins,energyMin,energyMax,topoBins,topoMin,topoMax);
  chronoVSenergy_EE[cleaned] = fTFileService_->make<TH2F>("chronoVSenergy_EE_cleaned","chronoVSenergy_EE_cleaned",energyBins,energyMin,energyMax,chronoBins,chronoMin,chronoMax);

  energy_EE[good] = fTFileService_->make<TH1F>("energy_EE_good","energy_EE_good;energy[GeV];events",energyBins,energyMin,energyMax);
  topo_EE[good] = fTFileService_->make<TH1F>("topo_EE_good","topo_EE_good;E1/E3x3;events",topoBins,topoMin,topoMax);
  chrono_EE[good] = fTFileService_->make<TH1F>("chrono_EE_good","chrono_EE_good;time[ns];events",chronoBins,chronoMin,chronoMax);
  shape_it_EE[good] = fTFileService_->make<TH1F>("shape_it_EE_good","shape_it_EE_good;log[3(#chi2/ndf) + 7];events",shapeBins,shapeMin,shapeMax);
  shape_oot_EE[good] = fTFileService_->make<TH1F>("shape_oot_EE_good","shape_oot_EE_good;log[3(#chi2/ndf) + 7];events",shapeBins,shapeMin,shapeMax);
  shape_itVSenergy_EE[good] = fTFileService_->make<TH2F>("shape_itVSenergy_EE_good","shape_itVSenergy_EE_good",energyBins,energyMin,energyMax,shapeBins,shapeMin,shapeMax);
  shape_itVStopo_EE[good] = fTFileService_->make<TH2F>("shape_itVStopo_EE_good","shape_itVStopo_EE_good",topoBins,topoMin,topoMax,shapeBins,shapeMin,shapeMax);
  shape_itVSchrono_EE[good] = fTFileService_->make<TH2F>("shape_itVSchrono_EE_good","shape_itVSchrono_EE_good",chronoBins,chronoMin,chronoMax,shapeBins,shapeMin,shapeMax);
  shape_ootVSenergy_EE[good] = fTFileService_->make<TH2F>("shape_ootVSenergy_EE_good","shape_ootVSenergy_EE_good",energyBins,energyMin,energyMax,shapeBins,shapeMin,shapeMax);
  shape_ootVStopo_EE[good] = fTFileService_->make<TH2F>("shape_ootVStopo_EE_good","shape_ootVStopo_EE_good",topoBins,topoMin,topoMax,shapeBins,shapeMin,shapeMax);
  shape_ootVSchrono_EE[good] = fTFileService_->make<TH2F>("shape_ootVSchrono_EE_good","shape_ootVSchrono_EE_good",chronoBins,chronoMin,chronoMax,shapeBins,shapeMin,shapeMax);
  topoVSenergy_EE[good] = fTFileService_->make<TH2F>("topoVSenergy_EE_good","topoVSenergy_EE_good",energyBins,energyMin,energyMax,topoBins,topoMin,topoMax);
  chronoVSenergy_EE[good] = fTFileService_->make<TH2F>("chronoVSenergy_EE_good","chronoVSenergy_EE_good",energyBins,energyMin,energyMax,chronoBins,chronoMin,chronoMax);


  for(int i=0;i<3;i++)
  {
     shape_itVSenergy_EB[i]->GetXaxis()->SetTitle("energy[GeV]");
     shape_itVStopo_EB[i]->GetXaxis()->SetTitle("E1/E3x3");
     shape_itVSchrono_EB[i]->GetXaxis()->SetTitle("time[ns]");
     shape_ootVSenergy_EB[i]->GetXaxis()->SetTitle("energy[GeV]");
     shape_ootVStopo_EB[i]->GetXaxis()->SetTitle("E1/E3x3");
     shape_ootVSchrono_EB[i]->GetXaxis()->SetTitle("time[ns]");
     topoVSenergy_EB[i]->GetXaxis()->SetTitle("energy[GeV]");
     chronoVSenergy_EB[i]->GetXaxis()->SetTitle("energy[GeV]");

     shape_itVSenergy_EB[i]->GetYaxis()->SetTitle("log[3(#chi2/ndf) + 7]");
     shape_itVStopo_EB[i]->GetYaxis()->SetTitle("log[3(#chi2/ndf) + 7]");
     shape_itVSchrono_EB[i]->GetYaxis()->SetTitle("log[3(#chi2/ndf) + 7]");
     shape_ootVSenergy_EB[i]->GetYaxis()->SetTitle("log[3(#chi2/ndf) + 7]");
     shape_ootVStopo_EB[i]->GetYaxis()->SetTitle("log[3(#chi2/ndf) + 7]");
     shape_ootVSchrono_EB[i]->GetYaxis()->SetTitle("log[3(#chi2/ndf) + 7]");
     topoVSenergy_EB[i]->GetYaxis()->SetTitle("E1/E3x3");
     chronoVSenergy_EB[i]->GetYaxis()->SetTitle("time[ns]");

  }

}




// ------------ method called once each job just after ending the event loop  ------------
void 
CompareAmplitudes::endJob() 
{
	//delete myTree;
	//delete global_data_;
}

//define this as a plug-in
DEFINE_FWK_MODULE(CompareAmplitudes);
