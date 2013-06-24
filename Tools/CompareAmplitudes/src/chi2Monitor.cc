// -*- C++ -*-
//
// Package:    chi2Monitor
// Class:      chi2Monitor
// 
/**\class chi2Monitor chi2Monitor.cc Tools/CompareAmplitudes/src/chi2Monitor.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Konstantinos Theofilatos
//         Created:  Thu Jul  22 18:38:05 CET 2010
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

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"

#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h"


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




#include <sstream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TProfile2D.h"
#include "TStyle.h"


class chi2Monitor : public edm::EDAnalyzer {
   public:
      explicit chi2Monitor(const edm::ParameterSet&);
      ~chi2Monitor();


      float R9var(const DetId & det, const edm::Handle<EcalRecHitCollection> &hits, const edm::ESHandle<CaloTopology>  &caloTopo);
      float R9PD(const DetId & det, const edm::Handle<EcalRecHitCollection> &hits, const edm::ESHandle<CaloTopology>  &caloTopo);
    



   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      template<class T> std::string any2string(T i);

      enum {EB,EE}; 
      enum {inTime,outOfTime}; 

      edm::Service<TFileService> fTFileService_; 


      edm::InputTag EBURecHitCollection_;
      edm::InputTag EEURecHitCollection_;
      edm::Handle<EcalUncalibratedRecHitCollection> EBURecHits_;
      edm::Handle<EcalUncalibratedRecHitCollection> EEURecHits_;

      edm::InputTag EBRecHitCollection_;
      edm::InputTag EERecHitCollection_;
      edm::Handle<EcalRecHitCollection> EBRecHits_;
      edm::Handle<EcalRecHitCollection> EERecHits_;


      TProfile *p_chi2_amp[3][2][2];
      TProfile *p_chi2_amp_calib[3][2][2];
      TH1F *h_chi2[19][2][2];
      TH1F *h_time[7][2][2];
      TH1F *h_topo[10][2][2];
      TH1F *h_chi2Plain[2][2][2];
      TH1F *h_energy[7][2][2];
      TH1F *h_nenergy[1][2][2];
      TH2F *h_chi2_amp[1][2][2];
      TH2F *h_chi2_top[6][2][2];
      TH2F *h_chi2_jitter[4][2][2];
      TH2F *h_chi2_amp_calib[1][2][2];
      TH2F *h_amp_jitter[1][2][2];
      TH2F *h_chi2_nenergy[1][2][2];
  
      TProfile2D *h_chi2_map[14][2][2];

      TH2F *h_occupancies[17][2][2]; 

      TH1F *h_met[16][2][2];
      TH1F *h_sumEt[16][2][2];
      TH1F *h_metOverSumEt[16][2][2];
  



//----------------------------
//-----------------------------      

	


};

chi2Monitor::chi2Monitor(const edm::ParameterSet& iConfig)
{

    EBURecHitCollection_ = iConfig.getParameter<edm::InputTag>("EBURecHitCollection");
    EEURecHitCollection_ = iConfig.getParameter<edm::InputTag>("EEURecHitCollection");

    EBRecHitCollection_ = iConfig.getParameter<edm::InputTag>("EBRecHitCollection");
    EERecHitCollection_ = iConfig.getParameter<edm::InputTag>("EERecHitCollection");


    //ebEnergyMin_ = iConfig.getUntrackedParameter<double>("ebEnergyMin",0.3);

}


chi2Monitor::~chi2Monitor()
{
}









// ------------ method called to for each event  ------------
void
chi2Monitor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;


//   std::cout << "::: Trying to fetch the collections " << std::endl;

  iEvent.getByLabel(EBURecHitCollection_, EBURecHits_);
  iEvent.getByLabel(EEURecHitCollection_, EEURecHits_);

  iEvent.getByLabel(EBRecHitCollection_, EBRecHits_);
  iEvent.getByLabel(EERecHitCollection_, EERecHits_);



  ESHandle<CaloTopology> caloTopo;
  iSetup.get<CaloTopologyRecord>().get(caloTopo);
  

  EcalUncalibratedRecHitCollection::const_iterator urechitIt;
  EcalRecHitCollection::const_iterator rechitIt;


// Get Geometry
  edm::ESHandle<CaloGeometry> caloGeom;
  iSetup.get<CaloGeometryRecord>().get(caloGeom);
  const EcalBarrelGeometry *ebGeom = (const EcalBarrelGeometry*)caloGeom->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);



  TVector3 metVector[20][2][2];
  float sumEt[20][2][2];
  for(int m=0;m<2;m++)
  {
    for (int n=0;n<2;n++)
    {
      for(int i =0;i<20;i++)
      {
	  metVector[i][m][n].SetXYZ(0.0,0.0,0.0);
	  sumEt[i][m][n] = 0.0;
      }
    }
  }


  for(EBRecHitCollection::const_iterator It = EBRecHits_->begin(); It != EBRecHits_->end(); ++It)
  {

    EBDetId ebDetId  = (*It).id();
    int iphi = ebDetId.iphi();
    int ieta = ebDetId.ieta();
    //int iz = ebDetId.zside();
    //bool gainSwitch = false;
    //int subDet = 0;

    urechitIt = EBURecHits_->find(ebDetId);
    rechitIt = EBRecHits_->find(ebDetId);


    bool urechitFound = false;
    if(urechitIt!=EBURecHits_->end())urechitFound=true;

    bool rechitFound = false;
    if(rechitIt!=EBRecHits_->end())rechitFound=true;




    if(urechitFound && rechitFound)
    {

      //int kFakeNeighbours = 0; if(rechitIt->recoFlag()== EcalRecHit::kFakeNeighbours)kFakeNeighbours=1;
      float r9 = R9var(ebDetId, EBRecHits_, caloTopo);
      float r9pd = R9PD(ebDetId, EBRecHits_, caloTopo);

      int kGood = 0 ; if(rechitIt->recoFlag()== EcalRecHit::kGood)kGood=1;
      int kOutOfTime = 0 ; if(rechitIt->recoFlag()== EcalRecHit::kOutOfTime)kOutOfTime=1;
      //float swissCross = 	 EcalSeverityLevelAlgo::swissCross( It->id(), *EBRecHits_, 0.0 , true);
      float swissCross = -1.0;
//      float swissCrossAfterChi2 = 	 EcalSeverityLevelAlgo::swissCrossAfterChi2( It->id(), *EBRecHits_, 0.0, 33 , true);
      float topology = swissCross;

      // ==== global rechit
      float energyWeights = rechitIt->energy();
      float energyRatio =  rechitIt->outOfTimeEnergy();
      float ampWeights = urechitIt->amplitude();
      float calibConstant = energyWeights/ampWeights;
      float ampRatio = energyRatio/calibConstant;
      float jitter = urechitIt->jitter();
      float chi2it = urechitIt->chi2();
      float time = rechitIt->time();
      float chi2oot = rechitIt->outOfTimeChi2();

      // compute physics 
      const CaloCellGeometry *cell = ebGeom->getGeometry(ebDetId);
      GlobalPoint cellPos = cell->getPosition();
      TVector3 hitPos(cellPos.x(), cellPos.y(), cellPos.z());
      hitPos *= 1.0/hitPos.Mag();
      hitPos *= rechitIt->energy();
      float hitEt = hitPos.Perp();


      // fill for each channel
      for(int chi2type=0;chi2type<2;chi2type++)
      {
	float chi2=0;
	float energy=0;
	float amp=0;
	if(chi2type==0)
	{
		chi2 = chi2it;
		amp = ampWeights;
		energy = energyWeights;
	}
	if(chi2type==1)
	{
		chi2=chi2oot;
		amp = ampRatio;
		energy = energyRatio;
	}

	bool spike = (topology>0.95 && energy>3)? true:false;
//	bool spikeChi2 = (swissCrossAfterChi2>0.95 && energy>3)? true:false;
	bool spikeChi2 = spike; // same as above
	bool outOfTime = kOutOfTime==1 ? true:false;
	bool badShape = chi2 >30 ? true:false;

	// fill histograms
	int cut;

	cut=0; p_chi2_amp[cut][chi2type][EB]->Fill(amp, chi2);
	if(amp>60 && !spike){cut=1; p_chi2_amp[cut][chi2type][EB]->Fill(amp, chi2);}
	if(amp>60 && !(spike||outOfTime)){cut=2; p_chi2_amp[cut][chi2type][EB]->Fill(amp, chi2);}

	if(amp>5){cut=0; p_chi2_amp_calib[cut][chi2type][EB]->Fill(amp, chi2/(amp*amp));}
	if(amp>60 && !spike){cut=1; p_chi2_amp_calib[cut][chi2type][EB]->Fill(amp, chi2/(amp*amp));}
	if(amp>60 && !(spike||outOfTime)){cut=2; p_chi2_amp_calib[cut][chi2type][EB]->Fill(amp, chi2/(amp*amp));}

	cut=0; h_chi2[cut][chi2type][EB]->Fill(chi2);
	if(energy>3){cut=1; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>3 && !spike){cut=2; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>3 && spike){cut=3; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>3 && !(spike||outOfTime)){cut=4; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>3 && (spike||outOfTime)){cut=5; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>0.5){cut=6; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>1){cut=7; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>2){cut=8; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>5){cut=9; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>7){cut=10; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>10){cut=11; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>15){cut=12; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>20){cut=13; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>30){cut=14; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>3 && r9>0.95){cut=15; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>3 && r9<0.95){cut=16; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>3 && r9pd>0.95){cut=17; h_chi2[cut][chi2type][EB]->Fill(chi2);}
	if(energy>3 && r9pd<0.95){cut=18; h_chi2[cut][chi2type][EB]->Fill(chi2);}

	if(true){cut=0; h_chi2_amp[cut][chi2type][EB]->Fill(amp,chi2);}
	if(energy>3){cut=0; h_chi2_top[cut][chi2type][EB]->Fill(topology,chi2);}
	if(energy>3 && !outOfTime){cut=1; h_chi2_top[cut][chi2type][EB]->Fill(topology,chi2);}
	if(energy>3){cut=2; h_chi2_top[cut][chi2type][EB]->Fill(r9,chi2);}
	if(energy>3 && !outOfTime){cut=3; h_chi2_top[cut][chi2type][EB]->Fill(r9,chi2);}
	if(energy>3){cut=4; h_chi2_top[cut][chi2type][EB]->Fill(r9pd,chi2);}
	if(energy>3 && !outOfTime){cut=5; h_chi2_top[cut][chi2type][EB]->Fill(r9pd,chi2);}

	if(amp>20){cut=0; h_chi2_jitter[cut][chi2type][EB]->Fill(jitter*25.0,chi2);}
	if(amp>20 && jitter!=0.0){cut=1; h_chi2_jitter[cut][chi2type][EB]->Fill(jitter*25.0,chi2);}
	if(energy>1 && jitter!=0.0){cut=2; h_chi2_jitter[cut][chi2type][EB]->Fill(jitter*25.0,chi2);}
	if(energy>3 && jitter!=0.0){cut=3; h_chi2_jitter[cut][chi2type][EB]->Fill(jitter*25.0,chi2);}

	if(amp>5){cut=0; h_chi2_amp_calib[cut][chi2type][EB]->Fill(amp,chi2/(amp*amp));}
	if(true){cut=0; h_amp_jitter[cut][chi2type][EB]->Fill(amp,jitter*25.0);}

	if(true){cut=0; h_chi2_nenergy[cut][chi2type][EB]->Fill(energy,chi2);}


	cut=0; h_chi2_map[cut][chi2type][EB]->Fill(iphi,ieta,chi2);
	if(energy>1){cut=1; h_chi2_map[cut][chi2type][EB]->Fill(iphi,ieta,chi2);}
	if(energy>3 && spike){cut=2; h_chi2_map[cut][chi2type][EB]->Fill(iphi,ieta,chi2);}
	if(energy>3 && outOfTime){cut=3; h_chi2_map[cut][chi2type][EB]->Fill(iphi,ieta,chi2);}
	if(energy>3 && badShape){cut=4; h_chi2_map[cut][chi2type][EB]->Fill(iphi,ieta,chi2);}
	if(energy>3 && (spike || outOfTime)){cut=5; h_chi2_map[cut][chi2type][EB]->Fill(iphi,ieta,chi2);}
	if(energy>3 && (spike || badShape)){cut=6; h_chi2_map[cut][chi2type][EB]->Fill(iphi,ieta,chi2);}
	if(energy>3 && (spike || outOfTime || badShape)){cut=7; h_chi2_map[cut][chi2type][EB]->Fill(iphi,ieta,chi2);}
	if(energy>3 && !spike){cut=8; h_chi2_map[cut][chi2type][EB]->Fill(iphi,ieta,chi2);}
	if(energy>3 && !outOfTime){cut=9; h_chi2_map[cut][chi2type][EB]->Fill(iphi,ieta,chi2);}
	if(energy>3 && !badShape){cut=10; h_chi2_map[cut][chi2type][EB]->Fill(iphi,ieta,chi2);}
	if(energy>3 && !(spike || outOfTime)){cut=11; h_chi2_map[cut][chi2type][EB]->Fill(iphi,ieta,chi2);}
	if(energy>3 && !(spike || badShape)){cut=12; h_chi2_map[cut][chi2type][EB]->Fill(iphi,ieta,chi2);}
	if(energy>3 && !(spike || outOfTime || badShape)){cut=13; h_chi2_map[cut][chi2type][EB]->Fill(iphi,ieta,chi2);}

	cut=0; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);
	if(energy>0.07){cut=1; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);}
	if(energy>3){cut=2; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);}
	if(energy>5){cut=3; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);}
	if(energy>20){cut=4; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);}
	if(energy>3 && spike){cut=5; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);}
	if(energy>3 && outOfTime){cut=6; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);}
	if(energy>3 && badShape){cut=7; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);}
	if(energy>3 && (spike || outOfTime)){cut=8; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);}
	if(energy>3 && (spike || badShape)){cut=9; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);}
	if(energy>3 && (spike || outOfTime || badShape)){cut=10; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);}
	if(energy>3 && !spike){cut=11; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);}
	if(energy>3 && !outOfTime){cut=12; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);}
	if(energy>3 && !badShape){cut=13; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);}
	if(energy>3 && !(spike || outOfTime)){cut=14; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);}
	if(energy>3 && !(spike || badShape)){cut=15; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);}
	if(energy>3 && !(spike || outOfTime || badShape)){cut=16; h_occupancies[cut][chi2type][EB]->Fill(iphi,ieta);}

	cut=0; h_energy[cut][chi2type][EB]->Fill(energy);
	if(!spike){cut=1; h_energy[cut][chi2type][EB]->Fill(energy);}
	if(!outOfTime){cut=2; h_energy[cut][chi2type][EB]->Fill(energy);}
	if(!badShape){cut=3; h_energy[cut][chi2type][EB]->Fill(energy);}
	if(!(spike || outOfTime)){cut=4; h_energy[cut][chi2type][EB]->Fill(energy);}
	if(!(spike || badShape)){cut=5; h_energy[cut][chi2type][EB]->Fill(energy);}
	if(!(spike || outOfTime || badShape)){cut=6; h_energy[cut][chi2type][EB]->Fill(energy);}

	cut=0; h_nenergy[cut][chi2type][EB]->Fill(energy);

	{cut=0; h_chi2Plain[cut][chi2type][EB]->Fill(chi2);}
	{cut=1; h_chi2Plain[cut][chi2type][EB]->Fill(exp(chi2/7 - 3));}


	if(energy>1){cut=0; h_time[cut][chi2type][EB]->Fill(time);}
	if(energy>2){cut=1; h_time[cut][chi2type][EB]->Fill(time);}
	if(energy>3){cut=2; h_time[cut][chi2type][EB]->Fill(time);}
	if(energy>5){cut=3; h_time[cut][chi2type][EB]->Fill(time);}
	if(energy>10){cut=4; h_time[cut][chi2type][EB]->Fill(time);}
	if(energy>15){cut=5; h_time[cut][chi2type][EB]->Fill(time);}
	if(energy>30){cut=6; h_time[cut][chi2type][EB]->Fill(time);}

	if(energy>1){cut=0; h_topo[cut][chi2type][EB]->Fill(topology);}
	if(energy>2){cut=1; h_topo[cut][chi2type][EB]->Fill(topology);}
	if(energy>3){cut=2; h_topo[cut][chi2type][EB]->Fill(topology);}
	if(energy>5){cut=3; h_topo[cut][chi2type][EB]->Fill(topology);}
	if(energy>10){cut=4; h_topo[cut][chi2type][EB]->Fill(topology);}
	if(energy>15){cut=5; h_topo[cut][chi2type][EB]->Fill(topology);}
	if(energy>30){cut=6; h_topo[cut][chi2type][EB]->Fill(topology);}
	if(energy>3 && !outOfTime){cut=7; h_topo[cut][chi2type][EB]->Fill(topology);}
	if(energy>3 && !badShape){cut=8; h_topo[cut][chi2type][EB]->Fill(topology);}
	if(energy>3 && !(badShape || outOfTime)){cut=9; h_topo[cut][chi2type][EB]->Fill(topology);}


  	if(energy>0.07){cut=0;metVector[cut][chi2type][EB]+=hitPos; sumEt[cut][chi2type][EB]+=hitEt;} 
    	if(energy>0.07 && !spike){cut=1;metVector[cut][chi2type][EB]+=hitPos; sumEt[cut][chi2type][EB]+=hitEt;} 
      	if(energy>0.07 && !outOfTime){cut=2;metVector[cut][chi2type][EB]+=hitPos; sumEt[cut][chi2type][EB]+=hitEt;} 
	if(energy>0.07 && !badShape){cut=3;metVector[cut][chi2type][EB]+=hitPos; sumEt[cut][chi2type][EB]+=hitEt;} 
	if(energy>0.07 && !(spike || outOfTime)){cut=4;metVector[cut][chi2type][EB]+=hitPos; sumEt[cut][chi2type][EB]+=hitEt;} 
	if(energy>0.07 && !(spike || badShape)){cut=5;metVector[cut][chi2type][EB]+=hitPos; sumEt[cut][chi2type][EB]+=hitEt;} 
	if(energy>0.07 && !(spike || badShape || outOfTime)){cut=6;metVector[cut][chi2type][EB]+=hitPos; sumEt[cut][chi2type][EB]+=hitEt;} 
	if(energy>0.07 && !(spikeChi2)){cut=7;metVector[cut][chi2type][EB]+=hitPos; sumEt[cut][chi2type][EB]+=hitEt;} 
	if(energy>0.07 && !(spikeChi2 || outOfTime)){cut=8;metVector[cut][chi2type][EB]+=hitPos; sumEt[cut][chi2type][EB]+=hitEt;} 
	if(energy>0.07 && !(spikeChi2 || badShape)){cut=9;metVector[cut][chi2type][EB]+=hitPos; sumEt[cut][chi2type][EB]+=hitEt;} 
	if(energy>0.07 && !(spikeChi2 || outOfTime || badShape)){cut=10;metVector[cut][chi2type][EB]+=hitPos; sumEt[cut][chi2type][EB]+=hitEt;} 
	if(energy>0.07 && !(chi2>28)){cut=11;metVector[cut][chi2type][EB]+=hitPos; sumEt[cut][chi2type][EB]+=hitEt;} 
	if(energy>0.07 && !(chi2>25)){cut=12;metVector[cut][chi2type][EB]+=hitPos; sumEt[cut][chi2type][EB]+=hitEt;} 
	if(energy>0.07 && !(chi2>33)){cut=13;metVector[cut][chi2type][EB]+=hitPos; sumEt[cut][chi2type][EB]+=hitEt;} 
	if(energy>0.07 && !(chi2>30 || outOfTime)){cut=14;metVector[cut][chi2type][EB]+=hitPos; sumEt[cut][chi2type][EB]+=hitEt;} 
	if(energy>0.07 && !(chi2>38 || outOfTime)){cut=15;metVector[cut][chi2type][EB]+=hitPos; sumEt[cut][chi2type][EB]+=hitEt;} 

      }  // 2 chi2 loop

    } // if rechit found
  } // end of rechit loop 

  // ======
  // Fill event histograms
  for(int cut=0;cut<16;cut++)
  {
    for(int chi2type=0;chi2type<2;chi2type++)
    {
      h_met[cut][chi2type][EB]->Fill( metVector[cut][chi2type][EB].Perp() );
      h_sumEt[cut][chi2type][EB]->Fill( sumEt[cut][chi2type][EB] );
      h_metOverSumEt[cut][chi2type][EB]->Fill( metVector[cut][chi2type][EB].Perp()/sumEt[cut][chi2type][EB] );
    }
  }

/*
   for(EcalRecHitCollection::const_iterator It = EERecHits_->begin(); It != EERecHits_->end(); ++It)
   {
      EEDetId eeDetId  = (*It).id();
      urechitIt = EEURecHits_->find(eeDetId);
    
   } 
*/

}

template<class T>
std::string chi2Monitor::any2string(T i)
{
        std::ostringstream buffer;
        buffer << i;
        return buffer.str();
}



// ------------ method called once each job just before starting event loop  ------------
void chi2Monitor::beginJob()
{
  std::string title="";

  for(int i=0;i<2;i++)
  for(int j=0;j<1;j++) // do not create EE histograms
  {
    std::string subDetName = "_EB"; 
    if(j==1)subDetName = "_EE";
    std::string chi2Type = "_it";
    if(i==1)chi2Type = "_oot"; 

    std::string suffix = subDetName+chi2Type;

    for (int cut=0;cut<1;cut++)
    {
      title = "h_chi2_amp"+suffix+"_c"+any2string(cut);
      h_chi2_amp[cut][i][j] = fTFileService_->make<TH2F>(title.c_str(),title.c_str(),2020,-20.,2000.,1000,0.,101.);
      h_chi2_amp[cut][i][j]->GetXaxis()->SetTitle("amplitude (ADC)");
      h_chi2_amp[cut][i][j]->GetYaxis()->SetTitle("chi2");
    }

    for (int cut=0;cut<3;cut++)
    {
      title = "p_chi2_amp"+suffix+"_c"+any2string(cut);
      p_chi2_amp[cut][i][j] = fTFileService_->make<TProfile>(title.c_str(),title.c_str(),5000,5., 2000.);
      p_chi2_amp[cut][i][j]->GetXaxis()->SetTitle("amplitude (ADC)");
      p_chi2_amp[cut][i][j]->GetYaxis()->SetTitle("chi2");
    }

    for (int cut=0;cut<3;cut++)
    {
      title = "p_chi2_amp_calib"+suffix+"_c"+any2string(cut);
      p_chi2_amp_calib[cut][i][j] = fTFileService_->make<TProfile>(title.c_str(),title.c_str(),5000,5., 2000.);
      p_chi2_amp_calib[cut][i][j]->GetXaxis()->SetTitle("amplitude (ADC)");
      p_chi2_amp_calib[cut][i][j]->GetYaxis()->SetTitle("#Sigma Ri^{2}/A^{2}");
    }
    for (int cut=0;cut<19;cut++)
    {
      title = "h_chi2"+suffix+"_c"+any2string(cut);
      h_chi2[cut][i][j] = fTFileService_->make<TH1F>(title.c_str(),title.c_str(),1000,0.,101.);
      h_chi2[cut][i][j]->GetXaxis()->SetTitle("chi2");
      h_chi2[cut][i][j]->GetYaxis()->SetTitle("events");
    }


    for (int cut=0;cut<1;cut++)
    {
      title = "h_chi2_nenergy"+suffix+"_c"+any2string(cut);
      h_chi2_nenergy[cut][i][j] = fTFileService_->make<TH2F>(title.c_str(),title.c_str(),5000,-20.,20.,1000,0.,101.);
      h_chi2_nenergy[cut][i][j]->GetXaxis()->SetTitle("energy (GeV)");
      h_chi2_nenergy[cut][i][j]->GetYaxis()->SetTitle("chi2");
    }


    for (int cut=0;cut<1;cut++)
    {
      title = "h_chi2_amp_calib"+suffix+"_c"+any2string(cut);
      h_chi2_amp_calib[cut][i][j] = fTFileService_->make<TH2F>(title.c_str(),title.c_str(),1000,5.,2000.,3000,0,0.1);
      h_chi2_amp_calib[cut][i][j]->GetXaxis()->SetTitle("amplitude (ADC)");
      h_chi2_amp_calib[cut][i][j]->GetYaxis()->SetTitle("#Sigma Ri^{2}/A^{2}");
    }

    for (int cut=0;cut<6;cut++)
    {
      title = "h_chi2_top"+suffix+"_c"+any2string(cut);
      h_chi2_top[cut][i][j] = fTFileService_->make<TH2F>(title.c_str(),title.c_str(),1000,0.,1.2,1000,0.,101.);
      h_chi2_top[cut][i][j]->GetXaxis()->SetTitle("1-S4/S1");
      h_chi2_top[cut][i][j]->GetYaxis()->SetTitle("chi2");
    }

    for (int cut=0;cut<4;cut++)
    {
      title = "h_chi2_jitter"+suffix+"_c"+any2string(cut);
      h_chi2_jitter[cut][i][j] = fTFileService_->make<TH2F>(title.c_str(),title.c_str(),1000,-100.,100.,1000,0.,101.);
      h_chi2_jitter[cut][i][j]->GetXaxis()->SetTitle("jitter*25 (ns)");
      h_chi2_jitter[cut][i][j]->GetYaxis()->SetTitle("chi2");
    }

    for (int cut=0;cut<1;cut++)
    {
      title = "h_amp_jitter"+suffix+"_c"+any2string(cut);
      h_amp_jitter[cut][i][j] = fTFileService_->make<TH2F>(title.c_str(),title.c_str(),2000,5.,2000,1000,-100.,100.);
      h_amp_jitter[cut][i][j]->GetXaxis()->SetTitle("jitter*25 (ns)");
      h_amp_jitter[cut][i][j]->GetYaxis()->SetTitle("chi2");
    }

    for (int cut=0;cut<14;cut++)
    {
      title = "h_chi2_map"+suffix+"_c"+any2string(cut);
      h_chi2_map[cut][i][j] = fTFileService_->make<TProfile2D>(title.c_str(),title.c_str(),361 , 0.0, 361.0, 172, -86, 86);
      h_chi2_map[cut][i][j]->GetXaxis()->SetTitle("phi");
      h_chi2_map[cut][i][j]->GetYaxis()->SetTitle("eta");
    }

    for (int cut=0;cut<17;cut++)
    {
      title = "h_occupancies"+suffix+"_c"+any2string(cut);
      h_occupancies[cut][i][j] = fTFileService_->make<TH2F>(title.c_str(),title.c_str(),361 , 0.0, 361.0, 172, -86, 86);
      h_occupancies[cut][i][j]->GetXaxis()->SetTitle("phi");
      h_occupancies[cut][i][j]->GetYaxis()->SetTitle("eta");
    }

    for (int cut=0;cut<7;cut++)
    {
      title = "h_energy"+suffix+"_c"+any2string(cut);
      h_energy[cut][i][j] = fTFileService_->make<TH1F>(title.c_str(),title.c_str(),1800, -10.,800.);
      h_energy[cut][i][j]->GetXaxis()->SetTitle("energy (GeV)");
      h_energy[cut][i][j]->GetYaxis()->SetTitle("events");
    }

    for (int cut=0;cut<1;cut++)
    {
      title = "h_nenergy"+suffix+"_c"+any2string(cut);
      h_nenergy[cut][i][j] = fTFileService_->make<TH1F>(title.c_str(),title.c_str(),4000, -20.,20.);
      h_nenergy[cut][i][j]->GetXaxis()->SetTitle("energy (GeV)");
      h_nenergy[cut][i][j]->GetYaxis()->SetTitle("events");
    }

    for (int cut=0;cut<2;cut++)
    {
      title = "h_chi2Plain"+suffix+"_c"+any2string(cut);
      h_chi2Plain[cut][i][j] = fTFileService_->make<TH1F>(title.c_str(),title.c_str(),1000,0.,65.);
      h_chi2Plain[cut][i][j]->GetXaxis()->SetTitle("chi2");
      h_chi2Plain[cut][i][j]->GetYaxis()->SetTitle("events");
    }

    for (int cut=0;cut<16;cut++)
    {
      title = "h_met"+suffix+"_c"+any2string(cut);
      h_met[cut][i][j] = fTFileService_->make<TH1F>(title.c_str(),title.c_str(),1000,0.,500.);
      h_met[cut][i][j]->GetXaxis()->SetTitle("ecal met");
      h_met[cut][i][j]->GetYaxis()->SetTitle("events");
    }

    for (int cut=0;cut<16;cut++)
    {
      title = "h_sumEt"+suffix+"_c"+any2string(cut);
      h_sumEt[cut][i][j] = fTFileService_->make<TH1F>(title.c_str(),title.c_str(),1000,0.,500.);
      h_sumEt[cut][i][j]->GetXaxis()->SetTitle("ecal sumEt");
      h_sumEt[cut][i][j]->GetYaxis()->SetTitle("events");
    }

    for (int cut=0;cut<16;cut++)
    {
      title = "h_metOverSumEt"+suffix+"_c"+any2string(cut);
      h_metOverSumEt[cut][i][j] = fTFileService_->make<TH1F>(title.c_str(),title.c_str(),1000,0.,1.1);
      h_metOverSumEt[cut][i][j]->GetXaxis()->SetTitle("ecal met/sumEt");
      h_metOverSumEt[cut][i][j]->GetYaxis()->SetTitle("events");
    }

    for (int cut=0;cut<7;cut++)
    {
      title = "h_time"+suffix+"_c"+any2string(cut);
      h_time[cut][i][j] = fTFileService_->make<TH1F>(title.c_str(),title.c_str(),1000,-75.,75.);
      h_time[cut][i][j]->GetXaxis()->SetTitle("time (ns)");
      h_time[cut][i][j]->GetYaxis()->SetTitle("events");
    }


    for (int cut=0;cut<10;cut++)
    {
      title = "h_topo"+suffix+"_c"+any2string(cut);
      h_topo[cut][i][j] = fTFileService_->make<TH1F>(title.c_str(),title.c_str(),1000,0.,1.2);
      h_topo[cut][i][j]->GetXaxis()->SetTitle("1 - S4/S1");
      h_topo[cut][i][j]->GetYaxis()->SetTitle("events");
    }

  }


}




// ------------ method called once each job just after ending the event loop  ------------
void 
chi2Monitor::endJob() 
{
}

float chi2Monitor::R9var(const DetId & det, const edm::Handle<EcalRecHitCollection> &hits, const edm::ESHandle<CaloTopology>  &caloTopo)
{
  
  float R9var = -1;
  if(hits->find(det) !=hits->end())
  {
    
    CaloNavigator<DetId> cursor = CaloNavigator<DetId>(det,caloTopo->getSubdetectorTopology(det));
    int side_=3;
    
    EcalRecHit centralHit = *hits->find(det);

    float E1 =  0;
    float E3x3 = 0;
    
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
	    E3x3 += tmpHit.energy();
        } 
      } 
    } 
    if(E3x3!=0)R9var = E1/E3x3;
    
  }

  return R9var;
}

float chi2Monitor::R9PD(const DetId & det, const edm::Handle<EcalRecHitCollection> &hits, const edm::ESHandle<CaloTopology>  &caloTopo)
{
  
  float R9PD = -1;
  if(hits->find(det) !=hits->end())
  {
    
    CaloNavigator<DetId> cursor = CaloNavigator<DetId>(det,caloTopo->getSubdetectorTopology(det));
    int side_=3;
    
    EcalRecHit centralHit = *hits->find(det);

    float E1 =  0;
    float E3x3 = 0;
    
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
            float tmpE = tmpHit.energy();
            float tmpChi2 = tmpHit.chi2();

	    if(!(i==0 && j==0)) // do not apply chi2 in the central crystal
    	    {
    		if(tmpE>0 && tmpChi2<33) E3x3 += tmpHit.energy();
	    }
	    else{if(tmpE>0)E3x3 += tmpHit.energy();} // for the central crystal
        } 
      } 
    } 
    if(E3x3!=0)R9PD = E1/E3x3;
    
  }

  return R9PD;
}


//define this as a plug-in
DEFINE_FWK_MODULE(chi2Monitor);
