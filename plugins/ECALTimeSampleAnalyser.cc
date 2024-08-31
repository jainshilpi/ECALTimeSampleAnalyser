// -*- C++ -*-
//
// Package:    Analyser/ECALTimeSampleAnalyser
// Class:      ECALTimeSampleAnalyser
//
/**\class ECALTimeSampleAnalyser ECALTimeSampleAnalyser.cc Analyser/ECALTimeSampleAnalyser/plugins/ECALTimeSampleAnalyser.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Shilpi Jain
//         Created:  Fri, 12 May 2023 19:29:28 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRawData/interface/EcalRawDataCollections.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/EcalObjects/interface/EcalPFRecHitThresholds.h"
#include "CondFormats/DataRecord/interface/EcalPFRecHitThresholdsRcd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include <TTree.h>
#include <TLorentzVector.h>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class ECALTimeSampleAnalyser : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ECALTimeSampleAnalyser(const edm::ParameterSet&);
      ~ECALTimeSampleAnalyser();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      static const int NMAXSAMPLES = 10;
      static const int NMAXCRYS = 100;
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  //virtual std::vector<std::array<double, NMAXSAMPLES> > getTimeSamplesAroundEle(std::vector<DetId> v_id, edm::Handle<EBDigiCollection> pEBDigi, edm::Handle<EEDigiCollection> pEEDigi);
  virtual std::vector<std::vector<double> > getTimeSamplesAroundEle(const CaloGeometry* geo, std::vector<DetId> v_id, edm::Handle<EBDigiCollection> pEBDigi, edm::Handle<EEDigiCollection> pEEDigi, const EcalRecHitCollection* EBRecHits, const EcalRecHitCollection* EERecHits,  const EcalPFRecHitThresholds* thresholds, std::vector<double> &hitsEnergy, std::vector<double> &hitsThr, std::vector<double> &hitsEta, std::vector<double> &hitsPhi);

  virtual std::vector<reco::GenParticle>::const_iterator  getGenMatch(std::vector<std::vector<reco::GenParticle>::const_iterator> genLep, reco::GsfElectron gsfele, double &dRmin);
  

  virtual bool getSimHitMatch(DetId id, edm::Handle<edm::PCaloHitContainer> pEBSim, edm::Handle<edm::PCaloHitContainer> pEESim, double &matchSimEnergyEB, double &matchSimEnergyEE);
  
      // ----------member data ---------------------------
      edm::ESGetToken<CaloTopology, CaloTopologyRecord> caloTopoToken_;
      edm::ESHandle<CaloTopology> caloTopology_;
      edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<EBDigiCollection> ebDigiToken_;
      edm::EDGetTokenT<EEDigiCollection> eeDigiToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken_;

      const std::string ebdigiCollection_;
  //edm::EDGetTokenT<reco::PhotonCollection>         recophotonCollection_;
      edm::EDGetTokenT<reco::GsfElectronCollection>         recoelectronCollection_;
  
      edm::EDGetTokenT<EcalRecHitCollection>           ebRecHitCollection_;
      edm::EDGetTokenT<EcalRecHitCollection>           eeRecHitCollection_;
  //edm::EDGetTokenT<EcalRecHitCollection>           esRecHitCollection_;
  
      edm::ESHandle<EcalPedestals> peds;
      edm::ESGetToken<EcalPedestals, EcalPedestalsRcd> pedsToken_;
      edm::ESHandle<EcalGainRatios> gains;
      edm::ESGetToken<EcalGainRatios, EcalGainRatiosRcd> gainsToken_;

      edm::EDGetTokenT<edm::PCaloHitContainer>           simHitEBCollection_;
      edm::EDGetTokenT<edm::PCaloHitContainer>           simHitEECollection_;
      
  //const std::string ebdigiProducer_;

  bool runMinBias_;
  bool debugL1_ = true;
  
  TTree   *treeEle;
  TTree   *treePho;
  TTree   *treeNoiseEB;
  TTree   *treeNoiseEE;

  double eleEta_, elePhi_, elePt_, eleE_;
  double e5x5_; 
  int nCrys_;
  int nsamples_;
  double genPt_, genEta_, genPhi_, gendR_, genE_, genStatus_;
  //std::vector<std::array<double, NMAXSAMPLES> > hitsAmplitudes_;
  //std::vector<std::array<double, 10> > hitsAmplitudes_;
  //double hitsAmplitudes_[NMAXCRYS][NMAXSAMPLES];
  std::vector<std::vector<double>> hitsAmplitudes_;
  std::vector<double> hitsEnergy_;
  std::vector<double> hitsThr_;
  std::vector<double> hitsEta_;
  std::vector<double> hitsPhi_;
  std::vector<int> hitsIsSimMatch_;
  std::vector<double> hitsSimEnEB_;
  std::vector<double> hitsSimEnEE_;


  ///hits noise
  std::vector<std::vector<double>> hitsNoiseEBAmplitudes_;
  std::vector<double> hitsNoiseEBEnergy_;
  std::vector<double> hitsNoiseEBThr_;
  std::vector<double> hitsNoiseEBEta_;
  std::vector<double> hitsNoiseEBPhi_;

  std::vector<std::vector<double>> hitsNoiseEEAmplitudes_;
  std::vector<double> hitsNoiseEEEnergy_;
  std::vector<double> hitsNoiseEEThr_;
  std::vector<double> hitsNoiseEEEta_;
  std::vector<double> hitsNoiseEEPhi_;


  edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticlesCollection_;



};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ECALTimeSampleAnalyser::ECALTimeSampleAnalyser(const edm::ParameterSet& iConfig):
      geometryToken_(esConsumes())
{
   //now do what ever initialization is needed
  ebDigiToken_ = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("EBdigiCollection"));
  eeDigiToken_ = consumes<EEDigiCollection>(iConfig.getParameter<edm::InputTag>("EEdigiCollection"));
  //recophotonCollection_       = consumes<reco::PhotonCollection>        (iConfig.getParameter<edm::InputTag>("recoPhotonSrc"));
  recoelectronCollection_       = consumes<reco::GsfElectronCollection>        (iConfig.getParameter<edm::InputTag>("recoEleSrc"));
  ebRecHitCollection_ = consumes<EcalRecHitCollection>          (iConfig.getParameter<edm::InputTag>("ebRecHitCollection"));
  eeRecHitCollection_ = consumes<EcalRecHitCollection>          (iConfig.getParameter<edm::InputTag>("eeRecHitCollection"));
  //esRecHitCollection_ = consumes<EcalRecHitCollection>          (iConfig.getParameter<edm::InputTag>("esRecHitCollection"));
  

  caloTopoToken_ = esConsumes();
  pedsToken_ = esConsumes<EcalPedestals, EcalPedestalsRcd>();
  gainsToken_ = esConsumes<EcalGainRatios, EcalGainRatiosRcd>();
  
  genParticlesCollection_   = consumes<std::vector<reco::GenParticle> >    (iConfig.getParameter<edm::InputTag>("genParticleSrc"));

  simHitEBCollection_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("simHitEBCollection"));
  simHitEECollection_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("simHitEECollection"));

  runMinBias_         = iConfig.getParameter<bool>("runMinBias");
}


ECALTimeSampleAnalyser::~ECALTimeSampleAnalyser()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ECALTimeSampleAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ///https://cmssdt.cern.ch/lxr/source/CalibCalorimetry/EcalLaserAnalyzer/plugins/EcalTestPulseAnalyzer.cc
  /// https://cmssdt.cern.ch/lxr/source/RecoEcal/EgammaClusterProducers/src/EcalDigiSelector.cc ---> selected digis
  ////seems to save 3x3 around the max hit only: https://cmssdt.cern.ch/lxr/source/RecoEcal/EgammaClusterProducers/src/EcalDigiSelector.cc#0164
  ///pedestal subtraction ---> 
  ///https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/EcalRecAlgos/src/EcalUncalibRecHitMultiFitAlgo.cc#0064

   using namespace edm;

   eleEta_ = -99;
   elePhi_ = -99;
   elePt_ = -99.;
   eleE_ = -99;
   nCrys_ = -99;
   e5x5_ = -99;
   hitsAmplitudes_.clear();
   hitsEnergy_.clear();
   hitsThr_.clear();
   hitsEta_.clear();
   hitsPhi_.clear();
   hitsIsSimMatch_.clear();
   hitsSimEnEB_.clear();
   hitsSimEnEE_.clear();

   hitsNoiseEBAmplitudes_.clear();
   hitsNoiseEBEnergy_.clear();
   hitsNoiseEBThr_.clear();
   hitsNoiseEBEta_.clear();
   hitsNoiseEBPhi_.clear();

   hitsNoiseEEAmplitudes_.clear();
   hitsNoiseEEEnergy_.clear();
   hitsNoiseEEThr_.clear();
   hitsNoiseEEEta_.clear();
   hitsNoiseEEPhi_.clear();

   
  // get geometry
  const CaloGeometry* geo = &iSetup.getData(geometryToken_);

   //hitsAmplitudes_.clear();
   nsamples_ = NMAXSAMPLES;

   //noZS::EcalClusterLazyTools lazyToolnoZS(iEvent, iSetup, ebRecHitCollection_, eeRecHitCollection_, esRecHitCollection_);
   
   caloTopology_ = iSetup.getHandle(caloTopoToken_);
   gains = iSetup.getHandle(gainsToken_);
   peds = iSetup.getHandle(pedsToken_);
   
   edm::Handle<EcalRecHitCollection> barrelRecHitsHandle;
   edm::Handle<EcalRecHitCollection> endcapRecHitsHandle;
   //edm::Handle<EcalRecHitCollection> esRecHitsHandle;
   iEvent.getByToken(ebRecHitCollection_,barrelRecHitsHandle);
   iEvent.getByToken(eeRecHitCollection_,endcapRecHitsHandle);
   //iEvent.getByToken(esRecHitCollection_,esRecHitsHandle);

   const EcalRecHitCollection* EBRecHits = nullptr;
   const EcalRecHitCollection* EERecHits = nullptr;
   //const EcalRecHitCollection* ESRecHits = nullptr;
  
   if ( !barrelRecHitsHandle.isValid() ){
     LogDebug("") << "Error! EB rechits can't get product!" << std::endl;
   } else{
     EBRecHits = barrelRecHitsHandle.product();
   }

   if ( !endcapRecHitsHandle.isValid() ){
     LogDebug("") << "Error! EE rechits can't get product!" << std::endl;
   } else{
     EERecHits = endcapRecHitsHandle.product();
   }
   
   /*if ( !esRecHitsHandle.isValid() ){
     LogDebug("") << "Error! ES rechits can't get product!" << std::endl;
   } else{
     ESRecHits = esRecHitsHandle.product();
   }
   */

   // retrieving crystal data from Event
   edm::Handle<EBDigiCollection> pEBDigi;
   edm::Handle<EEDigiCollection> pEEDigi;

   
   iEvent.getByToken(ebDigiToken_, pEBDigi);
   iEvent.getByToken(eeDigiToken_, pEEDigi);
   
   if (!pEBDigi.isValid()) {
     std::cout<<"Error! can't get the product retrieving EB crystal data, i.e. EBDigiCollection " <<std::endl;
   } 
   
   if (!pEEDigi.isValid()) {
     std::cout<<"Error! can't get the product retrieving EE crystal data, i.e. EEDigiCollection " <<std::endl;
   }


   edm::Handle<edm::PCaloHitContainer> pEBSim;
   edm::Handle<edm::PCaloHitContainer> pEESim;
   
   iEvent.getByToken(simHitEBCollection_, pEBSim);
   iEvent.getByToken(simHitEECollection_, pEESim);

   if (!pEBSim.isValid()) {
     std::cout<<"Error! can't get the product retrieving EB crystal data, i.e. PCaloHitContainer " <<std::endl;
   } 
   
   if (!pEESim.isValid()) {
     std::cout<<"Error! can't get the product retrieving EE crystal data, i.e. PCaloHitContainer " <<std::endl;
   }


   //https://github.com/swagata87/OldLocalCovMiniAOD/blob/main/plugins/OldLocalCovMiniAOD.cc
   edm::ESHandle<EcalPFRecHitThresholds> pThresholds;
   iSetup.get<EcalPFRecHitThresholdsRcd>().get(pThresholds);
   const EcalPFRecHitThresholds* thresholds = pThresholds.product();


   /// gen collection
   edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle;
   iEvent.getByToken(genParticlesCollection_, genParticlesHandle);

   if (!genParticlesHandle.isValid()) {
     edm::LogWarning("ggNtuplizer") << "no reco::GenParticles in event";
     return;
   }
   
     std::vector<std::vector<reco::GenParticle>::const_iterator> genLep;

     
     for (std::vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
     

     

     Int_t status = ip->status();
     bool photonOrLepton =
       (ip->pdgId() == 22 && (ip->isPromptFinalState() || ip->isLastCopy() || status == 1)) ||
       (status == 1 && abs(ip->pdgId()) == 11 && (ip->isPromptFinalState() || ip->isLastCopy())) ||
       (status == 1 && abs(ip->pdgId()) == 13 && (ip->isPromptFinalState() || ip->isLastCopy())) ||
       (status == 1 && (abs(ip->pdgId()) == 12 || abs(ip->pdgId()) == 14 || abs(ip->pdgId()) == 16)) ||
       (status == 1 && ( abs(ip->pdgId()) >= 11 && abs(ip->pdgId()) <= 16 ) && ip->pt() > 3.0)  ||
       (status < 10 && abs(ip->pdgId()) == 15 && ip->pt() > 3.0);

     bool isEle = ( abs(ip->pdgId()) == 11 && (ip->isPromptFinalState() || ip->isLastCopy()) );
     
     if(photonOrLepton){
       /*mcPt_     .push_back(ip->pt());
       mcMass_   .push_back(ip->mass());
       mcEta_    .push_back(ip->eta());
       mcPhi_    .push_back(ip->phi());
       mcE_      .push_back(ip->energy());
       mcEt_     .push_back(ip->et());
       mcStatus_ .push_back(ip->status());
       */
       genLep.push_back(ip);
     }
     
   }

     ///to get the noisy crystal collection, go in opp. eta and same phi
     std::vector<DetId> seedIDvec;
     std::vector<double> eta;
     std::vector<double> phi;

     if(!runMinBias_){
       ///grab electron collection
       Handle<reco::GsfElectronCollection> theRecoEleCollection;
       iEvent.getByToken(recoelectronCollection_, theRecoEleCollection);
       const reco::GsfElectronCollection theRecoEl = *(theRecoEleCollection.product());
       
       if (theRecoEleCollection.isValid()) {
	 for (uint j = 0; j < theRecoEl.size(); j++) {
	   
	   DetId seedDetId = (theRecoEl[j].superCluster()->seed()->hitsAndFractions())[0].first;
	   bool isBarrel = (seedDetId.subdetId() == EcalBarrel);
	   eleE_         = theRecoEl[j].energy();
	   elePt_        = theRecoEl[j].pt();
	   eleEta_       = theRecoEl[j].eta();
	   elePhi_       = theRecoEl[j].phi();
	   
	   ///for getting noise crystals
	   seedIDvec.push_back(seedDetId);
	   eta.push_back(eleEta_);
	   phi.push_back(elePhi_);
	   
	   ///find the matrix of crystals in 5x5 array around the crystal
	   //https://cmssdt.cern.ch/lxr/source/RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h#0869
	   //std::vector<DetId> v_id = noZS::EcalClusterTools::matrixDetId(caloTopology_.product(), seedDetId, 2);
	   
	   /// SJ increased this to 7x7 array on 11th Dec, 2023
	   //std::vector<DetId> v_id = noZS::EcalClusterTools::matrixDetId(caloTopology_.product(), seedDetId, 3);
	   
	   //// store 5x5 only
	   std::vector<DetId> v_id = noZS::EcalClusterTools::matrixDetId(caloTopology_.product(), seedDetId, 2);
	   
	   nCrys_ = (int) v_id.size();
	   
	   ////match each crystal with simHit
	   for (const auto& id : v_id) {
	     
	     double matchSimEnergyEB = -99, matchSimEnergyEE = -99;
	     bool isMatch = getSimHitMatch(id, pEBSim, pEESim, matchSimEnergyEB, matchSimEnergyEE);
	     
	     hitsIsSimMatch_.push_back(isMatch);
	     hitsSimEnEB_.push_back(matchSimEnergyEB);
	     hitsSimEnEE_.push_back(matchSimEnergyEE);
	   }
	   
	   ////end of match each crystal with simHit
	   
	   std::vector<double> hitsEnergy;
	   std::vector<double> hitsThr;
	   std::vector<double> hitsEta;
	   std::vector<double> hitsPhi;
	   
	   std::vector<std::vector<double>> hitsAmplitudes = getTimeSamplesAroundEle(geo, v_id, pEBDigi, pEEDigi, EBRecHits, EERecHits, thresholds, hitsEnergy, hitsThr, hitsEta, hitsPhi);
	     hitsAmplitudes_ = hitsAmplitudes;
	     hitsEnergy_ = hitsEnergy;
	     hitsThr_ = hitsThr;
	     e5x5_ = theRecoEl[j].e5x5();
	     hitsEta_ = hitsEta;
	     hitsPhi_ = hitsPhi;
	   
	     
	     gendR_ = 999;
	     std::vector<reco::GenParticle>::const_iterator genMatch = getGenMatch(genLep, theRecoEl[j], gendR_);
	     
	     if(gendR_ < 999){
	       genPt_ = genMatch->pt();
	       genEta_ = genMatch->eta();
	       genPhi_ = genMatch->phi();
	       genE_ = genMatch->energy();
	       genStatus_ = genMatch->status();
	     }
	     else{
	       genPt_ = -99;
	       genEta_ = -99;
	       genPhi_ = -99;
	       genE_ = -99;
	       genStatus_ = -99;
	     }
	     
	     
	     /*
	       int icrys = 0;
	       int isample = 0;
	       for (const auto& array : hitsAmplitudes) {
	       isample = 0;
	       for (const auto& element : array) {
	       
	       //std::cout<<"Conent of hitsAmplitudes_["<<icrys<<"]["<<isample<<"] is "<<hitsAmplitudes_[icrys][isample]<<std::endl;
	       
	       isample++;
	       }
	       icrys++;
	       }
	     */
	     treeEle->Fill();
	 }///end of electron loop
       }//if (electronHandle.isValid())
     }//if(!runMinBias_)

   /// Get time samples from the region where there are no electrons or jets

   /// EB
   if (barrelRecHitsHandle.isValid()) {

     std::vector<DetId> v_id_EB;
     if(debugL1_){
       std::cout<<"EB rehitsize "<<EBRecHits->size()<<std::endl;
     }
     
     for (uint j = 0; j < EBRecHits->size(); j++) {
       
       DetId id = (*EBRecHits)[j].id();

       double en = (*EBRecHits)[j].energy();

       if(en<=0) continue;
       
       /*auto it = std::find(v_id_EB.begin(), v_id_EB.end(), id);
       if (it != v_id_EB.end()) continue; //it means that the ID is present in the vector so dont use it. 
       */
       
       const GlobalPoint & rechitPoint = geo->getPosition(id);
       double rheta = rechitPoint.eta();
       double rhphi = rechitPoint.phi();

       if(fabs(rheta) > 1.4) continue; ///not taking border hits
       bool foundNoiseHit = true;
       for(uint k = 0; k < seedIDvec.size(); k++){
	 
	 const GlobalPoint & seedPoint = geo->getPosition(seedIDvec[k]);
	 double seedeta = seedPoint.eta();
	 double seedphi = seedPoint.phi();

	 double dEta = fabs(rheta-seedeta);
	 double dPhi = fabs(rhphi-seedphi);
	 if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
	 if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
	 

  
	 if(debugL1_){
	   std::cout<<"Eta : Phi : Energy "<<rheta<<" "<<rhphi<<" "<<(*EBRecHits)[j].energy()<<std::endl;
	   std::cout<<"Seed eta and phi "<<seedeta<<" "<<seedphi<<std::endl;
	   std::cout<<"dEta and dPhi "<<dEta<<" "<<dPhi<<std::endl;
	 }

	 ///opposite eta and same phi --> most probably noise
	 if( dEta<2. || dPhi>1.5 ){
	   foundNoiseHit = false;
	   break;
	 }//if( dEta<3.14 || dPhi>3 )
       }//for(uint k = 0; k < seedIDvec.size(); k++)

       if(debugL1_){
	 std::cout<<"FoundHit "<<foundNoiseHit<<std::endl;
       }

       if(foundNoiseHit==false) continue;
       ///now for this rechit, get 7x7 rechits around it
       //v_id_EB = noZS::EcalClusterTools::matrixDetId(caloTopology_.product(), id, 3);

       ///store 5x5 
       v_id_EB = noZS::EcalClusterTools::matrixDetId(caloTopology_.product(), id, 2);

       if(debugL1_){

	 std::cout<<"size of 7x7 crystals around this noisy crystal "<<v_id_EB.size()<<std::endl;
       }
       
       std::vector<double> hitsEnergy;
       std::vector<double> hitsThr;
       std::vector<double> hitsEta;
       std::vector<double> hitsPhi;
       
       std::vector<std::vector<double>> hitsAmplitudes = getTimeSamplesAroundEle(geo, v_id_EB, pEBDigi, pEEDigi, EBRecHits, EERecHits, thresholds, hitsEnergy, hitsThr, hitsEta, hitsPhi);
       hitsNoiseEBAmplitudes_ = hitsAmplitudes;
       hitsNoiseEBEnergy_ = hitsEnergy;
       hitsNoiseEBThr_ = hitsThr;
       hitsNoiseEBEta_ = hitsEta;
       hitsNoiseEBPhi_ = hitsPhi;

       ///print what is inside the hitsNoiseEBAmplitudes_
       int icrys = 0;
       int isample = 0;
       for (const auto& array : hitsNoiseEBAmplitudes_) {
	 isample = 0;
	 for (const auto& element : array) {

	   std::cout<<"Conent of hitsNoiseEBAmplitudes_["<<icrys<<"]["<<isample<<"] is "<<hitsNoiseEBAmplitudes_[icrys][isample]<<std::endl;
	   
	   isample++;
	 }
	 icrys++;
       }
       ///print what is inside the hitsNoiseEBAmplitudes_

       treeNoiseEB->Fill();
       
     }//for (uint j = 0; j < EBRecHits.size(); j++)
   }///if (ebRecHitCollection_.isValid())


   /// EE rechit collection
   if (endcapRecHitsHandle.isValid()) {

     std::vector<DetId> v_id_EE;
     
     for (uint j = 0; j < EERecHits->size(); j++) {
       
       DetId id = (*EERecHits)[j].id();

       double en = (*EERecHits)[j].energy();
       if(en<=0) continue;
       
       auto it = std::find(v_id_EE.begin(), v_id_EE.end(), id);
       if (it != v_id_EE.end()) continue; //it means that the ID is present in the vector so dont use it. 

       
       const GlobalPoint & rechitPoint = geo->getPosition(id);
       double rheta = rechitPoint.eta();
       double rhphi = rechitPoint.phi();

       if(fabs(rheta) < 1.56 || fabs(rheta) > 2.4 ) continue; ///not taking border hits
       bool foundNoiseHit = true;
       for(uint k = 0; k < seedIDvec.size(); k++){
	 
	 const GlobalPoint & seedPoint = geo->getPosition(seedIDvec[k]);
	 double seedeta = seedPoint.eta();
	 double seedphi = seedPoint.phi();

	 double dEta = fabs(rheta-seedeta);
	 double dPhi = fabs(rhphi-seedphi);

	 ///opposite eta and same phi --> most probably noise
	 if( dEta<2. || dPhi>1.5 ){
	   foundNoiseHit = false;
	   break;
	 }//if( dEta<3.14 || dPhi>3 )
       }//for(uint k = 0; k < seedIDvec.size(); k++)

       if(foundNoiseHit==false) continue;
       ///now for this rechit, get 7x7 rechits around it
       //v_id_EE = noZS::EcalClusterTools::matrixDetId(caloTopology_.product(), id, 3);

       /// store 5x5 only
       v_id_EE = noZS::EcalClusterTools::matrixDetId(caloTopology_.product(), id, 2);
       
       std::vector<double> hitsEnergy;
       std::vector<double> hitsThr;
       std::vector<double> hitsEta;
       std::vector<double> hitsPhi;
       
       std::vector<std::vector<double>> hitsAmplitudes = getTimeSamplesAroundEle(geo, v_id_EE, pEBDigi, pEEDigi, EBRecHits, EERecHits, thresholds, hitsEnergy, hitsThr, hitsEta, hitsPhi);
       hitsNoiseEEAmplitudes_ = hitsAmplitudes;
       hitsNoiseEEEnergy_ = hitsEnergy;
       hitsNoiseEEThr_ = hitsThr;
       hitsNoiseEEEta_ = hitsEta;
       hitsNoiseEEPhi_ = hitsPhi;

       treeNoiseEE->Fill();
     }//for (uint j = 0; j < EERecHits.size(); j++)
   }///if (ebRecHitCollection_.isValid())


   
   ////cant run on photons from ALCARECO - will have to run those modules

   /*Handle<reco::PhotonCollection> theRecoPhotonCollection;
   iEvent.getByToken(recophotonCollection_, theRecoPhotonCollection);
   const reco::PhotonCollection theRecoPh = *(theRecoPhotonCollection.product());
   
   if (theRecoPhotonCollection.isValid()) {
     for (uint j = 0; j < theRecoPh.size(); j++){
       
       DetId seedDetId = (theRecoPh[j].superCluster()->seed()->hitsAndFractions())[0].first;
       bool isBarrel = (seedDetId.subdetId() == EcalBarrel);
       
       if(isBarrel){
	 //https://cmssdt.cern.ch/lxr/source/CaloOnlineTools/EcalTools/plugins/EcalCosmicsHists.cc#0109
	 EcalDigiCollection::const_iterator thisdigi = pEBDigi->find(seedDetId);
	 if (thisdigi == pEBDigi->end()){
	   std::cout<<"ECALTimeSampleAnalyser!!!  WARNING: seedDetId not found in the pEBDigi collection!"<<std::endl;
	 }
	 else{
	   EBDataFrame df(*thisdigi);
	   
	   for (unsigned int i = 0; i < (*thisdigi).size(); ++i) {
	     EcalMGPASample samp_crystal(df.sample(i));
	   }//loop over time samples
	 }//else when the id is found
       }//if(isBarrel)
       
       if(!isBarrel){
	 //https://cmssdt.cern.ch/lxr/source/CaloOnlineTools/EcalTools/plugins/EcalCosmicsHists.cc#0109
	 EcalDigiCollection::const_iterator thisdigi = pEEDigi->find(seedDetId);
	 if (thisdigi == pEEDigi->end()){
	   std::cout<<"ECALTimeSampleAnalyser!!!  WARNING: seedDetId not found in the EEDigi collection!"<<std::endl;
	 }
	 else{
	   EEDataFrame df(*thisdigi);
	   
	   for (unsigned int i = 0; i < (*thisdigi).size(); ++i) {
	     EcalMGPASample samp_crystal(df.sample(i));
	   }//loop over time samples
	 }//else when the id is found
       }//if(isBarrel)
       

     }///end of loop over photons
   }//if (photonHandle.isValid()) 
   /////////////////
   */



   /*
   if (pEBDigi) {
     for (EBDigiCollection::const_iterator digiItr = pEBDigi->begin(); digiItr != pEBDigi->end();
	  ++digiItr) {  // Loop on EB crystals
       EBDetId id_crystal(digiItr->id());
       EBDataFrame df(*digiItr);

       
       int etaG = id_crystal.ieta();  // global
       int phiG = id_crystal.iphi();  // global

       int etaL;  // local
       int phiL;  // local
       std::pair<int, int> LocalCoord = MEEBGeom::localCoord(etaG, phiG);
       etaL = LocalCoord.first;
       phiL = LocalCoord.second;
       eta = etaG;
       phi = phiG;
       side = MEEBGeom::side(etaG, phiG);
       EcalElectronicsId elecid_crystal = TheMapping.getElectronicsId(id_crystal);
       towerID = elecid_crystal.towerId();
       int strip = elecid_crystal.stripId();
       int xtal = elecid_crystal.xtalId();
       channelID = 5 * (strip - 1) + xtal - 1;  // FIXME
       int module = MEEBGeom::lmmod(etaG, phiG);
       int iMod = module - 1;
       assert(module >= *min_element(modules.begin(), modules.end()) &&
	      module <= *max_element(modules.begin(), modules.end()));
     

       std::cout<<"Size of hte digi sample "<<(*digiItr).size()<<std::endl;
       for (unsigned int i = 0; i < (*digiItr).size(); ++i) {
	 EcalMGPASample samp_crystal(df.sample(i));
	 std::cout<<"ADC of "<<i<<"th sample is "<<samp_crystal.adc()<<" and gain is "<<samp_crystal.gainId()<<std::endl;
	 
	 
	 adc[i] = samp_crystal.adc();
	 adcG[i] = samp_crystal.gainId();
	 if (i == 0)
	 adcgain = adcG[i];
	 if (i > 0)
	 adcgain = TMath::Max(adcG[i], adcgain);
	
       }
       
       // Remove pedestal
       //====================
       for (dsum = 0., dsum1 = 0., k = 0; k < _presample; k++) {
	 dsum += adc[k];
	 if (k < _presample - 1)
	   dsum1 += adc[k];
       }
       bl = dsum / ((double)_presample);
       for (val_max = 0., k = 0; k < _nsamples; k++) {
	 yrange[k] = adc[k] - bl;
	 if (yrange[k] > val_max) {
	   val_max = yrange[k];
	 }
       }
       
     }
   }
   */
}

   
// ------------ method called once each job just before starting event loop  ------------
void
ECALTimeSampleAnalyser::beginJob()
{
  edm::Service<TFileService> fs;
  treeEle    = fs->make<TTree>("EventTreeEle", "Event data");

  //
  treeEle->Branch("eleE",                    &eleE_);
  treeEle->Branch("elePt",                   &elePt_);
  treeEle->Branch("eleEta",                  &eleEta_);
  treeEle->Branch("elePhi",                  &elePhi_);
  treeEle->Branch("hitsAmplitudes",         &hitsAmplitudes_);
  treeEle->Branch("hitsEnergy",         &hitsEnergy_);
  treeEle->Branch("hitsThr",         &hitsThr_);
  treeEle->Branch("hitsEta",         &hitsEta_);
  treeEle->Branch("hitsPhi",         &hitsPhi_);

  
  treeEle->Branch("hitsIsSimMatch",         &hitsIsSimMatch_);
  treeEle->Branch("hitsSimEnEB",         &hitsSimEnEB_);
  treeEle->Branch("hitsSimEnEE",         &hitsSimEnEE_);

  treeEle->Branch("nsamples",         &nsamples_);
  treeEle->Branch("e5x5",         &e5x5_);
  treeEle->Branch("genPt",         &genPt_);
  treeEle->Branch("genEta",         &genEta_);
  treeEle->Branch("genPhi",         &genPhi_);
  treeEle->Branch("genE",         &genE_);
  treeEle->Branch("genStatus",         &genStatus_);
  treeEle->Branch("gendR",         &gendR_);


  ////Noise tree
  treeNoiseEB    = fs->make<TTree>("EventTreeNoiseEB", "Event data");
  treeNoiseEB->Branch("hitsNoiseEBAmplitudes",         &hitsNoiseEBAmplitudes_);
  treeNoiseEB->Branch("hitsNoiseEBEnergy",         &hitsNoiseEBEnergy_);
  treeNoiseEB->Branch("hitsNoiseEBThr",         &hitsNoiseEBThr_);
  treeNoiseEB->Branch("hitsNoiseEBEta",         &hitsNoiseEBEta_);
  treeNoiseEB->Branch("hitsNoiseEBPhi",         &hitsNoiseEBPhi_);

  treeNoiseEE    = fs->make<TTree>("EventTreeNoiseEE", "Event data");
  treeNoiseEE->Branch("hitsNoiseEEAmplitudes",         &hitsNoiseEEAmplitudes_);
  treeNoiseEE->Branch("hitsNoiseEEEnergy",         &hitsNoiseEEEnergy_);
  treeNoiseEE->Branch("hitsNoiseEEThr",         &hitsNoiseEEThr_);
  treeNoiseEE->Branch("hitsNoiseEEEta",         &hitsNoiseEEEta_);
  treeNoiseEE->Branch("hitsNoiseEEPhi",         &hitsNoiseEEPhi_);


}

// ------------ method called once each job just after ending the event loop  ------------
void
ECALTimeSampleAnalyser::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ECALTimeSampleAnalyser::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}


std::vector<std::vector<double>> ECALTimeSampleAnalyser::getTimeSamplesAroundEle(
										 const CaloGeometry* geo,					 std::vector<DetId> v_id, 
										 edm::Handle<EBDigiCollection> pEBDigi, 
										 edm::Handle<EEDigiCollection> pEEDigi, 
										 const EcalRecHitCollection* EBRecHits, 
										 const EcalRecHitCollection* EERecHits, 
										 const EcalPFRecHitThresholds* thresholds,
										 std::vector<double> &hitsEnergy,
										 std::vector<double> &hitsThr,
										 std::vector<double> &hitsEta,
										 std::vector<double> &hitsPhi){
  
  //// follow ecalmultifit algo as linked in teh analyse function above and do the pedestal subtraction. link also given below
  ///https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/EcalRecAlgos/src/EcalUncalibRecHitMultiFitAlgo.cc#0064
  
  ////https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitWorkerMultiFit.cc#0512

  //using sampleType = std::array<double, NMAXSAMPLES>;
  std::vector<std::vector<double>> hitTimeSamples;

  
    for (const auto& id : v_id) {
    
      float rhThres = -99.0;
      if (thresholds != nullptr) {
        rhThres = (*thresholds)[id];  // access PFRechit thresholds for noise cleaning
      }
      
      hitsThr.push_back(rhThres);
      
      bool isBarrel = (id.subdetId() == EcalBarrel); 
      const EcalPedestals::Item* aped = nullptr;
      const EcalMGPAGainRatio* aGain = nullptr;
      
      std::vector<double> amplitudes(NMAXSAMPLES);
      for(int isample=0; isample<NMAXSAMPLES; isample++){
	amplitudes[isample] = -99.;
      }
      
      double rechitEn = -99;
      
      double rheta = -99;
      double rhphi = -99;
	

      if(isBarrel){
	
	//https://cmssdt.cern.ch/lxr/source/CaloOnlineTools/EcalTools/plugins/EcalCosmicsHists.cc#

	if(debugL1_){
	  std::cout<<"Inside getRecHitsAroundElectrons"<<std::endl;
	}

	
	EcalRecHitCollection::const_iterator it = EBRecHits->find(id);

	if(it == EBRecHits->end()){
	  if(debugL1_) std::cout<<"This rechit not found in the rechit collection"<<std::endl;
	}
	
	if (it != EBRecHits->end()) {
	  if ((it->checkFlag(EcalRecHit::kTowerRecovered) || it->checkFlag(EcalRecHit::kWeird) ||
	       (it->detid().subdetId() == EcalBarrel && it->checkFlag(EcalRecHit::kDiWeird)))){

	    if(debugL1_) std::cout<<"This hit is a problematic one - checkFlag"<<std::endl;
	    rechitEn = 0.0;
	  }
	  else{
	    rechitEn = it->energy();
	  }
	} else {
	  rechitEn =  0;
	}
      
	///eta phi
	//EBDetId det = it->id();
	//EBDetId det = it->rawId();
	if(debugL1_){
	  std::cout << "is this null ID "<< id.null()<<" raw ID " << id.rawId()<<" For this hit, detID got whose energy is "<<rechitEn<<std::endl;
	}
	
	EBDetId det(id.rawId());
	
	const GlobalPoint & rechitPoint = geo->getPosition(det);
	rheta = rechitPoint.eta();
	rhphi = rechitPoint.phi();
	
	if(debugL1_){
	  std::cout<<"For this hit, eta and phi "<<rheta<<" "<<rhphi<<std::endl;
	}
	
	unsigned int hashedIndex = EBDetId(det).hashedIndex();
	aped = &peds->barrel(hashedIndex);
	aGain = &gains->barrel(hashedIndex);

	if(debugL1_) std::cout<<"aped : aGain : "<<aped<<" "<<aGain<<std::endl;
	//EcalDigiCollection::const_iterator thisdigi = pEBDigi->find(id);
	EcalDigiCollection::const_iterator thisdigi = pEBDigi->find(det);
	if (thisdigi == pEBDigi->end()){

	  if(debugL1_) std::cout<<"ECALTimeSampleAnalyser!!!  WARNING: this ID not found in the pEBDigi collection!"<<std::endl;
	
	  for(int isample=0; isample<NMAXSAMPLES; isample++){
	    amplitudes[isample] = -99.;
	  }
	  hitTimeSamples.push_back(amplitudes);
	}
	
	else{
	  EBDataFrame df(*thisdigi);
	  
	  for (unsigned int i = 0; i < (*thisdigi).size(); i++) {
	    EcalMGPASample samp_crystal(df.sample(i));
	    if(debugL1_) std::cout<<"ADC of "<<i<<"th sample is "<<samp_crystal.adc()<<" and gain is "<<samp_crystal.gainId()<<std::endl;
	    
	    double amplitude = 0.;
	    int gainId = samp_crystal.gainId();
	    
	    double pedestal = 0.;
	    double gainratio = 1.;
	    
	    
	    if (gainId == 0 || gainId == 3) {
	      pedestal = aped->mean_x1;
	      gainratio = aGain->gain6Over1() * aGain->gain12Over6();
	    } else if (gainId == 1) {
	      pedestal = aped->mean_x12;
	      gainratio = 1.;
	    } else if (gainId == 2) {
	      pedestal = aped->mean_x6;
	      gainratio = aGain->gain12Over6();
	    }
	  /// static pedestal here - no multifit at the moment
	    amplitude = ((double)(samp_crystal.adc()) - pedestal) * gainratio;
	    if (gainId == 0) {
	      //saturation
	      amplitude = (4095. - pedestal) * gainratio;
	    }
	    
	    amplitudes[i] = amplitude;
	    //std::cout<<"amplitude after sub "<<amplitude<<std::endl;
	  }///loop over time samples
	  hitTimeSamples.push_back(amplitudes);
	}//else when the id is found
	
      }//if(isBarrel)
      
      if(!isBarrel){
	//https://cmssdt.cern.ch/lxr/source/CaloOnlineTools/EcalTools/plugins/EcalCosmicsHists.cc#
	
	EcalRecHitCollection::const_iterator it = EERecHits->find(id);
	if (it != EERecHits->end()) {
	  if ((it->checkFlag(EcalRecHit::kTowerRecovered) || it->checkFlag(EcalRecHit::kWeird) ||
	       (it->detid().subdetId() == EcalBarrel && it->checkFlag(EcalRecHit::kDiWeird)))){
	    rechitEn = 0.0;
	  }
	  else{
	    rechitEn = it->energy();
	  }
	} else {
	  rechitEn =  0;
	}

	///eta phi
	//EEDetId det = it->id();
	EEDetId det(id.rawId());
	const GlobalPoint & rechitPoint = geo->getPosition(det);
	rheta = rechitPoint.eta();
	rhphi = rechitPoint.phi();
	
	unsigned int hashedIndex = EEDetId(det).hashedIndex();
	aped = &peds->endcap(hashedIndex);
	aGain = &gains->endcap(hashedIndex);
	
	EcalDigiCollection::const_iterator thisdigi = pEEDigi->find(det);
	if (thisdigi == pEEDigi->end()){
	  //std::cout<<"ECALTimeSampleAnalyser!!!  WARNING: seedDetId not found in the pEEDigi collection!"<<std::endl;
	}
	else{
	  EEDataFrame df(*thisdigi);
	  
	  for (unsigned int i = 0; i < (*thisdigi).size(); i++) {
	    EcalMGPASample samp_crystal(df.sample(i));
	    //std::cout<<"ADC of "<<i<<"th sample is "<<samp_crystal.adc()<<" and gain is "<<samp_crystal.gainId()<<std::endl;
	    
	    double amplitude = 0.;
	    int gainId = samp_crystal.gainId();
	    
	    double pedestal = 0.;
	    double gainratio = 1.;
	    
	    
	    if (gainId == 0 || gainId == 3) {
	      pedestal = aped->mean_x1;
	      gainratio = aGain->gain6Over1() * aGain->gain12Over6();
	    } else if (gainId == 1) {
	      pedestal = aped->mean_x12;
	      gainratio = 1.;
	    } else if (gainId == 2) {
	      pedestal = aped->mean_x6;
	      gainratio = aGain->gain12Over6();
	    }
	    /// static pedestal here - no multifit at the moment
	    amplitude = ((double)(samp_crystal.adc()) - pedestal) * gainratio;
	    if (gainId == 0) {
	      //saturation
	    amplitude = (4095. - pedestal) * gainratio;
	    }
	    
	    amplitudes[i] = amplitude;
	    //std::cout<<"amplitude after sub "<<amplitude<<std::endl;
	    
	  }//loop over time samples
	  hitTimeSamples.push_back(amplitudes);
	}//else when the id is found
      }//if(!isBarrel)
      hitsEnergy.push_back(rechitEn);
      hitsEta.push_back(rheta);
      hitsPhi.push_back(rhphi);
    }//for (const auto& id : v_id)
    
    return hitTimeSamples;
    
}



//genmatch 
std::vector<reco::GenParticle>::const_iterator  ECALTimeSampleAnalyser::getGenMatch(std::vector<std::vector<reco::GenParticle>::const_iterator> genLep, reco::GsfElectron gsfele, double &dRmin){

  std::vector<reco::GenParticle>::const_iterator genMatch;
    TLorentzVector *ele = new TLorentzVector();
    ele->SetPtEtaPhiM(gsfele.pt(), gsfele.eta(), gsfele.phi(), 0.511/1000.);

    dRmin = 999;
    
  for(auto& lep : genLep){
    TLorentzVector *gen = new TLorentzVector(); 
    gen->SetPtEtaPhiM(lep->pt(), lep->eta(), lep->phi(), lep->mass());
    double deltaR = ele->DeltaR(*gen);
    
    if(deltaR < dRmin){
      //std::cout<<"dR at the moment "<<deltaR<<std::endl;
      dRmin = deltaR;
      genMatch = lep;
    }
  }
  
  //std::cout<<"===final dR "<<dRmin<<std::endl;
  return genMatch;
}



bool ECALTimeSampleAnalyser::getSimHitMatch(
					      DetId crys_id,
					      edm::Handle<edm::PCaloHitContainer> pEBSim,
					      edm::Handle<edm::PCaloHitContainer> pEESim,
					      double &matchSimEnergyEB,
					      double &matchSimEnergyEE
					    ){


  bool foundHit = false;
  bool isBarrel = (crys_id.subdetId() == EcalBarrel);

  double totEn = 0;
      
  if(isBarrel){

    //EBDetId id = pEBSim->find(id);
    //double energy = pEBSim->energy();
    for (edm::PCaloHitContainer::const_iterator simItr = pEBSim->begin(); simItr != pEBSim->end(); ++simItr) {

      EBDetId id = simItr->id();
      if(id == crys_id){
	foundHit = true;
	totEn += simItr->energy();
      }
      
    }// for (PCaloHitContainer::const_iter....)
    matchSimEnergyEB = totEn;
    
  }//if(isBarrel)


  if(!isBarrel){

    //EBDetId id = pEBSim->find(id);
    //double energy = pEBSim->energy();
    for (edm::PCaloHitContainer::const_iterator simItr = pEESim->begin(); simItr != pEESim->end(); ++simItr) {
      
      EEDetId id = simItr->id();
      if(id == crys_id){
	foundHit = true;
	totEn += simItr->energy();
      }
      
    }// for (PCaloHitContainer::const_iter....)

    matchSimEnergyEE = totEn;
  }//if(!isBarrel)
  
  return foundHit;

}

//define this as a plug-in
DEFINE_FWK_MODULE(ECALTimeSampleAnalyser);
