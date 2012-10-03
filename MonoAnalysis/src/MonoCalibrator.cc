// -*- C++ -*-
//
// Package:    MonoCalibrator
// Class:      MonoCalibrator
// 
/**\class MonoCalibrator MonoCalibrator.cc Monopoles/MonoCalibrator/src/MonoCalibrator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christopher Cowden
//         Created:  Tue Feb  7 16:21:08 CST 2012
// $Id: MonoCalibrator.cc,v 1.5 2012/09/23 21:04:43 cowden Exp $
//
//


// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"


// data formats
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"


// Ecal includes
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"


// Monopole analysis includes
#include "Monopoles/MonoAlgorithms/interface/NPVHelper.h"
#include "Monopoles/MonoAlgorithms/interface/MonoEcalObs0.h"
#include "Monopoles/MonoAlgorithms/interface/EnergyFlowFunctor.h"
#include "Monopoles/MonoAlgorithms/interface/ClustCategorizer.h"


// ROOT includes
#include "TTree.h"
#include "TBranch.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF2.h"


//
// class declaration
//

class MonoCalibrator : public edm::EDAnalyzer {

   typedef std::vector<reco::Photon> PhotonCollection;
   //typedef std::vector<reco::Electron> ElectronCollection;
   typedef std::vector<reco::GsfElectron> ElectronCollection;
   typedef std::vector<reco::BasicCluster> BasicClusterCollection;


   public:
      explicit MonoCalibrator(const edm::ParameterSet&);
      ~MonoCalibrator();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

    // clear tree variables
    void clear();

    // fill delta R maps 
    template <class S, class T>
    void fillDRMap( const S &, const T &, std::vector<double> *);




      // ----------member data ---------------------------

    // input tags
    edm::InputTag m_TagEcalEB_RecHits;
    edm::InputTag m_Tag_Jets;
    edm::InputTag m_Tag_Photons;
    edm::InputTag m_Tag_Electrons;

    // Monopole Ecal Observables
    //Mono::MonoEcalObs0 m_ecalObs;
    Mono::MonoEcalObs0Calibrator m_ecalCalib;

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
MonoCalibrator::MonoCalibrator(const edm::ParameterSet& iConfig)
  :m_TagEcalEB_RecHits(iConfig.getParameter<edm::InputTag>("EcalEBRecHits") )
  ,m_Tag_Jets(iConfig.getParameter<edm::InputTag>("JetTag") )
  ,m_Tag_Photons(iConfig.getParameter<edm::InputTag>("PhotonTag") )
  ,m_Tag_Electrons(iConfig.getParameter<edm::InputTag>("ElectronTag") )
  //,m_ecalObs(iConfig)
  ,m_ecalCalib(iConfig)
{

}


MonoCalibrator::~MonoCalibrator()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MonoCalibrator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  clear();

  // execute calibration for event
  //m_ecalCalib.calculateMijn(iSetup,iEvent);
  m_ecalCalib.fillClust(iSetup,iEvent);

}


// ------------ method called once each job just before starting event loop  ------------
void 
MonoCalibrator::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MonoCalibrator::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MonoCalibrator::beginRun(edm::Run const&, edm::EventSetup const&)
{

}



void MonoCalibrator::clear()
{ 

}


template <class S, class T>
inline void MonoCalibrator::fillDRMap(const S &a, const T &bcoll, std::vector<double> *map )
{

  assert(map);

  map->resize(bcoll->size(),0.);

  for ( unsigned i=0; i != bcoll->size(); i++ ) 
    (*map)[i] = reco::deltaR(a,bcoll->at(i));
  

  std::sort(map->begin(),map->end() ); 

}


// ------------ method called when ending the processing of a run  ------------
void 
MonoCalibrator::endRun(edm::Run const&, edm::EventSetup const&)
{
  m_ecalCalib.calculateHij();
  m_ecalCalib.dumpCalibration();

}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MonoCalibrator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MonoCalibrator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MonoCalibrator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MonoCalibrator);
