// -*- C++ -*-
//
// Package:    MonoAnalysis
// Class:      MonoAnalysis
// 
/**\class MonoAnalysis MonoAnalysis.cc Monopoles/MonoAnalysis/src/MonoAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christopher Cowden
//         Created:  Tue Feb  7 16:21:08 CST 2012
// $Id: MonoAnalysis.cc,v 1.1 2012/04/25 17:55:19 cowden Exp $
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

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"


// Ecal includes
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

// Monopole analysis includes
#include "Monopoles/MonoAlgorithms/interface/NPVHelper.h"


// ROOT includes
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TH2D.h"


//
// class declaration
//

class MonoAnalysis : public edm::EDAnalyzer {

   typedef std::vector<reco::Photon> PhotonCollection;


   public:
      explicit MonoAnalysis(const edm::ParameterSet&);
      ~MonoAnalysis();

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


    // TFileService
    edm::Service<TFileService> m_fs;


    TTree * m_tree;

    // Event information
    unsigned m_run;
    unsigned m_lumi;
    unsigned m_event;

    unsigned m_NPV;

    // Jet information
    unsigned m_jet_N;
    std::vector<double> m_jet_E;
    std::vector<double> m_jet_p;
    std::vector<double> m_jet_pt;
    std::vector<double> m_jet_px;
    std::vector<double> m_jet_py;
    std::vector<double> m_jet_pz;
    std::vector<double> m_jet_eta;
    std::vector<double> m_jet_phi;


    // Photon information
    unsigned m_pho_N;
    std::vector<double> m_pho_E;
    std::vector<double> m_pho_p;
    std::vector<double> m_pho_pt;
    std::vector<double> m_pho_px;
    std::vector<double> m_pho_py;
    std::vector<double> m_pho_pz;
    std::vector<double> m_pho_eta;
    std::vector<double> m_pho_phi;


    // Ecal RecHits
    std::vector<double> m_ehit_eta;
    std::vector<double> m_ehit_phi;
    std::vector<double> m_ehit_time;
    std::vector<double> m_ehit_energy;
    std::vector<double> m_ehit_otEnergy;
    std::vector<double> m_ehit_jetIso;
    std::vector<double> m_ehit_phoIso;






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
MonoAnalysis::MonoAnalysis(const edm::ParameterSet& iConfig)
  :m_TagEcalEB_RecHits(iConfig.getParameter<edm::InputTag>("EcalEBRecHits") )
  ,m_Tag_Jets(iConfig.getParameter<edm::InputTag>("JetTag") )
  ,m_Tag_Photons(iConfig.getParameter<edm::InputTag>("PhotonTag") )

{
   //now do what ever initialization is needed

}


MonoAnalysis::~MonoAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MonoAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  clear();

  m_run = iEvent.id().run();
  m_lumi = iEvent.id().luminosityBlock();
  m_event = iEvent.id().event();

  // get NPV for this event
  m_NPV = Mono::getNPV(iEvent,iSetup);


  // get jet collection
  Handle<reco::PFJetCollection> jets;
  iEvent.getByLabel(m_Tag_Jets,jets);
  

  // fill jet branches
  reco::PFJetCollection::const_iterator itJet = jets->begin();
  for ( ; itJet != jets->end(); itJet++ ) {

    m_jet_E.push_back( (*itJet).energy() );
    m_jet_p.push_back( (*itJet).p() );
    m_jet_pt.push_back( (*itJet).pt() );
    m_jet_px.push_back( (*itJet).px() );
    m_jet_py.push_back( (*itJet).py() );
    m_jet_pz.push_back( (*itJet).pz() );
    m_jet_eta.push_back( (*itJet).eta() );
    m_jet_phi.push_back( (*itJet).phi() );

    m_jet_N++;

  }




  // get photon collection
  Handle<PhotonCollection> photons;
  iEvent.getByLabel(m_Tag_Photons,photons);

  // fill photon branches
  PhotonCollection::const_iterator itPho = photons->begin();
  for ( ; itPho != photons->end(); itPho++ ) {

    m_pho_E.push_back( (*itPho).energy() );
    m_pho_p.push_back( (*itPho).p() );
    m_pho_pt.push_back( (*itPho).pt() );
    m_pho_px.push_back( (*itPho).px() );
    m_pho_py.push_back( (*itPho).py() );
    m_pho_pz.push_back( (*itPho).pz() );
    m_pho_eta.push_back( (*itPho).eta() );
    m_pho_phi.push_back( (*itPho).phi() );

    m_pho_N++;

  }


  // get RecHit collection
  Handle<EBRecHitCollection > ecalRecHits;
  iEvent.getByLabel(m_TagEcalEB_RecHits,ecalRecHits);
  assert( ecalRecHits->size() > 0 );


  ESHandle<CaloGeometry> calo;
  iSetup.get<CaloGeometryRecord>().get(calo);
  const CaloGeometry *m_caloGeo = (const CaloGeometry*)calo.product();
  const CaloSubdetectorGeometry *geom = m_caloGeo->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);

  std::vector<double> dRPhotons;
  std::vector<double> dRJets;

  // fill RecHit branches
  EBRecHitCollection::const_iterator itHit = ecalRecHits->begin();
  for ( ; itHit != ecalRecHits->end(); itHit++ ) {

    EBDetId detId( (*itHit).id() );
    const CaloCellGeometry *cell = geom->getGeometry( detId );

    m_ehit_eta.push_back( cell->getPosition().eta() );
    m_ehit_phi.push_back( cell->getPosition().phi() );
    m_ehit_energy.push_back( (*itHit).energy() );
    m_ehit_time.push_back( (*itHit).time() );
    m_ehit_otEnergy.push_back( (*itHit).outOfTimeEnergy() );

    fillDRMap(cell->getPosition(),photons,&dRPhotons);
    fillDRMap(cell->getPosition(),jets,&dRJets);

    m_ehit_jetIso.push_back( dRJets.size() > 0 ? dRJets[0] : -1. );
    m_ehit_phoIso.push_back( dRPhotons.size() > 0 ? dRPhotons[0] : -1. );

  }




  m_tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
MonoAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MonoAnalysis::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MonoAnalysis::beginRun(edm::Run const&, edm::EventSetup const&)
{


  m_tree = m_fs->make<TTree>("tree","tree");

  m_tree->Branch("run",&m_run,"run/i");
  m_tree->Branch("lumiBlock",&m_lumi,"lumiBlock/i");
  m_tree->Branch("event",&m_event,"evnet/i");

  m_tree->Branch("NPV",&m_NPV,"NPV/i");

  m_tree->Branch("jet_N",&m_jet_N,"jet_N/i");
  m_tree->Branch("jet_E",&m_jet_E);
  m_tree->Branch("jet_p",&m_jet_p);
  m_tree->Branch("jet_pt",&m_jet_pt);
  m_tree->Branch("jet_px",&m_jet_px);
  m_tree->Branch("jet_py",&m_jet_py);
  m_tree->Branch("jet_pz",&m_jet_pz);
  m_tree->Branch("jet_eta",&m_jet_eta);
  m_tree->Branch("jet_phi",&m_jet_phi);


  m_tree->Branch("pho_N",&m_pho_N,"pho_N/i");
  m_tree->Branch("pho_E",&m_pho_E);
  m_tree->Branch("pho_p",&m_pho_p);
  m_tree->Branch("pho_pt",&m_pho_pt);
  m_tree->Branch("pho_px",&m_pho_px);
  m_tree->Branch("pho_py",&m_pho_py);
  m_tree->Branch("pho_pz",&m_pho_pz);
  m_tree->Branch("pho_eta",&m_pho_eta);
  m_tree->Branch("pho_phi",&m_pho_phi);

  m_tree->Branch("ehit_eta",&m_ehit_eta);
  m_tree->Branch("ehit_phi",&m_ehit_phi);
  m_tree->Branch("ehit_time",&m_ehit_time);
  m_tree->Branch("ehit_energy",&m_ehit_energy);
  m_tree->Branch("ehit_otEnergy",&m_ehit_otEnergy);
  m_tree->Branch("ehit_jetIso",&m_ehit_jetIso);
  m_tree->Branch("ehit_phoIso",&m_ehit_phoIso);
  


}


void MonoAnalysis::clear()
{ 


     m_run = 0;
     m_lumi = 0;
     m_event = 0;
  
    m_NPV = 0;

    // Jet information
    m_jet_N = 0;
    m_jet_E.clear();
    m_jet_p.clear();
    m_jet_pt.clear();
    m_jet_px.clear();
    m_jet_py.clear();
    m_jet_pz.clear();
    m_jet_eta.clear();
    m_jet_phi.clear();


    // Photon information
    m_pho_N = 0;
    m_pho_E.clear();
    m_pho_p.clear();
    m_pho_pt.clear();
    m_pho_px.clear();
    m_pho_py.clear();
    m_pho_pz.clear();
    m_pho_eta.clear();
    m_pho_phi.clear();


    // Ecal RecHits
    m_ehit_eta.clear();
    m_ehit_phi.clear();
    m_ehit_time.clear();
    m_ehit_energy.clear();
    m_ehit_otEnergy.clear();
    m_ehit_jetIso.clear();
    m_ehit_phoIso.clear();


}


template <class S, class T>
inline void MonoAnalysis::fillDRMap(const S &a, const T &bcoll, std::vector<double> *map )
{

  assert(map);

  map->resize(bcoll->size(),0.);

  for ( unsigned i=0; i != bcoll->size(); i++ ) 
    (*map)[i] = reco::deltaR(a,bcoll->at(i));
  

  std::sort(map->begin(),map->end() ); 

}


// ------------ method called when ending the processing of a run  ------------
void 
MonoAnalysis::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MonoAnalysis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MonoAnalysis::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MonoAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MonoAnalysis);
