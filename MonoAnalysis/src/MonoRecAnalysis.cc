// -*- C++ -*-
//
// Package:    MonoRecAnalysis
// Class:      MonoRecAnalysis
// 
/**\class MonoRecAnalysis MonoRecAnalysis.cc Monopoles/MonoRecAnalysis/src/MonoRecAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christopher Cowden
//         Created:  Tue Feb  14 16:21:08 CST 2012
// $Id: MonoRecAnalysis.cc,v 1.1 2012/04/25 17:55:19 cowden Exp $
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// data formats
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/Common/interface/SortedCollection.h"


// Ecal includes
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"


// For trigger
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"


// Monopole analysis includes
#include "Monopoles/MonoAlgorithms/interface/EcalMapper.h"
#include "Monopoles/MonoAlgorithms/interface/MonoTruthSnooper.h"
#include "Monopoles/MonoAlgorithms/interface/MonoGenTrackExtrapolator.h"


// ROOT includes
#include "TH2D.h"
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TLatex.h"


//
// class declaration
//


class MonoRecAnalysis : public edm::EDAnalyzer {



   public:
      explicit MonoRecAnalysis(const edm::ParameterSet&);
      ~MonoRecAnalysis();

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

    void triggerAnalysis(const edm::Event &, const edm::EventSetup &);
    bool triggerDecision(const std::string &trigger, const edm::Event &);

      // ----------member data ---------------------------
    
    // edm InputTags to retrieve reco collections
    edm::InputTag m_TagEcalEB_RecHits;


   std::vector<std::string> m_triggers;

    // trigger utils
    L1GtUtils m_L1Utils;


    // TFileService
    edm::Service<TFileService> m_fs;


    TH2D * m_h_2_tVsE;

    TTree * m_tree;

    unsigned m_run;
    unsigned m_lumi;
    unsigned m_event; 
    std::vector<double> m_EcalTime;
    std::vector<double> m_EcalE; 
    std::vector<double> m_EcalEta;
    std::vector<double> m_EcalPhi;
    double m_m1_eta;
    double m_m1_phi;
    double m_m2_eta;
    double m_m2_phi;
    std::vector<double> m_d1_eta;
    std::vector<double> m_d1_phi;
    std::vector<int> m_d1_pdgId;
    std::vector<double> m_d2_eta;
    std::vector<double> m_d2_phi;
    std::vector<int> m_d2_pdgId;

    double m_monoExp_eta;
    double m_monoExp_phi;
    double m_monoExp_g;
    double m_monoExp_m;
    double m_monoExp_z;

    double m_amonExp_eta;
    double m_amonExp_phi;
    double m_amonExp_g;
    double m_amonExp_m;
    double m_amonExp_z;

    std::vector<int> m_triggerDecision;


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
MonoRecAnalysis::MonoRecAnalysis(const edm::ParameterSet& iConfig)
  :m_TagEcalEB_RecHits(iConfig.getParameter<edm::InputTag>("EcalEBRecHits") )
  ,m_triggers(iConfig.getParameter<std::vector<std::string> >("triggerNames") )
{
   //now do what ever initialization is needed

  m_triggerDecision.resize(m_triggers.size(),0);

}


MonoRecAnalysis::~MonoRecAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MonoRecAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


  clear();

  // put event information in TTree branch variables
  m_run = iEvent.id().run();
  m_lumi = iEvent.id().luminosityBlock();
  m_event = iEvent.id().event();


  triggerAnalysis(iEvent,iSetup);
  

  // get handle on EcalEB RecHits
  Handle<EBRecHitCollection > EcalRecHits;
  iEvent.getByLabel(m_TagEcalEB_RecHits,EcalRecHits);

  assert( EcalRecHits->size() > 0 );


  Mono::MonoTruthSnoop snoopy(iEvent,iSetup);
  const HepMC::GenParticle *mono = snoopy.mono(Mono::monopole);
  const HepMC::GenParticle *anti = snoopy.mono(Mono::anti_monopole);


  m_m1_eta = mono ? mono->momentum().eta() : -999.;
  m_m1_phi = mono ? mono->momentum().phi() : -999.;
  m_m2_eta = anti ? anti->momentum().eta() : -999.;
  m_m2_phi = anti ? anti->momentum().phi() : -999.;
  
  std::vector<const HepMC::GenParticle *> d1 = snoopy.monoDaughter(Mono::monopole);
  std::vector<const HepMC::GenParticle *> d2 = snoopy.monoDaughter(Mono::anti_monopole);

  for ( unsigned i=0; i < d1.size(); i++ ) {
    std::cout << "Found monopole daughter!" << std::endl;
    m_d1_eta.push_back( d1[i]->momentum().eta() );
    m_d1_phi.push_back( d1[i]->momentum().phi() );
    m_d1_pdgId.push_back( d1[i]->pdg_id() );
  }

  for ( unsigned i=0; i < d2.size(); i++ ) {
    std::cout << "Found monopole daughter!" << std::endl;
    m_d2_eta.push_back( d2[i]->momentum().eta() );
    m_d2_phi.push_back( d2[i]->momentum().phi() );
    m_d2_pdgId.push_back( d2[i]->pdg_id() );
  }

  Mono::MonoGenTrackExtrapolator extrap;

  if ( mono ) {
    extrap.setMonopole(*mono);
    m_monoExp_eta = extrap.etaVr(1.29);
    m_monoExp_phi = extrap.phi();
    m_monoExp_g = extrap.charge();
    m_monoExp_m = extrap.mass();
    m_monoExp_z = extrap.zVr(1.29);
  }

  if ( anti ) {
    extrap.setMonopole(*anti);
    m_amonExp_eta = extrap.etaVr(1.29);
    m_amonExp_phi = extrap.phi();
    m_amonExp_g = extrap.charge();
    m_amonExp_m = extrap.mass(); 
    m_amonExp_z = extrap.zVr(1.29);
  }



  Mono::EcalMapper ecalMapper(iSetup);

  //ecalMapper.setMarker(mono->momentum().eta(),mono->momentum().phi(),"M");
  //ecalMapper.setMarker(anti->momentum().eta(),anti->momentum().phi(),"A");

  ecalMapper.fillMap(iEvent);

  edm::ESHandle<CaloGeometry> calo;
  iSetup.get<CaloGeometryRecord>().get(calo);
  const CaloGeometry *m_caloGeo = (const CaloGeometry*)calo.product();
  const CaloSubdetectorGeometry *geom = m_caloGeo->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);



  // loop over Hits
  EBRecHitCollection::const_iterator itHit = EcalRecHits->begin();
  for ( ; itHit != EcalRecHits->end(); itHit++ ) {

    EBDetId detId( (*itHit).id() );
    const CaloCellGeometry *cell = geom->getGeometry( detId );

    m_h_2_tVsE->Fill((*itHit).energy(),(*itHit).time());

    m_EcalTime.push_back( (*itHit).time() );
    m_EcalE.push_back( (*itHit).energy() );  
    m_EcalEta.push_back( cell->getPosition().eta() );
    m_EcalPhi.push_back( cell->getPosition().phi() );

  }


  m_tree->Fill();


}


// triggerAnalysis method
void
MonoRecAnalysis::triggerAnalysis( const edm::Event &ev, const edm::EventSetup &es )
{
  

  m_L1Utils.retrieveL1EventSetup(es);
  std::cout << "Trigger Menu: " << m_L1Utils.l1TriggerMenu() << std::endl;
  std::cout << "\timplementation " << m_L1Utils.l1TriggerMenuImplementation() << std::endl;


  for ( unsigned i=0; i != m_triggers.size(); i++ ) {

    bool decision = triggerDecision( m_triggers[i], ev );
    m_triggerDecision[i] = (int)decision;

  }



}

bool
MonoRecAnalysis::triggerDecision( const std::string &triggerName, const edm::Event &ev )
{

  int errorCode = -1;


  bool decision = m_L1Utils.decision(ev,triggerName,errorCode);



  if (errorCode == 0) {
   
    return decision; 

  } else if (errorCode == 1) {
        
    // algorithm / technical trigger  does not exist in the L1 menu
    // do something
    throw cms::Exception("TriggerNotIncluded")
	<< triggerName << " is not included in the current menu: " 
	<< m_L1Utils.l1TriggerMenu() << std::endl;
    

  } else {

    // error - see error code
    // do whatever needed
    throw cms::Exception("UnkownTriggerAnalysisError")
	<< "Encountered unknown trigger error with the following trigger: "
	<< triggerName << std::endl;

  }

  return false;

}


// ------------ method called once each job just before starting event loop  ------------
void 
MonoRecAnalysis::beginJob()
{

 m_h_2_tVsE = m_fs->make<TH2D>("tVsE","tVsE",50,0.,20.,50,-15.,30.); 

  m_tree = m_fs->make<TTree>("tree","tree");

  m_tree->Branch("EcalTime",&m_EcalTime);
  m_tree->Branch("EcalEnergy",&m_EcalE);
  m_tree->Branch("EcalEta",&m_EcalEta);
  m_tree->Branch("EcalPhi",&m_EcalPhi);
  m_tree->Branch("mono_eta",&m_m1_eta);
  m_tree->Branch("mono_phi",&m_m1_phi);
  m_tree->Branch("monoExp_eta",&m_monoExp_eta);
  m_tree->Branch("monoExp_z",&m_monoExp_z);
  m_tree->Branch("monoExp_phi",&m_monoExp_phi);
  m_tree->Branch("monoExp_g",&m_monoExp_g);
  m_tree->Branch("monoExp_m",&m_monoExp_m);
  m_tree->Branch("anti_eta",&m_m2_eta);
  m_tree->Branch("anti_phi",&m_m2_phi);
  m_tree->Branch("antiExp_eta",&m_amonExp_eta);
  m_tree->Branch("antiExp_z",&m_amonExp_z);
  m_tree->Branch("antiExp_phi",&m_amonExp_phi);
  m_tree->Branch("antiExp_g",&m_amonExp_g);
  m_tree->Branch("antiExp_m",&m_amonExp_m); 
  m_tree->Branch("md_eta",&m_d1_eta);
  m_tree->Branch("md_phi",&m_d1_phi);
  m_tree->Branch("md_pdgId",&m_d1_pdgId);
  m_tree->Branch("ad_eta",&m_d2_eta);
  m_tree->Branch("ad_phi",&m_d2_phi);
  m_tree->Branch("ad_pdgId",&m_d2_pdgId);
  m_tree->Branch("run",&m_run);
  m_tree->Branch("lumi",&m_lumi);
  m_tree->Branch("event",&m_event);


  char name[50];
  for ( unsigned i=0; i != m_triggers.size(); i++ ) {
    sprintf(name,"%s/I",m_triggers[i].c_str()); 
    m_tree->Branch(m_triggers[i].c_str(),&m_triggerDecision[i],name);
  }


}

// ------------ method called once each job just after ending the event loop  ------------
void 
MonoRecAnalysis::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MonoRecAnalysis::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MonoRecAnalysis::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MonoRecAnalysis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MonoRecAnalysis::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MonoRecAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


void MonoRecAnalysis::clear() 
{

  m_m1_eta = 0.;
  m_m1_phi = 0.;
  m_m2_eta = 0.;
  m_m2_phi = 0.;

  m_monoExp_eta = 0.;
  m_monoExp_phi = 0.;
  m_monoExp_g = 0.;
  m_monoExp_m = 0.;
  m_amonExp_eta = 0.;
  m_amonExp_phi = 0.;
  m_amonExp_g = 0.;
  m_amonExp_m = 0.;

  m_run = 0.;
  m_lumi = 0.;
  m_event = 0.;

  m_EcalTime.clear();
  m_EcalE.clear();
  m_EcalEta.clear();
  m_EcalPhi.clear();

  m_d1_eta.clear();
  m_d1_phi.clear();
  m_d1_pdgId.clear();
  m_d2_eta.clear();
  m_d2_phi.clear();
  m_d2_pdgId.clear();


  m_triggerDecision.resize(m_triggers.size(),0);

}

//define this as a plug-in
DEFINE_FWK_MODULE(MonoRecAnalysis);
