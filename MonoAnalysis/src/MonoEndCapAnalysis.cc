// -*- C++ -*-
//
// Package:    MonoEndCapAnalysis
// Class:      MonoEndCapAnalysis
// 
/**\class MonoEndCapAnalysis MonoEndCapAnalysis.cc Monopoles/MonoEndCapAnalysis/src/MonoEndCapAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christopher Cowden
//         Created:  Tue Feb  7 16:21:08 CST 2012
// $Id: MonoEndCapAnalysis.cc,v 1.2 2013/06/08 04:24:43 cowden Exp $
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

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

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
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

// Monopole algorithms includes
#include "Monopoles/MonoAlgorithms/interface/MonoEcalObs0.h"

// ROOT includes
#include "TTree.h"
#include "TBranch.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF2.h"


//
// class declaration
//

class MonoEndCapAnalysis : public edm::EDAnalyzer {

   typedef std::vector<reco::Photon> PhotonCollection;
   //typedef std::vector<reco::Electron> ElectronCollection;
   typedef std::vector<reco::GsfElectron> ElectronCollection;
   typedef std::vector<reco::BasicCluster> BasicClusterCollection;


   public:
      explicit MonoEndCapAnalysis(const edm::ParameterSet&);
      ~MonoEndCapAnalysis();

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

    // get the number of primary vertices
    inline unsigned getNPV(const edm::Event &, const edm::EventSetup &);

    // fill delta R maps 
    template <class S, class T>
    void fillDRMap( const S &, const T &, std::vector<double> *);



      // ----------member data ---------------------------

    // EndCap map limits
    static const unsigned s_nBins=100;
    static constexpr double  s_xMax=100.5;
    static constexpr double  s_xMin=0.5;

    static constexpr double s_EEz=3.144;

    // input tags
    edm::InputTag m_TagEcalEE_RecHits;
    bool m_isData;

    /////////////////////////////
    // map between detId,x,y,eta,phi
    struct EEDetIdMap {
	EEDetId id;
	double x,y,eta,phi;
    };
    std::vector<EEDetIdMap> m_idMap;


    // TFileService
    edm::Service<TFileService> m_fs;

    TTree * m_tree;

    // Event information
    unsigned m_run;
    unsigned m_lumi;
    unsigned m_event;

    unsigned m_NPV;

    // Ecal RecHits
    std::vector<double> m_ehit_eta;
    std::vector<double> m_ehit_phi;
    std::vector<double> m_ehit_time;
    std::vector<double> m_ehit_energy;
    std::vector<double> m_ehit_otEnergy;
    std::vector<double> m_ehit_kWeird;
    std::vector<double> m_ehit_kDiWeird;
    std::vector<double> m_ehit_jetIso;
    std::vector<double> m_ehit_phoIso;

    std::vector<double> m_ehit_x;
    std::vector<double> m_ehit_y;
    std::vector<double> m_ehit_z;

    std::vector<double> m_check_plus;
    std::vector<double> m_check_minus;


    // Ecal hybrid clusters
    unsigned m_nClusterEgamma;
    std::vector<double> m_egClust_E;
    std::vector<double> m_egClust_size;
    std::vector<double> m_egClust_eta;
    std::vector<double> m_egClust_phi;
    std::vector<double> m_egClust_frac51;
    std::vector<double> m_egClust_frac15;
    std::vector<double> m_egClust_e55;
    std::vector<double> m_egClust_eMax;
    std::vector<double> m_egClust_matchDR;
    std::vector<double> m_egClust_tagged;
    std::vector<double> m_egClust_matchPID;

    unsigned m_nCleanClusterEgamma;
    std::vector<double> m_cleanClust_E;
    std::vector<double> m_cleanClust_size;
    std::vector<double> m_cleanClust_eta;
    std::vector<double> m_cleanClust_phi;
    std::vector<double> m_cleanClust_frac51;
    std::vector<double> m_cleanClust_frac15;
    std::vector<double> m_cleanClust_e55;
    std::vector<double> m_cleanClust_eMax;
    std::vector<double> m_cleanClust_matchDR;
    std::vector<double> m_cleanClust_tagged;
    std::vector<double> m_cleanClust_matchPID;

    // generator monopoles
    double m_mono_eta;
    double m_mono_phi;
    double m_mono_m;
    double m_mono_pT;

    double m_amon_eta;
    double m_amon_phi;
    double m_amon_m;
    double m_amon_pT;

    // extrapolated monopoles
    double m_monoExp_x;
    double m_monoExp_y;
    double m_monoExp_z;
    double m_monoExp_t;
    double m_monoExp_r;
    double m_monoExp_phi;
    double m_monoExp_eta;
    double m_monoExp_ix;
    double m_monoExp_iy;
    
    double m_amonExp_x;
    double m_amonExp_y;
    double m_amonExp_z;
    double m_amonExp_t;
    double m_amonExp_r;
    double m_amonExp_phi;
    double m_amonExp_eta;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
const unsigned MonoEndCapAnalysis::s_nBins;
constexpr double MonoEndCapAnalysis::s_xMax;
constexpr double MonoEndCapAnalysis::s_xMin;
constexpr double MonoEndCapAnalysis::s_EEz;


//
// constructors and destructor
//
MonoEndCapAnalysis::MonoEndCapAnalysis(const edm::ParameterSet& iConfig)
  :m_TagEcalEE_RecHits(iConfig.getParameter<edm::InputTag>("EcalEERecHits") )
  ,m_isData(iConfig.getParameter<bool>("isData") )
{ 
  if(CLHEP::electron_charge==0) std::cout << "asdf" << std::endl;  
  
}



MonoEndCapAnalysis::~MonoEndCapAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MonoEndCapAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  clear();

  m_run = iEvent.id().run();
  m_lumi = iEvent.id().luminosityBlock();
  m_event = iEvent.id().event();

  // get NPV for this event
  m_NPV = getNPV(iEvent,iSetup);

  //////////////////////////////////////////////////////
  // create histogram maps for end-caps
  //////////////////////////////////////////////////////

  EEDetId detId();

  char histName[50];
  sprintf(histName,"EEplus_%d",m_event);
  //TH2D * h_plus = m_fs->make<TH2D>(histName,histName,nBinsEta,etaMin,etaMax,nBinsPhi,-M_PI,M_PI);

  sprintf(histName,"EEminus_%d",m_event);
  //TH2D * h_minus = m_fs->make<TH2D>(histName,histName,nBinsEta,etaMin,etaMax,nBinsPhi,-M_PI,M_PI);


  ///////////////////////////////////////////////////////
  // get RecHit collection
  //////////////////////////////////////////////////////
  Handle<EERecHitCollection > ecalRecHits;
  iEvent.getByLabel(m_TagEcalEE_RecHits,ecalRecHits);
  assert( ecalRecHits->size() > 0 );

  // create histograms for End cap maps
//  char name[50];
  /*sprintf(name,"eMapp_%d",m_event);
  TH2D * hEmapp = m_fs->make<TH2D>(name,"",s_nBins,s_xMin,s_xMax,s_nBins,s_xMin,s_xMax);

  sprintf(name,"tMapp_%d",m_event);
  TH2D * hTmapp = m_fs->make<TH2D>(name,"",s_nBins,s_xMin,s_xMax,s_nBins,s_xMin,s_xMax);

  sprintf(name,"eMapn_%d",m_event);
  TH2D * hEmapn = m_fs->make<TH2D>(name,"",s_nBins,s_xMin,s_xMax,s_nBins,s_xMin,s_xMax);

  sprintf(name,"tMapn_%d",m_event);
  TH2D * hTmapn = m_fs->make<TH2D>(name,"",s_nBins,s_xMin,s_xMax,s_nBins,s_xMin,s_xMax);

  sprintf(name,"checkp_%d",m_event);
  TH2D * hCheckp = m_fs->make<TH2D>(name,"",s_nBins,s_xMin,s_xMax,s_nBins,s_xMin,s_xMax);
  sprintf(name,"checkn_%d",m_event);
  TH2D * hCheckn = m_fs->make<TH2D>(name,"",s_nBins,s_xMin,s_xMax,s_nBins,s_xMin,s_xMax);*/

  ESHandle<CaloGeometry> calo;
  iSetup.get<CaloGeometryRecord>().get(calo);
  const CaloGeometry *m_caloGeo = (const CaloGeometry*)calo.product();
  const CaloSubdetectorGeometry *geom = m_caloGeo->getSubdetectorGeometry(DetId::Ecal,EcalEndcap);

  // fill RecHit branches
  EERecHitCollection::const_iterator itHit = ecalRecHits->begin();
  for ( ; itHit != ecalRecHits->end(); itHit++ ) {

    EEDetId detId( (*itHit).id() );
    const CaloCellGeometry *cell = geom->getGeometry( detId );

    m_ehit_eta.push_back( cell->getPosition().eta() );
    m_ehit_phi.push_back( cell->getPosition().phi() );
    m_ehit_energy.push_back( (*itHit).energy() );
    m_ehit_time.push_back( (*itHit).time() );
    // the outOfTimeEnergy method has been removed from the EcalRecHit class
    // in CMSSW_7.  I leave this commment here is a note/reminder this is something
    // I don't immediately know how to fix, but this analyzer is not used much
    // so it doesn't need to be fixed at the moment.
    //m_ehit_otEnergy.push_back( (*itHit).outOfTimeEnergy() );

    m_ehit_kWeird.push_back( (*itHit).checkFlag(EcalRecHit::kWeird) );
    m_ehit_kDiWeird.push_back( (*itHit).checkFlag(EcalRecHit::kDiWeird) );

    m_ehit_x.push_back( cell->getPosition().x() );
    m_ehit_y.push_back( cell->getPosition().y() );
    m_ehit_z.push_back( cell->getPosition().z() );

    /*if (detId.zside() > 0 ) {
      int bin = hEmapp->FindBin(detId.ix(),detId.iy());
      hEmapp->SetBinContent(bin,(*itHit).energy());
      hTmapp->SetBinContent(bin,(*itHit).time());
      hCheckp->SetBinContent(bin,hCheckp->GetBinContent(bin)+1);
    } else {
    
      int bin = hEmapn->FindBin(detId.ix(),detId.iy());
      hEmapn->SetBinContent(bin,(*itHit).energy());
      hTmapn->SetBinContent(bin,(*itHit).time());
      hCheckn->SetBinContent(bin,hCheckn->GetBinContent(bin)+1);

    }*/
  }

  /*for ( unsigned b=1; b <= s_nBins*s_nBins; b++ ){
    m_check_plus.push_back( hCheckp->GetBinContent(b) );
    m_check_minus.push_back( hCheckn->GetBinContent(b) );
  }*/

  // get a handle on topology
  ESHandle<CaloTopology> topo;
  iSetup.get<CaloTopologyRecord>().get(topo);
  const CaloTopology * topology = (const CaloTopology*)topo.product();


  // get BasicCluster Collection
  Handle<BasicClusterCollection> bClusters;
  edm::InputTag bcClusterTag("multi5x5SuperClusters","uncleanOnlyMulti5x5EndcapBasicClusters"); 
  iEvent.getByLabel(bcClusterTag,bClusters);
  const unsigned nbClusters = bClusters->size();

  EcalClusterTools ecalTool;

  // cluster tagger
  Mono::GenMonoClusterTagger tagger(0.3);
  if ( !m_isData ) {
    tagger.initialize(iEvent,iSetup);
    if ( nbClusters) tagger.tag(nbClusters,&(*bClusters)[0]);
  } 

  for ( unsigned i=0; i != nbClusters; i++ ) {
    m_egClust_E.push_back( (*bClusters)[i].energy() );
    m_egClust_size.push_back( (*bClusters)[i].size() );
    m_egClust_eta.push_back( (*bClusters)[i].eta() );
    m_egClust_phi.push_back( (*bClusters)[i].phi() );

    const float e55 = ecalTool.e5x5((*bClusters)[i],ecalRecHits.product(),topology);
    const float e51 = ecalTool.e5x1((*bClusters)[i],ecalRecHits.product(),topology);
    const float e15 = ecalTool.e1x5((*bClusters)[i],ecalRecHits.product(),topology);
    const float eMax = ecalTool.eMax((*bClusters)[i],ecalRecHits.product());
    m_egClust_frac51.push_back( e51/e55 );
    m_egClust_frac15.push_back( e15/e55 );
    m_egClust_e55.push_back(e55);
    m_egClust_eMax.push_back(eMax/e55);

    if ( !m_isData ) {
      m_egClust_matchDR.push_back(tagger.matchDR()[i]);
      m_egClust_tagged.push_back(tagger.tagResult()[i]);
      m_egClust_matchPID.push_back(tagger.matchPID()[i]);
    }
  }
  m_nClusterEgamma = nbClusters;



  // get BasicCluster Collection
  Handle<BasicClusterCollection> cleanClusters;
  edm::InputTag cleanClusterTag("multi5x5SuperClusters","multi5x5EndcapBasicClusters"); 
  iEvent.getByLabel(cleanClusterTag,cleanClusters);
  const unsigned nCleanClusters = cleanClusters->size();

  // cluster tagger
  tagger.clearTags();
  if ( !m_isData && nCleanClusters ) tagger.tag(nCleanClusters,&(*cleanClusters)[0]);

  for ( unsigned i=0; i != nCleanClusters; i++ ) {
    m_cleanClust_E.push_back( (*cleanClusters)[i].energy() );
    m_cleanClust_size.push_back( (*cleanClusters)[i].size() );
    m_cleanClust_eta.push_back( (*cleanClusters)[i].eta() );
    m_cleanClust_phi.push_back( (*cleanClusters)[i].phi() );

    const float e55 = ecalTool.e5x5((*cleanClusters)[i],ecalRecHits.product(),topology);
    const float e51 = ecalTool.e5x1((*cleanClusters)[i],ecalRecHits.product(),topology);
    const float e15 = ecalTool.e1x5((*cleanClusters)[i],ecalRecHits.product(),topology);
    const float eMax = ecalTool.eMax((*cleanClusters)[i],ecalRecHits.product());
    m_cleanClust_frac51.push_back( e51/e55 );
    m_cleanClust_frac15.push_back( e15/e55 );
    m_cleanClust_e55.push_back(e55);
    m_cleanClust_eMax.push_back(eMax/e55);

    if ( !m_isData ) {
      m_cleanClust_matchDR.push_back(tagger.matchDR()[i]);
      m_cleanClust_tagged.push_back(tagger.tagResult()[i]);
      m_cleanClust_matchPID.push_back(tagger.matchPID()[i]);
    }
  }
  m_nCleanClusterEgamma = nCleanClusters;


  //////////////////////////////////////
  // fill in monopole MC information
  if ( !m_isData ) {
    Mono::MonoTruthSnoop snoopy(iEvent,iSetup);
    const HepMC::GenParticle *mono = snoopy.mono(Mono::monopole);
    const HepMC::GenParticle *amon = snoopy.mono(Mono::anti_monopole);

    Mono::MonoGenTrackExtrapolator extrap;

    if ( mono ) {
      m_mono_eta = mono->momentum().eta();
      m_mono_phi = mono->momentum().phi();
      m_mono_pT = mono->momentum().perp();
      m_mono_m = mono->momentum().m();

      extrap.setMonopole(*mono);
      double tmp = extrap.rVz(s_EEz);
      m_monoExp_z = s_EEz;
      m_monoExp_phi = extrap.phi();
      if ( tmp == -1. ) {
	m_monoExp_z = -s_EEz;
	tmp = extrap.rVz(-s_EEz);
      }
      m_monoExp_r = tmp;
      m_monoExp_t = extrap.tVz(m_monoExp_z);
      m_monoExp_eta = extrap.eta(m_monoExp_z,m_monoExp_r);
      m_monoExp_x = m_monoExp_r*cos(m_monoExp_phi);
      m_monoExp_y = m_monoExp_r*sin(m_monoExp_phi);

      const GlobalPoint p(m_monoExp_x,m_monoExp_y,m_monoExp_z);
      const DetId monoId =  geom->getClosestCell( p );
      const EEDetId eeId(monoId);
      m_monoExp_ix = eeId.ix();
      m_monoExp_iy = eeId.iy();
    }

    if ( amon ) {
      m_amon_eta = amon->momentum().eta();
      m_amon_phi = amon->momentum().phi();
      m_amon_pT = amon->momentum().perp();
      m_amon_m = amon->momentum().m();

      extrap.setMonopole(*amon);
      double tmp = extrap.rVz(s_EEz);
      m_amonExp_phi = extrap.phi();
      m_amonExp_z = s_EEz;
      if ( tmp == -1 ) {
	m_amonExp_z = -s_EEz;
	tmp = extrap.rVz(-s_EEz);
      }
      m_amonExp_r = tmp;
      m_amonExp_t = extrap.tVz(m_amonExp_z);
      m_amonExp_eta = extrap.eta(m_amonExp_z,m_amonExp_r);
      m_amonExp_x = m_amonExp_r*cos(m_amonExp_phi);
      m_amonExp_y = m_amonExp_r*sin(m_amonExp_phi);
    }

  }

  ////////////////////////////////////////
  //
  m_tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
MonoEndCapAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MonoEndCapAnalysis::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MonoEndCapAnalysis::beginRun(edm::Run const& run, edm::EventSetup const& es)
{

  // setup geometry stuff
  edm::ESHandle<CaloGeometry> calo;
  es.get<CaloGeometryRecord>().get(calo);
  const CaloGeometry *m_caloGeo = (const CaloGeometry*)calo.product();
  const CaloSubdetectorGeometry *geom = m_caloGeo->getSubdetectorGeometry(DetId::Ecal,EcalEndcap);

  const std::vector<DetId> & valIds = geom->getValidDetIds();
  TH1D *hix = m_fs->make<TH1D>("hix","ix",100,0.5,100.5);
  TH1D *hiy = m_fs->make<TH1D>("hiy","iy",100,0.5,100.5);
  TH2D *hi = m_fs->make<TH2D>("hi","ix iy",100,0.5,100.5,100,0.5,100.5);
  for ( unsigned i=0; i != valIds.size(); i++ ) {
    EEDetId eeId(valIds[i]);
    hix->Fill(eeId.ix());
    hiy->Fill(eeId.iy());
    hi->Fill(eeId.ix(),eeId.iy());

    const CaloCellGeometry *cell = geom->getGeometry( eeId );
    const GlobalPoint & pos = cell->getPosition();
    
    EEDetIdMap map;
    map.id = eeId;
    map.x = pos.x();
    map.y = pos.y();
    map.eta = pos.eta();
    map.phi = pos.phi();

    m_idMap.push_back(map);
  }
  

  // setup tree variables
  m_tree = m_fs->make<TTree>("tree","tree");

  m_tree->Branch("run",&m_run,"run/i");
  m_tree->Branch("lumiBlock",&m_lumi,"lumiBlock/i");
  m_tree->Branch("event",&m_event,"evnet/i");

  m_tree->Branch("NPV",&m_NPV,"NPV/i");

  m_tree->Branch("ehit_eta",&m_ehit_eta);
  m_tree->Branch("ehit_phi",&m_ehit_phi);
  m_tree->Branch("ehit_time",&m_ehit_time);
  m_tree->Branch("ehit_E",&m_ehit_energy);
  m_tree->Branch("ehit_kWeird",&m_ehit_kWeird);
  m_tree->Branch("ehit_kDiWeird",&m_ehit_kDiWeird);

  m_tree->Branch("ehit_x",&m_ehit_x);
  m_tree->Branch("ehit_y",&m_ehit_y);
  m_tree->Branch("ehit_z",&m_ehit_z);

  m_tree->Branch("check_plus",&m_check_plus);
  m_tree->Branch("check_minus",&m_check_minus);

  m_tree->Branch("egClust_N",&m_nClusterEgamma,"egClust_N/i");
  m_tree->Branch("egClust_E",&m_egClust_E);
  m_tree->Branch("egClust_size",&m_egClust_size);
  m_tree->Branch("egClust_eta",&m_egClust_eta);
  m_tree->Branch("egClust_phi",&m_egClust_phi);
  m_tree->Branch("egClust_frac51",&m_egClust_frac51);
  m_tree->Branch("egClust_frac15",&m_egClust_frac15);
  m_tree->Branch("egClust_eMax",&m_egClust_eMax);
  m_tree->Branch("egClust_e55",&m_egClust_e55);
  m_tree->Branch("egClust_matchDR",&m_egClust_matchDR);
  m_tree->Branch("egClust_matchPID",&m_egClust_matchPID);
  m_tree->Branch("egClust_tagged",&m_egClust_tagged);

  m_tree->Branch("cleanClust_N",&m_nCleanClusterEgamma,"cleanClust_N/i");
  m_tree->Branch("cleanClust_E",&m_cleanClust_E);
  m_tree->Branch("cleanClust_size",&m_cleanClust_size);
  m_tree->Branch("cleanClust_eta",&m_cleanClust_eta);
  m_tree->Branch("cleanClust_phi",&m_cleanClust_phi);
  m_tree->Branch("cleanClust_frac51",&m_cleanClust_frac51);
  m_tree->Branch("cleanClust_frac15",&m_cleanClust_frac15);
  m_tree->Branch("cleanClust_eMax",&m_cleanClust_eMax);
  m_tree->Branch("cleanClust_e55",&m_cleanClust_e55);
  m_tree->Branch("cleanClust_matchDR",&m_cleanClust_matchDR);
  m_tree->Branch("cleanClust_matchPID",&m_cleanClust_matchPID);
  m_tree->Branch("cleanClust_tagged",&m_cleanClust_tagged);

  m_tree->Branch("mono_eta",&m_mono_eta,"mono_eta/D");
  m_tree->Branch("mono_phi",&m_mono_phi,"mono_phi/D");
  m_tree->Branch("mono_m",&m_mono_m,"mono_m/D");
  m_tree->Branch("mono_pT",&m_mono_pT,"mono_pT/D");
  
  m_tree->Branch("amon_eta",&m_amon_eta,"amon_eta/D");
  m_tree->Branch("amon_phi",&m_amon_phi,"amon_phi/D");
  m_tree->Branch("amon_m",&m_amon_m,"amon_m/D");
  m_tree->Branch("amon_pT",&m_amon_pT,"amon_pT/D");

  m_tree->Branch("monoExp_x",&m_monoExp_x,"monoExp_x/D");
  m_tree->Branch("monoExp_y",&m_monoExp_y,"monoExp_y/D");
  m_tree->Branch("monoExp_z",&m_monoExp_z,"monoExp_z/D");
  m_tree->Branch("monoExp_t",&m_monoExp_t,"monoExp_t/D");
  m_tree->Branch("monoExp_r",&m_monoExp_r,"monoExp_r/D");
  m_tree->Branch("monoExp_phi",&m_monoExp_phi,"monoExp_phi/D");
  m_tree->Branch("monoExp_eta",&m_monoExp_eta,"monoExp_eta/D");
  m_tree->Branch("monoExp_ix",&m_monoExp_ix,"monoExp_ix/D");
  m_tree->Branch("monoExp_iy",&m_monoExp_iy,"monoExp_iy/D");

  m_tree->Branch("amonExp_x",&m_amonExp_x,"amonExp_x/D");
  m_tree->Branch("amonExp_y",&m_amonExp_y,"amonExp_y/D");
  m_tree->Branch("amonExp_z",&m_amonExp_z,"amonExp_z/D");
  m_tree->Branch("amonExp_t",&m_amonExp_t,"amonExp_t/D");
  m_tree->Branch("amonExp_r",&m_amonExp_r,"amonExp_r/D");
  m_tree->Branch("amonExp_phi",&m_amonExp_phi,"amonExp_phi/D");
  m_tree->Branch("amonExp_eta",&m_amonExp_eta,"amonExp_eta/D");

}



void MonoEndCapAnalysis::clear()
{ 


     m_run = 0;
     m_lumi = 0;
     m_event = 0;
  
    m_NPV = 0;



    m_nClusterEgamma = 0;
    m_egClust_E.clear();
    m_egClust_size.clear();
    m_egClust_eta.clear();
    m_egClust_phi.clear();
    m_egClust_frac51.clear();
    m_egClust_frac15.clear();
    m_egClust_e55.clear();
    m_egClust_eMax.clear();
    m_egClust_matchDR.clear();
    m_egClust_matchPID.clear();
    m_egClust_tagged.clear();

    m_nCleanClusterEgamma = 0;
    m_cleanClust_E.clear();
    m_cleanClust_size.clear();
    m_cleanClust_eta.clear();
    m_cleanClust_phi.clear();
    m_cleanClust_frac51.clear();
    m_cleanClust_frac15.clear();
    m_cleanClust_e55.clear();
    m_cleanClust_eMax.clear();
    m_cleanClust_matchDR.clear();
    m_cleanClust_matchPID.clear();
    m_cleanClust_tagged.clear();

    m_ehit_x.clear();
    m_ehit_y.clear();
    m_ehit_z.clear();

    m_check_plus.clear();
    m_check_minus.clear();

    // Ecal RecHits
    m_ehit_eta.clear();
    m_ehit_phi.clear();
    m_ehit_time.clear();
    m_ehit_energy.clear();
    m_ehit_otEnergy.clear();
    m_ehit_kWeird.clear();
    m_ehit_kDiWeird.clear();
    m_ehit_jetIso.clear();
    m_ehit_phoIso.clear();

  m_mono_eta = 0.;
  m_mono_phi = 0.;
  m_mono_m = 0.;
  m_mono_pT = 0.;

  m_amon_eta = 0.;
  m_amon_phi = 0.;
  m_amon_m = 0.;
  m_amon_pT = 0.;

  m_monoExp_x = 0.;
  m_monoExp_y = 0.;
  m_monoExp_z = 0.;
  m_monoExp_t = 0.;
  m_monoExp_r = 0.;
  m_monoExp_phi = 0.;
  m_monoExp_eta = 0.;
  m_monoExp_ix = 0.;
  m_monoExp_iy = 0.;

  m_amonExp_x = 0.;
  m_amonExp_y = 0.;
  m_amonExp_z = 0.;
  m_amonExp_t = 0.;
  m_amonExp_r = 0.;
  m_amonExp_phi = 0.;
  m_amonExp_eta = 0.;

}


template <class S, class T>
inline void MonoEndCapAnalysis::fillDRMap(const S &a, const T &bcoll, std::vector<double> *map )
{

  assert(map);

  map->resize(bcoll->size(),0.);

  for ( unsigned i=0; i != bcoll->size(); i++ ) 
    (*map)[i] = reco::deltaR(a,bcoll->at(i));
  

  std::sort(map->begin(),map->end() ); 

}


// ------------ method called when ending the processing of a run  ------------
void 
MonoEndCapAnalysis::endRun(edm::Run const&, edm::EventSetup const&)
{

}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MonoEndCapAnalysis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MonoEndCapAnalysis::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MonoEndCapAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


inline unsigned int MonoEndCapAnalysis::getNPV(const edm::Event & iEvent, const edm::EventSetup &iSetup )
{
  // PrimaryVertex analysis
  edm::InputTag m_PVTag = edm::InputTag("offlinePrimaryVertices","");

  edm::Handle<reco::VertexCollection> handlePV;
  iEvent.getByLabel(m_PVTag,handlePV);

  int totalNPV = 0.;

  reco::VertexCollection::const_iterator pv = handlePV->begin();
  for ( ; pv != handlePV->end(); pv++ ) {
    if ( !pv->isFake() && pv->ndof() > 4.0 ) {
      ++totalNPV;
    }
  }

  return totalNPV;
}



//define this as a plug-in
DEFINE_FWK_MODULE(MonoEndCapAnalysis);
