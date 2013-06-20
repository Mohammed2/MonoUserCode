// -*- C++ -*-
//
// Package:    MonoNtupleDumper
// Class:      MonoNtupleDumper
// 
/**\class MonoNtupleDumper MonoNtupleDumper.cc Monopoles/MonoNtupleDumper/src/MonoNtupleDumper.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christopher Cowden
//         Created:  Tue Feb  7 16:21:08 CST 2012
// $Id: MonoNtupleDumper.cc,v 1.7 2013/06/13 21:45:27 cowden Exp $
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
#include "DataFormats/METReco/interface/PFMET.h"


// Ecal includes
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"


// Monopole analysis includes
#include "Monopoles/MonoAlgorithms/interface/NPVHelper.h"
#include "Monopoles/MonoAlgorithms/interface/MonoEcalObs0.h"
#include "Monopoles/MonoAlgorithms/interface/ClustCategorizer.h"
#include "Monopoles/MonoAlgorithms/interface/MonoTrackMatcher.h"

#include "Monopoles/TrackCombiner/interface/MplTracker.h"

// ROOT includes
#include "TTree.h"
#include "TBranch.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF2.h"


//
// class declaration
//

class MonoNtupleDumper : public edm::EDAnalyzer {

   typedef std::vector<reco::BasicCluster> BasicClusterCollection;
   typedef std::vector<reco::Photon> PhotonCollection;
   //typedef std::vector<reco::Electron> ElectronCollection;
   typedef std::vector<reco::GsfElectron> ElectronCollection;
   typedef reco::GsfElectron Electron;
   typedef reco::Photon Photon;


   public:
      explicit MonoNtupleDumper(const edm::ParameterSet&);
      ~MonoNtupleDumper();

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


      // ----------member data ---------------------------

    std::string m_output;
    TFile *m_outputFile;

    // input tags
    edm::InputTag m_TagEcalEB_RecHits;
    edm::InputTag m_Tag_Jets;
    edm::InputTag m_Tag_Photons;
    edm::InputTag m_Tag_Electrons;
    edm::InputTag m_Tag_MET;
    bool m_isData;

    // Monopole Ecal Observables
    Mono::MonoEcalObs0 m_ecalObs;

    //Tracker
    MplTracker *_Tracker;

    // TFileService
    //edm::Service<TFileService> m_fs;

    // map cluster category (lengthxwidth thing) to a histogram
    // showing the average energy or time  in each cell
    //TFileDirectory *m_avgDir;
    std::map<Mono::ClustCategorizer,TH2D *> m_clustEMap;
    std::map<Mono::ClustCategorizer,TH2D *> m_clustTMap;
    std::map<Mono::ClustCategorizer,unsigned> m_clustCatCount;

    std::vector<double> m_betas;
    std::vector<double> m_betaTs;

    TTree * m_tree;

    // Event information
    unsigned m_run;
    unsigned m_lumi;
    unsigned m_event;

    unsigned m_NPV;

    // Ecal Observable information
    unsigned m_nClusters;
    std::vector<double> m_clust_E;
    std::vector<double> m_clust_eta;
    std::vector<double> m_clust_phi;
    std::vector<double> m_clust_L;
    std::vector<double> m_clust_W;
    std::vector<double> m_clust_N;
    std::vector<double> m_clust_sigEta;
    std::vector<double> m_clust_sigPhi;
    std::vector<double> m_clust_meanEta;
    std::vector<double> m_clust_meanPhi;
    std::vector<double> m_clust_skewEta;
    std::vector<double> m_clust_skewPhi;
    std::vector<double> m_clust_seedFrac;
    std::vector<double> m_clust_firstFrac;
    std::vector<double> m_clust_secondFrac;
    std::vector<double> m_clust_thirdFrac;
    std::vector<double> m_clust_phiStripFrac;
    std::vector<double> m_clust_matchDR;
    std::vector<double> m_clust_tagged;
    std::vector<double> m_clust_matchPID;
    std::vector<double> m_clust_matchTime;
    std::vector<double> m_clust_matchPt;
    std::vector<double> m_clust_hsE;
    std::vector<double> m_clust_hsTime;
    std::vector<int>    m_clust_hsInSeed;
    std::vector<int>    m_clust_hsWeird;
    std::vector<int>    m_clust_hsDiWeird;
    // Treat these arrays as 3D arrays
    // There is space for 15 clusters of 100 total elements in each cluster
    // One must use m_nClusters, m_clust_L, and m_clust_W when unpacking
    // the cluster from the TTree.
    static const unsigned WS = 100;
    static const unsigned SS = 15*100;
    double m_clust_Ecells[1500]; 
    double m_clust_Tcells[1500]; 



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

    // Ecal hybrid clusters (cleaned collection)
    unsigned m_nCleanEgamma;
    std::vector<double> m_egClean_E;
    std::vector<double> m_egClean_size;
    std::vector<double> m_egClean_eta;
    std::vector<double> m_egClean_phi;
    std::vector<double> m_egClean_frac51;
    std::vector<double> m_egClean_frac15;
    std::vector<double> m_egClean_e55;
    std::vector<double> m_egClean_eMax;
    std::vector<double> m_egClean_matchDR;
    std::vector<double> m_egClean_tagged;
    std::vector<double> m_egClean_matchPID;

    // Ecal RecHits
    std::vector<double> m_ehit_eta;
    std::vector<double> m_ehit_phi;
    std::vector<double> m_ehit_time;
    std::vector<double> m_ehit_energy;
    std::vector<double> m_ehit_otEnergy;
    std::vector<double> m_ehit_flag;
    std::vector<double> m_ehit_kWeird;
    std::vector<double> m_ehit_kDiWeird;


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
    std::vector<double> m_jet_matchDR;
    std::vector<int>  m_jet_tagged;
    std::vector<int>  m_jet_matchPID;


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
    std::vector<double> m_pho_matchDR;
    std::vector<int>  m_pho_tagged;
    std::vector<int>  m_pho_matchPID;

    // Electron information
    unsigned m_ele_N;
    std::vector<double> m_ele_E;
    std::vector<double> m_ele_p;
    std::vector<double> m_ele_pt;
    std::vector<double> m_ele_px;
    std::vector<double> m_ele_py;
    std::vector<double> m_ele_pz;
    std::vector<double> m_ele_eta;
    std::vector<double> m_ele_phi;
    std::vector<double> m_ele_matchDR;
    std::vector<int>  m_ele_tagged;
    std::vector<int>  m_ele_matchPID;

    // MET information
    double m_mpt;
    double m_mpPhi;

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
MonoNtupleDumper::MonoNtupleDumper(const edm::ParameterSet& iConfig)
  :m_TagEcalEB_RecHits(iConfig.getParameter<edm::InputTag>("EcalEBRecHits") )
  ,m_Tag_Jets(iConfig.getParameter<edm::InputTag>("JetTag") )
  ,m_Tag_Photons(iConfig.getParameter<edm::InputTag>("PhotonTag") )
  ,m_Tag_Electrons(iConfig.getParameter<edm::InputTag>("ElectronTag") )
  ,m_Tag_MET(iConfig.getParameter<edm::InputTag>("METTag") )
  ,m_isData(iConfig.getParameter<bool>("isData") )
  ,m_ecalObs(iConfig)
  ,m_output(iConfig.getParameter<std::string>("Output"))
{
  _Tracker = new MplTracker(iConfig);
}


MonoNtupleDumper::~MonoNtupleDumper()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MonoNtupleDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  clear();

  m_run = iEvent.id().run();
  m_lumi = iEvent.id().luminosityBlock();
  m_event = iEvent.id().event();

  /////////////////////////////////////
  // get NPV for this event
  m_NPV = Mono::getNPV(iEvent,iSetup);


  // execute observable calculations
  double monoObs = m_ecalObs.calculate(iSetup,iEvent,&m_betas,&m_betaTs);
  const Mono::EBmap & ebMap = m_ecalObs.ecalMap();

  // limits in eta and phi of ebmap
  const unsigned nEta = ebMap.nEta();
  const unsigned nPhi = ebMap.nPhi();

  /////////////////////////////////////
  // cluster analysis
  // retrieve the clusterBuilder from the ecal obs.
  const Mono::ClusterBuilder clusterBuilder = m_ecalObs.clusterBuilder();
  m_nClusters = clusterBuilder.nClusters();

  // tag clusters to gen level monopole extrapolation
  Mono::GenMonoClusterTagger tagger(0.3);
  if ( !m_isData ) {
    tagger.initialize(iEvent,iSetup);
    if ( m_nClusters ) tagger.tag(m_nClusters,clusterBuilder.clusters(),ebMap);
  }


  // cycle of all clusters found
  for ( unsigned i=0; i != m_nClusters; i++ ) {
    const Mono::MonoEcalCluster & cluster = clusterBuilder.clusters()[i];

    // get basic cluster info
    const unsigned width = cluster.clusterWidth();
    const unsigned length = cluster.clusterLength();
    const unsigned cEta = cluster.ieta();
    const unsigned cPhi = cluster.iphi();

    // get MC monopole tag info
    if ( !m_isData ) {
      m_clust_matchDR.push_back(tagger.matchDR()[i]);
      m_clust_tagged.push_back(tagger.tagResult()[i]);
      m_clust_matchPID.push_back(tagger.matchPID()[i]);
      m_clust_matchTime.push_back(tagger.matchTime()[i]);
      m_clust_matchPt.push_back(tagger.matchPt()[i]);
    }

    m_clust_E.push_back( cluster.clusterEnergy() );
    m_clust_eta.push_back( ebMap.eta(cEta) );
    m_clust_phi.push_back( ebMap.phi(cPhi) );
    m_clust_L.push_back( length );
    m_clust_W.push_back( width );


    // prepare 2D histograms to facilitate some basic analysis on the clusters
    // like skewness
    const unsigned wings = width/2U;
    char histName[50];
    sprintf(histName,"clustHist_%d_%d",iEvent.id().event(),i);
    TH2D hist(histName,histName,length,-(float)length/2.,(float)length/2.,width,-(int)wings,wings);
    sprintf(histName,"clustTHist_%d_%d",iEvent.id().event(),i);
    TH2D Thist(histName,histName,length,0,length,width,-(int)wings,wings);

    double histMax=0.;
    double hsTime=-1.;
    int phiBin=UINT_MAX;
    bool kWeird=false;
    bool kDiWeird=false;

    // fill in cluster energy and time maps
    const bool exceedsWS = length*width > WS;
    const bool exceedsSS = length*width+i*WS > SS;
    const bool excessive = exceedsWS || exceedsSS;
    for ( unsigned j=0; j != width; j++ ) {
      int ji = (int)j-(int)width/2;
      for ( unsigned k=0; k != length; k++ ) {
        const double energy = cluster.energy(k,ji,ebMap);
	if ( energy > histMax ) {
	  histMax = energy;
	  hsTime = cluster.time(k,ji,ebMap);
	  phiBin = ji;
       	  kWeird = cluster.getRecHit(k,ji,ebMap)->checkFlag( EcalRecHit::kWeird );
       	  kDiWeird = cluster.getRecHit(k,ji,ebMap)->checkFlag( EcalRecHit::kDiWeird );
	}
	hist.SetBinContent(k+1,j+1,energy);
	Thist.SetBinContent(k+1,j+1,cluster.time(k,ji,ebMap));
	if ( !excessive ) m_clust_Ecells[i*WS+j*length+k] = cluster.energy(k,ji,ebMap);
	if ( !excessive ) m_clust_Tcells[i*WS+j*length+k] = cluster.time(k,ji,ebMap);
      }
    }


    m_clust_sigEta.push_back( hist.GetRMS(1) );
    m_clust_sigPhi.push_back( hist.GetRMS(2) );
    m_clust_skewEta.push_back( hist.GetSkewness(1) );
    m_clust_skewPhi.push_back( hist.GetSkewness(2) );

    const double clustE = cluster.clusterEnergy();
    m_clust_seedFrac.push_back( cluster.clusterSeed().energy()/clustE );
    const unsigned center = length/2U;
    m_clust_firstFrac.push_back( cluster.energy(center,0,ebMap)/clustE );
    m_clust_secondFrac.push_back( (cluster.energy(center+1U,0,ebMap)+cluster.energy(center-1U,0,ebMap))/clustE );
    m_clust_thirdFrac.push_back( (cluster.energy(center+2U,0,ebMap)+cluster.energy(center-2U,0,ebMap))/clustE );

    m_clust_hsE.push_back(histMax/clustE);
    m_clust_hsTime.push_back(hsTime);
    m_clust_hsInSeed.push_back(phiBin);
    m_clust_hsWeird.push_back( kWeird );
    m_clust_hsDiWeird.push_back( kDiWeird );

    // perform Gaussian fit to cluster
    // normalise to total energy of cluster
    hist.Scale( 1./cluster.clusterEnergy() );


    // fill in aggregate cluster information maps
    Mono::ClustCategorizer cluCat(length,width);
    std::map<Mono::ClustCategorizer,TH2D*>::iterator enIter = m_clustEMap.find(cluCat);
    std::map<Mono::ClustCategorizer,TH2D*>::iterator tmIter = m_clustTMap.find(cluCat);
    bool foundEn = enIter == m_clustEMap.end();
    bool foundTm = tmIter == m_clustTMap.end();
    assert( foundEn == foundTm ); // assert maps are somewhat synchronous

    // if category not found create the histogram
    if ( foundEn ) {
      char name[50];
      sprintf(name,"avgEnClust_%d_%d",cluCat.length,cluCat.width);
      m_clustEMap[cluCat] = new TH2D(name,name,cluCat.length,0,cluCat.length,cluCat.width,-(int)cluCat.width/2,(int)cluCat.width/2);
      sprintf(name,"avgTmClust_%d_%d",cluCat.length,cluCat.width);
      m_clustTMap[cluCat] = new TH2D(name,name,cluCat.length,0,cluCat.length,cluCat.width,-(int)cluCat.width/2,(int)cluCat.width/2);
    }

    m_clustCatCount[cluCat]++;
    TH2D * avgEnMap = m_clustEMap[cluCat];
    TH2D * avgTmMap = m_clustTMap[cluCat];

    assert( avgEnMap );
    assert( avgTmMap );

    for ( int binx = 1; binx <= hist.GetNbinsX(); binx++ ) {
      for ( int biny = 1; biny <= hist.GetNbinsY(); biny++ ) {
      	avgEnMap->SetBinContent(binx,biny,avgEnMap->GetBinContent(binx,biny)+hist.GetBinContent(binx,biny));
      	avgTmMap->SetBinContent(binx,biny,avgTmMap->GetBinContent(binx,biny)+Thist.GetBinContent(binx,biny));
      }
    }
    

  }





  // get RecHit collection
  Handle<EBRecHitCollection > ecalRecHits;
  iEvent.getByLabel(m_TagEcalEB_RecHits,ecalRecHits);
  assert( ecalRecHits->size() > 0 );


  // get calo geometry and topology
  ESHandle<CaloGeometry> calo;
  iSetup.get<CaloGeometryRecord>().get(calo);
  const CaloGeometry *m_caloGeo = (const CaloGeometry*)calo.product();
  const CaloSubdetectorGeometry *geom = m_caloGeo->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);

  ESHandle<CaloTopology> topo;
  iSetup.get<CaloTopologyRecord>().get(topo);
  const CaloTopology * topology = (const CaloTopology*)topo.product();

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

    m_ehit_kWeird.push_back( (*itHit).checkFlag(EcalRecHit::kWeird) );
    m_ehit_kDiWeird.push_back( (*itHit).checkFlag(EcalRecHit::kDiWeird) );
    m_ehit_flag.push_back( (*itHit).recoFlag() );

  }


  // get BasicCluster Collection
  Handle<BasicClusterCollection> bClusters;
  edm::InputTag bcClusterTag("hybridSuperClusters","uncleanOnlyHybridBarrelBasicClusters"); 
  iEvent.getByLabel(bcClusterTag,bClusters);
  const unsigned nbClusters = bClusters->size();

  EcalClusterTools ecalTool;
  std::vector<int> exclFlags;
  std::vector<int> sevExcl;

  tagger.clearTags();
  if ( !m_isData && nbClusters ) tagger.tag(nbClusters,&(*bClusters)[0]);

  unsigned nClusterCount=0;
  for ( unsigned i=0; i != nbClusters; i++ ) {
    if ( (*bClusters)[i].energy() < 50. ) continue;

    nClusterCount++;
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
  m_nClusterEgamma = nClusterCount;

  // get BasicCluster Collection (cleaned)
  Handle<BasicClusterCollection> cClusters;
  edm::InputTag ccClusterTag("hybridSuperClusters","hybridBarrelBasicClusters"); 
  iEvent.getByLabel(ccClusterTag,cClusters);
  const unsigned ncClusters = cClusters->size();

  tagger.clearTags();
  if ( !m_isData && ncClusters ) tagger.tag(ncClusters,&(*cClusters)[0]);

  nClusterCount=0;
  for ( unsigned i=0; i != ncClusters; i++ ) {
    if ( (*cClusters)[i].energy() < 50. ) continue;

    nClusterCount++;
    m_egClean_E.push_back( (*cClusters)[i].energy() );
    m_egClean_size.push_back( (*cClusters)[i].size() );
    m_egClean_eta.push_back( (*cClusters)[i].eta() );
    m_egClean_phi.push_back( (*cClusters)[i].phi() );

    const float e55 = ecalTool.e5x5((*cClusters)[i],ecalRecHits.product(),topology);
    const float e51 = ecalTool.e5x1((*cClusters)[i],ecalRecHits.product(),topology);
    const float e15 = ecalTool.e1x5((*cClusters)[i],ecalRecHits.product(),topology);
    const float eMax = ecalTool.eMax((*cClusters)[i],ecalRecHits.product());
    m_egClean_frac51.push_back( e51/e55 );
    m_egClean_frac15.push_back( e15/e55 );
    m_egClean_e55.push_back(e55);
    m_egClean_eMax.push_back(eMax/e55);

    if ( !m_isData ) {
      m_egClean_matchDR.push_back(tagger.matchDR()[i]);
      m_egClean_tagged.push_back(tagger.tagResult()[i]);
      m_egClean_matchPID.push_back(tagger.matchPID()[i]);
    }
  }
  m_nCleanEgamma = nClusterCount;

  ////////////////////////////////
  // Tracking analysis
  _Tracker->analyze(iEvent, iSetup);
  const Mono::MonoEcalCluster * clusters = clusterBuilder.clusters(); 
  _Tracker->doMatch(m_nClusters,clusters,ebMap);


  // get jet collection
  Handle<reco::PFJetCollection> jets;
  iEvent.getByLabel(m_Tag_Jets,jets);
  const unsigned nJets = jets->size();

  // fill jet branches
  tagger.clearTags();
  if ( !m_isData && nJets ) tagger.tag(nJets,&(*jets)[0]);
  for ( unsigned i=0; i != nJets; i++ ) {

    const reco::PFJet & jet = (*jets)[i];

    m_jet_E.push_back( jet.energy() );
    m_jet_p.push_back( jet.p() );
    m_jet_pt.push_back( jet.pt() );
    m_jet_px.push_back( jet.px() );
    m_jet_py.push_back( jet.py() );
    m_jet_pz.push_back( jet.pz() );
    m_jet_eta.push_back( jet.eta() );
    m_jet_phi.push_back( jet.phi() );

    if ( !m_isData ) {
      m_jet_matchDR.push_back( tagger.matchDR()[i] );
      m_jet_tagged.push_back(tagger.tagResult()[i]);
      m_jet_matchPID.push_back(tagger.matchPID()[i]);
    }

  }
  m_jet_N = nJets;

  // get photon collection
  Handle<PhotonCollection> photons;
  iEvent.getByLabel(m_Tag_Photons,photons);
  const unsigned nPhotons = photons->size();

  // fill photon branches
  tagger.clearTags();
  if ( !m_isData && nPhotons ) tagger.tag(nPhotons,&(*photons)[0]); 
  for ( unsigned i=0; i != nPhotons; i++ ) {
    
    const Photon & pho = (*photons)[i];

    m_pho_E.push_back( pho.energy() );
    m_pho_p.push_back( pho.p() );
    m_pho_pt.push_back( pho.pt() );
    m_pho_px.push_back( pho.px() );
    m_pho_py.push_back( pho.py() );
    m_pho_pz.push_back( pho.pz() );
    m_pho_eta.push_back( pho.eta() );
    m_pho_phi.push_back( pho.phi() );

    if ( !m_isData ) {
      m_pho_matchDR.push_back( tagger.matchDR()[i] );
      m_pho_tagged.push_back( tagger.tagResult()[i] );
      m_pho_matchPID.push_back( tagger.matchPID()[i] );
    }

  }
  m_pho_N = nPhotons;

  // get electron collection
  Handle<ElectronCollection> electrons;
  iEvent.getByLabel(m_Tag_Electrons,electrons);
  const unsigned nElectrons = electrons->size();

  // fill electron branches
  tagger.clearTags();
  if ( !m_isData && nElectrons ) tagger.tag(nElectrons,&(*electrons)[0]);
  for (unsigned i=0; i != nElectrons; i++) {

    const Electron & ele = (*electrons)[i];

    m_ele_E.push_back( ele.energy() );
    m_ele_p.push_back( ele.p() );
    m_ele_pt.push_back( ele.pt() );
    m_ele_px.push_back( ele.px() );
    m_ele_py.push_back( ele.py() );
    m_ele_pz.push_back( ele.pz() );
    m_ele_eta.push_back( ele.eta() );
    m_ele_phi.push_back( ele.phi() );

    if ( !m_isData ) {
      m_ele_matchDR.push_back( tagger.matchDR()[i] );
      m_ele_tagged.push_back( tagger.tagResult()[i] );
      m_ele_matchPID.push_back( tagger.matchPID()[i] );
    }

  }
  m_ele_N = nElectrons;

  // get MET collection
  Handle<std::vector<reco::PFMET> > met;
  iEvent.getByLabel(m_Tag_MET,met);

  // fill MET branches
  m_mpt = (*met)[0].sumEt();
  m_mpPhi = (*met)[0].phi();

  // fill tree, must go last in this function
  m_tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
MonoNtupleDumper::beginJob()
{
  m_outputFile = new TFile(m_output.c_str(), "recreate");

  //m_avgDir = new TFileDirectory( m_fs->mkdir("avgClusterMaps"));

  m_tree = new TTree("monopoles","Monopole Variables");

  m_tree->Branch("run",&m_run,"run/i");
  m_tree->Branch("lumiBlock",&m_lumi,"lumiBlock/i");
  m_tree->Branch("event",&m_event,"evnet/i");

  m_tree->Branch("NPV",&m_NPV,"NPV/i");

  m_tree->Branch("clust_N",&m_nClusters,"clust_N/i");
  m_tree->Branch("clust_E",&m_clust_E);
  m_tree->Branch("clust_eta",&m_clust_eta);
  m_tree->Branch("clust_phi",&m_clust_phi);
  m_tree->Branch("clust_L",&m_clust_L);
  m_tree->Branch("clust_W",&m_clust_W);
  m_tree->Branch("clust_sigEta",&m_clust_sigEta);
  m_tree->Branch("clust_sigPhi",&m_clust_sigPhi);
  m_tree->Branch("clust_skewEta",&m_clust_skewEta);
  m_tree->Branch("clust_skewPhi",&m_clust_skewPhi);
  m_tree->Branch("clust_seedFrac",&m_clust_seedFrac);
  m_tree->Branch("clust_firstFrac",&m_clust_firstFrac);
  m_tree->Branch("clust_secondFrac",&m_clust_secondFrac);
  m_tree->Branch("clust_thirdFrac",&m_clust_thirdFrac);
  m_tree->Branch("clust_matchDR",&m_clust_matchDR);
  m_tree->Branch("clust_matchTime",&m_clust_matchTime);
  m_tree->Branch("clust_matchPt",&m_clust_matchPt);
  m_tree->Branch("clust_matchPID",&m_clust_matchPID);
  m_tree->Branch("clust_tagged",&m_clust_tagged);
  m_tree->Branch("clust_hsE",&m_clust_hsE);
  m_tree->Branch("clust_hsTime",&m_clust_hsTime);
  m_tree->Branch("clust_hsInSeed",&m_clust_hsInSeed);
  m_tree->Branch("clust_hsWeird",&m_clust_hsWeird);
  m_tree->Branch("clust_hsDiWeird",&m_clust_hsDiWeird);
  m_tree->Branch("clust_Ecells",&m_clust_Ecells,"clust_Ecells[1500]/D");
  m_tree->Branch("clust_Tcells",&m_clust_Tcells,"clust_Tcells[1500]/D");

  m_tree->Branch("egClust_N",&m_nCleanEgamma,"egClust_N/i");
  m_tree->Branch("egClust_E",&m_egClust_E);
  m_tree->Branch("egClust_size",&m_egClust_size);
  m_tree->Branch("egClust_eta",&m_egClust_eta);
  m_tree->Branch("egClust_phi",&m_egClust_phi);
  m_tree->Branch("egClust_frac51",&m_egClust_frac51);
  m_tree->Branch("egClust_frac15",&m_egClust_frac15);
  m_tree->Branch("egClust_e55",&m_egClust_e55);
  m_tree->Branch("egClust_eMax",&m_egClust_eMax);
  m_tree->Branch("egClust_matchDR",&m_egClust_matchDR);
  m_tree->Branch("egClust_matchPID",&m_egClust_matchPID);
  m_tree->Branch("egClust_tagged",&m_egClust_tagged);

  m_tree->Branch("egClean_N",&m_nClusterEgamma,"egClean_N/i");
  m_tree->Branch("egClean_E",&m_egClean_E);
  m_tree->Branch("egClean_size",&m_egClean_size);
  m_tree->Branch("egClean_eta",&m_egClean_eta);
  m_tree->Branch("egClean_phi",&m_egClean_phi);
  m_tree->Branch("egClean_frac51",&m_egClean_frac51);
  m_tree->Branch("egClean_frac15",&m_egClean_frac15);
  m_tree->Branch("egClean_e55",&m_egClean_e55);
  m_tree->Branch("egClean_eMax",&m_egClean_eMax);
  m_tree->Branch("egClean_matchDR",&m_egClean_matchDR);
  m_tree->Branch("egClean_matchPID",&m_egClean_matchPID);
  m_tree->Branch("egClean_tagged",&m_egClean_tagged);

  m_tree->Branch("ehit_eta",&m_ehit_eta);
  m_tree->Branch("ehit_phi",&m_ehit_phi);
  m_tree->Branch("ehit_time",&m_ehit_time);
  m_tree->Branch("ehit_E",&m_ehit_energy);
  m_tree->Branch("ehit_kWeird",&m_ehit_kWeird);
  m_tree->Branch("ehit_kDiWeird",&m_ehit_kDiWeird);
  m_tree->Branch("ehit_flag",&m_ehit_flag);

  _Tracker->beginJob(m_tree);

  m_tree->Branch("jet_N",&m_jet_N,"jet_N/i");
  m_tree->Branch("jet_E",&m_jet_E);
  m_tree->Branch("jet_p",&m_jet_p);
  m_tree->Branch("jet_pt",&m_jet_pt);
  m_tree->Branch("jet_px",&m_jet_px);
  m_tree->Branch("jet_py",&m_jet_py);
  m_tree->Branch("jet_pz",&m_jet_pz);
  m_tree->Branch("jet_eta",&m_jet_eta);
  m_tree->Branch("jet_phi",&m_jet_phi);
  m_tree->Branch("jet_matchDR",&m_jet_matchDR);
  m_tree->Branch("jet_tagged",&m_jet_tagged);
  m_tree->Branch("jet_matchPID",&m_jet_matchPID);


  m_tree->Branch("pho_N",&m_pho_N,"pho_N/i");
  m_tree->Branch("pho_E",&m_pho_E);
  m_tree->Branch("pho_p",&m_pho_p);
  m_tree->Branch("pho_pt",&m_pho_pt);
  m_tree->Branch("pho_px",&m_pho_px);
  m_tree->Branch("pho_py",&m_pho_py);
  m_tree->Branch("pho_pz",&m_pho_pz);
  m_tree->Branch("pho_eta",&m_pho_eta);
  m_tree->Branch("pho_phi",&m_pho_phi);
  m_tree->Branch("pho_matchDR",&m_pho_matchDR);
  m_tree->Branch("pho_tagged",&m_pho_tagged);
  m_tree->Branch("pho_matchPID",&m_pho_matchPID);

  m_tree->Branch("ele_N",&m_ele_N,"ele_N/i");
  m_tree->Branch("ele_E",&m_ele_E);
  m_tree->Branch("ele_p",&m_ele_p);
  m_tree->Branch("ele_pt",&m_ele_pt);
  m_tree->Branch("ele_px",&m_ele_px);
  m_tree->Branch("ele_py",&m_ele_py);
  m_tree->Branch("ele_pz",&m_ele_pz);
  m_tree->Branch("ele_eta",&m_ele_eta);
  m_tree->Branch("ele_phi",&m_ele_phi);
  m_tree->Branch("ele_matchDR",&m_ele_matchDR);
  m_tree->Branch("ele_tagged",&m_ele_tagged);
  m_tree->Branch("ele_matchPID",&m_ele_matchPID);

  m_tree->Branch("mpt_pt",&m_mpt,"mpt_pt/D");
  m_tree->Branch("mpt_phi",&m_mpPhi,"mpt_phi/D");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MonoNtupleDumper::endJob() 
{
  m_outputFile->cd();
  m_tree->Write();

  for(std::map<Mono::ClustCategorizer,TH2D*>::iterator i = m_clustEMap.begin(); i != m_clustEMap.end(); i++)
  {
    i->second->Write();
  }
  for(std::map<Mono::ClustCategorizer,TH2D*>::iterator i = m_clustTMap.begin(); i != m_clustTMap.end(); i++)
  {
    i->second->Write();
  }

  _Tracker->endJob();

  m_outputFile->Close();
}

// ------------ method called when starting to processes a run  ------------
void 
MonoNtupleDumper::beginRun(edm::Run const&, edm::EventSetup const&)
{
}



void MonoNtupleDumper::clear()
{ 


     m_run = 0;
     m_lumi = 0;
     m_event = 0;
  
    m_NPV = 0;

    m_nClusters = 0;
    m_clust_E.clear();
    m_clust_eta.clear();
    m_clust_phi.clear();
    m_clust_L.clear();
    m_clust_W.clear();
    m_clust_sigEta.clear();
    m_clust_sigPhi.clear();
    m_clust_skewEta.clear();
    m_clust_skewPhi.clear();
    m_clust_seedFrac.clear();
    m_clust_firstFrac.clear();
    m_clust_secondFrac.clear();
    m_clust_thirdFrac.clear();
    m_clust_matchDR.clear();
    m_clust_matchTime.clear();
    m_clust_matchPt.clear();
    m_clust_matchPID.clear();
    m_clust_tagged.clear();
    m_clust_hsE.clear();
    m_clust_hsTime.clear();
    m_clust_hsInSeed.clear();
    m_clust_hsWeird.clear();
    m_clust_hsDiWeird.clear();
    for ( unsigned i=0; i != SS; i++ ){
      m_clust_Ecells[i] = -999.;
      m_clust_Tcells[i] = -999.;
    }


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

    m_nCleanEgamma = 0;
    m_egClean_E.clear();
    m_egClean_size.clear();
    m_egClean_eta.clear();
    m_egClean_phi.clear();
    m_egClean_frac51.clear();
    m_egClean_frac15.clear();
    m_egClean_e55.clear();
    m_egClean_eMax.clear();
    m_egClean_matchDR.clear();
    m_egClean_matchPID.clear();
    m_egClean_tagged.clear();

    // Ecal RecHits
    m_ehit_eta.clear();
    m_ehit_phi.clear();
    m_ehit_time.clear();
    m_ehit_energy.clear();
    m_ehit_otEnergy.clear();
    m_ehit_kWeird.clear();
    m_ehit_kDiWeird.clear();
    m_ehit_flag.clear();

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
    m_jet_matchDR.clear();
    m_jet_tagged.clear();
    m_jet_matchPID.clear();


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

    // Electron information
    m_ele_N = 0;
    m_ele_E.clear();
    m_ele_p.clear();
    m_ele_pt.clear();
    m_ele_px.clear();
    m_ele_py.clear();
    m_ele_pz.clear();
    m_ele_eta.clear();
    m_ele_phi.clear();

    m_mpt = 0.;
    m_mpPhi = 0.;

}



// ------------ method called when ending the processing of a run  ------------
void 
MonoNtupleDumper::endRun(edm::Run const&, edm::EventSetup const&)
{
  //m_ecalCalib.calculateHij();
  //m_ecalCalib.dumpCalibration();

  std::map<Mono::ClustCategorizer,unsigned>::iterator counts = m_clustCatCount.begin();
  std::map<Mono::ClustCategorizer,unsigned>::iterator end = m_clustCatCount.end();
  for ( ; counts != end; counts++ ) {
    const Mono::ClustCategorizer & cat = counts->first;
    unsigned & count = counts->second;
    m_clustEMap[cat]->Scale(1.0/count);  
    m_clustTMap[cat]->Scale(1.0/count);
  }

}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MonoNtupleDumper::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MonoNtupleDumper::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MonoNtupleDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MonoNtupleDumper);
