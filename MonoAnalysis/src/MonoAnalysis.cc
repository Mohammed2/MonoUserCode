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
// $Id: MonoAnalysis.cc,v 1.2 2012/06/14 17:08:14 cowden Exp $
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


// ROOT includes
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF2.h"


//
// class declaration
//

class MonoAnalysis : public edm::EDAnalyzer {

   typedef std::vector<reco::Photon> PhotonCollection;
   //typedef std::vector<reco::Electron> ElectronCollection;
   typedef std::vector<reco::GsfElectron> ElectronCollection;
   typedef std::vector<reco::BasicCluster> BasicClusterCollection;


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
    edm::InputTag m_Tag_Electrons;

    // Monopole Ecal Observables
    Mono::MonoEcalObs0 m_ecalObs;


    // TFileService
    edm::Service<TFileService> m_fs;

    TF2 * m_func;
    double m_fitParams[5];
    double m_fitLimits[4][2];

    Mono::EnergyFlowFunctor m_functor;

    TTree * m_tree;

    // Event information
    unsigned m_run;
    unsigned m_lumi;
    unsigned m_event;

    unsigned m_NPV;

    // Ecal Observable information
    unsigned m_nSeeds;
    std::vector<double> m_seed_E;
    std::vector<double> m_seed_eta;
    std::vector<double> m_seed_phi;
    std::vector<double> m_seed_ieta;
    std::vector<double> m_seed_iphi;
    std::vector<double> m_seed_L;
    std::vector<double> m_seed_cell_E;
    std::vector<double> m_seed_cell_time;
    std::vector<std::vector<double> > m_seed_cell_eDist;
    std::vector<std::vector<double> > m_seed_cell_tDist;

    unsigned m_nClusters;
    std::vector<double> m_clust_E;
    std::vector<double> m_clust_L;
    std::vector<double> m_clust_W;
    std::vector<double> m_clust_N;
    std::vector<double> m_clust_sigEta;
    std::vector<double> m_clust_sigPhi;
    std::vector<double> m_clust_meanEta;
    std::vector<double> m_clust_meanPhi;
    std::vector<double> m_clust_chi2;
    std::vector<double> m_clust_NDF;
    std::vector<double> m_clust_diff;

    // Ecal hybrid clusters
    unsigned m_nClusterEgamma;
    std::vector<double> m_egClust_E;
    std::vector<double> m_egClust_size;
    std::vector<double> m_egClust_eta;
    std::vector<double> m_egClust_phi;
    


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
  ,m_Tag_Electrons(iConfig.getParameter<edm::InputTag>("ElectronTag") )
  ,m_ecalObs(iConfig)
{
   //now do what ever initialization is needed
   m_seed_cell_eDist.resize(10);
   m_seed_cell_tDist.resize(10);

   m_func = new TF2("myFunc","[0]*exp(-([1]-x)^2/[2]^2/2-([3]-y)^2/[4]^2/2)",-10,10,-10,10);
   m_func->SetParName(0,"N");
   m_func->SetParName(1,"Mean X");
   m_func->SetParName(2,"Sigma X");
   m_func->SetParName(3,"Mean Y");
   m_func->SetParName(4,"Sigma Y");

   m_fitParams[0] = 1.;
   m_fitParams[1] = 0.;
   m_fitParams[2] = 1.;
   m_fitParams[3] = 0.;
   m_fitParams[4] = 1.;  

   m_fitLimits[0][0] = -2.;
   m_fitLimits[0][1] = 2.;
   m_fitLimits[1][0] = 0.;
   m_fitLimits[1][1] = 6.;
   m_fitLimits[2][0] = -2.;
   m_fitLimits[2][1] = 2.;
   m_fitLimits[3][0] = 0.;
   m_fitLimits[3][1] = 6.;


   double pars[3] = {1.0,0.4,0.4};
   m_functor.setParameters(3,pars);  

}


MonoAnalysis::~MonoAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

   delete m_func;
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

  // execute observable calculations
  double monoObs = m_ecalObs.calculate(iSetup,iEvent);
  const Mono::StripSeedFinder & sFinder = m_ecalObs.finder();
  const Mono::EBmap & ebMap = m_ecalObs.ecalMap();

  //const unsigned seedLength = sFinder.seedLength();
  const unsigned nEta = ebMap.nEta();
  const unsigned nPhi = ebMap.nPhi();
  m_nSeeds = sFinder.nSeeds();
  for ( unsigned i=0; i != m_nSeeds; i++ ) {
    const unsigned seedLength = sFinder.seeds()[i].seedLength();
    m_seed_L.push_back( seedLength );
    m_seed_E.push_back( sFinder.seeds()[i].energy() );

    const unsigned iphi = sFinder.seeds()[i].iphi();
    const unsigned ieta = sFinder.seeds()[i].ieta();
    m_seed_ieta.push_back( ieta );
    m_seed_iphi.push_back( iphi );

    m_seed_eta.push_back( ebMap.eta(ieta) );
    m_seed_phi.push_back( ebMap.phi(iphi) );

    char histName[50];
    sprintf(histName,"seedHist_%d_%d",iEvent.id().event(),i);
    TH1D *hist = m_fs->make<TH1D>(histName,histName,seedLength,0,seedLength);
    sprintf(histName,"seedTHist_%d_%d",iEvent.id().event(),i);
    TH1D *Thist = m_fs->make<TH1D>(histName,histName,seedLength,0,seedLength);

    for (unsigned c=0; c != seedLength; c++ ) {
      //const unsigned loc = (iphi+c)*nEta+ieta;  // cross-check
      const unsigned loc = iphi*nEta+ieta+c;
      m_seed_cell_E.push_back( ebMap[loc] );
      m_seed_cell_time.push_back( ebMap.time(loc) );
      m_seed_cell_eDist[c].push_back( ebMap[loc] );
      m_seed_cell_tDist[c].push_back( ebMap.time(loc) );
      hist->SetBinContent(c+1,ebMap[loc]);
      Thist->SetBinContent(c+1,ebMap.time(loc));
    }
  }

  const Mono::ClusterBuilder clusterBuilder = m_ecalObs.clusterBuilder();
  m_nClusters = clusterBuilder.nClusters();
  for ( unsigned i=0; i != m_nClusters; i++ ) {
    const unsigned width = clusterBuilder.clusters()[i].clusterWidth();
    const unsigned length = clusterBuilder.clusters()[i].clusterLength();
    const unsigned cEta = clusterBuilder.clusters()[i].ieta();
    const unsigned cPhi = clusterBuilder.clusters()[i].iphi();

    m_clust_E.push_back( clusterBuilder.clusters()[i].clusterEnergy() );
    m_clust_L.push_back( length );
    m_clust_W.push_back( width );

    const unsigned wings = width/2U;
    char histName[50];
    sprintf(histName,"clustHist_%d_%d",iEvent.id().event(),i);
    TH2D *hist = m_fs->make<TH2D>(histName,histName,length,-(float)length/2.,(float)length/2.,width,-(int)wings,wings);
    sprintf(histName,"clustTHist_%d_%d",iEvent.id().event(),i);
    TH2D *Thist = m_fs->make<TH2D>(histName,histName,length,0,length,width,-(int)wings,wings);

    hist->GetXaxis()->SetTitle("#eta bin"); 
    hist->GetYaxis()->SetTitle("#phi bin"); 
    hist->GetZaxis()->SetTitle("energy"); 
    Thist->GetXaxis()->SetTitle("#eta bin"); 
    Thist->GetYaxis()->SetTitle("#phi bin"); 
    Thist->GetZaxis()->SetTitle("time"); 
    for ( unsigned j=wings; j != 0; j-- ) { 
      int newPhiS = cPhi-j; 
      unsigned newPhi = 0U; 
      if ( newPhiS < 0 ) newPhi = nPhi+newPhiS; 
      else newPhi = newPhiS; 
      assert( newPhi < nPhi ); 
      for ( unsigned k=0; k != length; k++ ) { 
	//hist->SetBinContent(k+1,wings+1-j,ebMap[newPhi*nEta+cEta+k]); 
        hist->SetBinContent(k+1,wings+1-j,clusterBuilder.clusters()[i].energy((int)k,newPhiS,ebMap));
	Thist->SetBinContent(k+1,j-wings+1,ebMap.time(newPhi*nEta+cEta+k)); 
      } 
    } 
    for ( unsigned j=0; j != wings+1; j++ ) { 
      unsigned newPhi = cPhi+j;
      if ( newPhi >= nPhi ) newPhi = newPhi-nPhi;
      assert( newPhi < nPhi );
      for ( unsigned k=0; k != length; k++ ) {
	//hist->SetBinContent(k+1,j+wings+1,ebMap[newPhi*nEta+cEta+k]); 
	hist->SetBinContent(k+1,j+wings+1,clusterBuilder.clusters()[i].energy((int)k,(int)newPhi,ebMap)); 
	Thist->SetBinContent(k+1,j+wings+1,ebMap.time(newPhi*nEta+cEta+k)); 
      }
    }

    // perform Gaussian fit to cluster
    // normalise to total energy of cluster
    hist->Scale( 1./clusterBuilder.clusters()[i].clusterEnergy() );
    m_func->SetParameters(m_fitParams);
    for ( unsigned i=0; i != 4; i++ )
      m_func->SetParLimits(i+1,m_fitLimits[i][0],m_fitLimits[i][1]);
    hist->Fit("myFunc","QB");
    const TF1 * theFit = hist->GetFunction("myFunc");
    //m_clust_N.push_back( theFit->GetParameter(0) );
    m_clust_meanEta.push_back( theFit->GetParameter(1) );
    m_clust_sigEta.push_back( theFit->GetParameter(2) );
    m_clust_meanPhi.push_back( theFit->GetParameter(3) );
    m_clust_sigPhi.push_back( theFit->GetParameter(4) );
    m_clust_chi2.push_back( theFit->GetChisquare() ); 
    m_clust_NDF.push_back( theFit->GetNDF() ); 
  }


  // get BasicCluster Collection
  Handle<BasicClusterCollection> bClusters;
  edm::InputTag bcClusterTag("hybridSuperClusters","hybridBarrelBasicClusters"); 
  iEvent.getByLabel(bcClusterTag,bClusters);
  const unsigned nbClusters = bClusters->size();
  for ( unsigned i=0; i != nbClusters; i++ ) {
    m_egClust_E.push_back( (*bClusters)[i].energy() );
    m_egClust_size.push_back( (*bClusters)[i].size() );
    m_egClust_eta.push_back( (*bClusters)[i].eta() );
    m_egClust_phi.push_back( (*bClusters)[i].phi() );
  }
  m_nClusterEgamma = nbClusters;


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

  // get electron collection
  Handle<ElectronCollection> electrons;
  iEvent.getByLabel(m_Tag_Electrons,electrons);

  // fill electron branches
  ElectronCollection::const_iterator itEle = electrons->begin();
  for ( ; itEle != electrons->end(); itEle++ ) {

    m_ele_E.push_back( (*itEle).energy() );
    m_ele_p.push_back( (*itEle).p() );
    m_ele_pt.push_back( (*itEle).pt() );
    m_ele_px.push_back( (*itEle).px() );
    m_ele_py.push_back( (*itEle).py() );
    m_ele_pz.push_back( (*itEle).pz() );
    m_ele_eta.push_back( (*itEle).eta() );
    m_ele_phi.push_back( (*itEle).phi() );

    m_ele_N++;

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

  m_tree->Branch("seed_N",&m_nSeeds,"seed_N/i");
  m_tree->Branch("seed_E",&m_seed_E);
  m_tree->Branch("seed_eta",&m_seed_eta);
  m_tree->Branch("seed_phi",&m_seed_phi);
  m_tree->Branch("seed_ieta",&m_seed_ieta);
  m_tree->Branch("seed_iphi",&m_seed_iphi);
  m_tree->Branch("seed_L",&m_seed_L);
  m_tree->Branch("seedCell_E",&m_seed_cell_E);
  m_tree->Branch("seedCell_time",&m_seed_cell_time);
  
  const unsigned nCells = m_seed_cell_eDist.size();
  char name[50];
  for ( unsigned i=0; i != nCells; i++ ) {
    sprintf(name,"seedCell_E_%d",i);
    m_tree->Branch(name,&m_seed_cell_eDist[i]);
    sprintf(name,"seedCell_time_%d",i);
    m_tree->Branch(name,&m_seed_cell_tDist[i]);
  }

  m_tree->Branch("clust_N",&m_nClusters,"clust_N/i");
  m_tree->Branch("clust_E",&m_clust_E);
  m_tree->Branch("clust_L",&m_clust_L);
  m_tree->Branch("clust_W",&m_clust_W);
  m_tree->Branch("clust_fitN",&m_clust_N);
  m_tree->Branch("clust_sigEta",&m_clust_sigEta);
  m_tree->Branch("clust_sigPhi",&m_clust_sigPhi);
  m_tree->Branch("clust_meanEta",&m_clust_meanEta);
  m_tree->Branch("clust_meanPhi",&m_clust_meanPhi);
  m_tree->Branch("clust_chi2",&m_clust_chi2);
  m_tree->Branch("clust_NDF",&m_clust_NDF);

  m_tree->Branch("egClust_N",&m_nClusterEgamma,"egClust_N/i");
  m_tree->Branch("egClust_E",&m_egClust_E);
  m_tree->Branch("egClust_size",&m_egClust_size);
  m_tree->Branch("egClust_eta",&m_egClust_eta);
  m_tree->Branch("egClust_phi",&m_egClust_phi);

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

  m_tree->Branch("ele_N",&m_ele_N,"ele_N/i");
  m_tree->Branch("ele_E",&m_ele_E);
  m_tree->Branch("ele_p",&m_ele_p);
  m_tree->Branch("ele_pt",&m_ele_pt);
  m_tree->Branch("ele_px",&m_ele_px);
  m_tree->Branch("ele_py",&m_ele_py);
  m_tree->Branch("ele_pz",&m_ele_pz);
  m_tree->Branch("ele_eta",&m_ele_eta);
  m_tree->Branch("ele_phi",&m_ele_phi);

  m_tree->Branch("ehit_eta",&m_ehit_eta);
  m_tree->Branch("ehit_phi",&m_ehit_phi);
  m_tree->Branch("ehit_time",&m_ehit_time);


}



void MonoAnalysis::clear()
{ 


     m_run = 0;
     m_lumi = 0;
     m_event = 0;
  
    m_NPV = 0;

    // obs information
    m_nSeeds = 0;
    m_seed_E.clear();
    m_seed_eta.clear();
    m_seed_phi.clear();
    m_seed_ieta.clear();
    m_seed_iphi.clear();
    m_seed_L.clear();
    m_seed_cell_E.clear();
    m_seed_cell_time.clear();
   
    const unsigned nCells = m_seed_cell_eDist.size();
    for ( unsigned i=0; i != nCells; i++ ) {
      m_seed_cell_eDist[i].clear();
      m_seed_cell_tDist[i].clear();
    }

    m_nClusters = 0;
    m_clust_E.clear();
    m_clust_L.clear();
    m_clust_W.clear();
    m_clust_N.clear();
    m_clust_sigEta.clear();
    m_clust_sigPhi.clear();
    m_clust_meanEta.clear();
    m_clust_meanPhi.clear();
    m_clust_chi2.clear();
    m_clust_NDF.clear();


    m_nClusterEgamma = 0;
    m_egClust_E.clear();
    m_egClust_size.clear();
    m_egClust_eta.clear();
    m_egClust_phi.clear();

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
