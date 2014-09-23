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
// $Id: MonoAnalysis.cc,v 1.11 2013/03/17 12:44:31 cowden Exp $
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
    bool m_isData;

    // Monopole Ecal Observables
    Mono::MonoEcalObs0 m_ecalObs;
    //Mono::MonoEcalObs0Calibrator m_ecalCalib;


    // TFileService
    edm::Service<TFileService> m_fs;

    TF2 * m_func;
    double m_fitParams[5];
    double m_fitLimits[4][2];

    Mono::EnergyFlowFunctor m_functor;

    // directory for cluster map histograms
    TFileDirectory *m_histDir;
    // map cluster category (lengthxwidth thing) to a histogram
    // showing the average energy or time  in each cell
    TFileDirectory *m_avgDir;
    std::map<Mono::ClustCategorizer,TH2D *> m_clustEMap;
    std::map<Mono::ClustCategorizer,TH2D *> m_clustTMap;
    std::map<Mono::ClustCategorizer,unsigned> m_clustCatCount;

    TTree * m_tree;

    // Event information
    unsigned m_run;
    unsigned m_lumi;
    unsigned m_event;

    unsigned m_NPV;

    // observables
    std::vector<double> m_betas;
    std::vector<double> m_betaTs;

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
    std::vector<double> m_clust_eta;
    std::vector<double> m_clust_phi;
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

    // MC gen monopoles
    double m_mono_eta;
    double m_mono_phi;
    double m_mono_pT;
    double m_mono_extrapEta;
    double m_amon_eta;
    double m_amon_phi;
    double m_amon_pT;
    double m_amon_extrapEta;



    // Ecal hybrid clusters
    unsigned m_nClusterEgamma;
    std::vector<double> m_egClust_E;
    std::vector<double> m_egClust_size;
    std::vector<double> m_egClust_eta;
    std::vector<double> m_egClust_phi;
    std::vector<double> m_egClust_matchDR;
    std::vector<double> m_egClust_tagged;
    std::vector<double> m_egClust_matchPID;


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
    std::vector<double> m_ehit_flag;
    std::vector<double> m_ehit_kWeird;
    std::vector<double> m_ehit_kDiWeird;
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
  ,m_isData(iConfig.getParameter<bool>("isData") )
  ,m_ecalObs(iConfig)
  //,m_ecalCalib(iConfig)
{
   if(CLHEP::electron_charge==0) std::cout << "asdf" << std::endl;

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

  // execute calibration for event
  //m_ecalCalib.calculateMijn(iSetup,iEvent);

  // execute observable calculations
  //double monoObs = m_ecalObs.calculate(iSetup,iEvent,&m_betas,&m_betaTs);
  const Mono::StripSeedFinder & sFinder = m_ecalObs.finder();
  const Mono::EBmap & ebMap = m_ecalObs.ecalMap();

  //const unsigned seedLength = sFinder.seedLength();
  const unsigned nEta = ebMap.nEta();
  //const unsigned nPhi = ebMap.nPhi();
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
    TH1D *hist = m_histDir->make<TH1D>(histName,histName,seedLength,0,seedLength);
    sprintf(histName,"seedTHist_%d_%d",iEvent.id().event(),i);
    TH1D *Thist = m_histDir->make<TH1D>(histName,histName,seedLength,0,seedLength);

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


  /////////////////////////////////////
  // cluster analysis
  const Mono::ClusterBuilder clusterBuilder = m_ecalObs.clusterBuilder();
  m_nClusters = clusterBuilder.nClusters();

  Mono::GenMonoClusterTagger tagger(0.3);
  if ( !m_isData ) {
    tagger.initialize(iEvent,iSetup);
    if ( m_nClusters ) tagger.tag(m_nClusters,clusterBuilder.clusters(),ebMap);
  }

  for ( unsigned i=0; i != m_nClusters; i++ ) {
    const Mono::MonoEcalCluster & cluster = clusterBuilder.clusters()[i];
    const unsigned width = cluster.clusterWidth();
    const unsigned length = cluster.clusterLength();
    const unsigned cEta = cluster.ieta();
    const unsigned cPhi = cluster.iphi();

    bool isTagged = false;

    if ( !m_isData ) {
      m_clust_matchDR.push_back(tagger.matchDR()[i]);
      m_clust_tagged.push_back(tagger.tagResult()[i]);
      if ( tagger.tagResult()[i] ) isTagged = true;
      m_clust_matchPID.push_back(tagger.matchPID()[i]);
      m_clust_matchTime.push_back(tagger.matchTime()[i]);
      m_clust_matchPt.push_back(tagger.matchPt()[i]);
    }

    m_clust_E.push_back( cluster.clusterEnergy() );
    m_clust_eta.push_back( ebMap.eta(cEta) );
    m_clust_phi.push_back( ebMap.phi(cPhi) );
    m_clust_L.push_back( length );
    m_clust_W.push_back( width );


    const unsigned wings = width/2U;
    char histName[50];
    sprintf(histName,"clustHist_%d_%d",iEvent.id().event(),i);
    TH2D *hist = m_histDir->make<TH2D>(histName,histName,length,-(float)length/2.,(float)length/2.,width,-(int)wings,wings);
    sprintf(histName,"clustTHist_%d_%d",iEvent.id().event(),i);
    TH2D *Thist = m_histDir->make<TH2D>(histName,histName,length,0,length,width,-(int)wings,wings);

    hist->GetXaxis()->SetTitle("#eta bin"); 
    hist->GetYaxis()->SetTitle("#phi bin"); 
    hist->GetZaxis()->SetTitle("energy"); 
    Thist->GetXaxis()->SetTitle("#eta bin"); 
    Thist->GetYaxis()->SetTitle("#phi bin"); 
    Thist->GetZaxis()->SetTitle("time"); 

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
	hist->SetBinContent(k+1,j+1,energy);
	Thist->SetBinContent(k+1,j+1,cluster.time(k,ji,ebMap));
	if ( !excessive ) m_clust_Ecells[i*WS+j*length+k] = cluster.energy(k,ji,ebMap);
	if ( !excessive ) m_clust_Tcells[i*WS+j*length+k] = cluster.time(k,ji,ebMap);
      }
    }


    m_clust_sigEta.push_back( hist->GetRMS(1) );
    m_clust_sigPhi.push_back( hist->GetRMS(2) );
    m_clust_skewEta.push_back( hist->GetSkewness(1) );
    m_clust_skewPhi.push_back( hist->GetSkewness(2) );

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

    char text[50];
    sprintf(text,"Cluster Energy beta=%.4f, T=%d",0.,isTagged);
    hist->SetTitle(text); 
    // perform Gaussian fit to cluster
    // normalise to total energy of cluster
    hist->Scale( 1./cluster.clusterEnergy() );
    m_func->SetParameters(m_fitParams);
    for ( unsigned j=0; j != 4; j++ )
      m_func->SetParLimits(j+1,m_fitLimits[j][0],m_fitLimits[j][1]);
    //hist->Fit("myFunc","QB");
    //const TF1 * theFit = hist->GetFunction("myFunc");
    //m_clust_N.push_back( theFit->GetParameter(0) );
    //m_clust_meanEta.push_back( theFit->GetParameter(1) );
    //m_clust_sigEta.push_back( theFit->GetParameter(2) );
    //m_clust_meanPhi.push_back( theFit->GetParameter(3) );
    //m_clust_sigPhi.push_back( theFit->GetParameter(4) );
    //m_clust_chi2.push_back( theFit->GetChisquare() ); 
    //m_clust_NDF.push_back( theFit->GetNDF() ); 


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
      m_clustEMap[cluCat] = m_avgDir->make<TH2D>(name,name,cluCat.length,0,cluCat.length,cluCat.width,-(int)cluCat.width/2,(int)cluCat.width/2);
      sprintf(name,"avgTmClust_%d_%d",cluCat.length,cluCat.width);
      m_clustTMap[cluCat] = m_avgDir->make<TH2D>(name,name,cluCat.length,0,cluCat.length,cluCat.width,-(int)cluCat.width/2,(int)cluCat.width/2);
    }

    m_clustCatCount[cluCat]++;
    TH2D * avgEnMap = m_clustEMap[cluCat];
    TH2D * avgTmMap = m_clustTMap[cluCat];

    assert( avgEnMap );
    assert( avgTmMap );

    if ( isTagged ) 
    for ( int binx = 1; binx <= hist->GetNbinsX(); binx++ ) {
      for ( int biny = 1; biny <= hist->GetNbinsY(); biny++ ) {
      	avgEnMap->SetBinContent(binx,biny,avgEnMap->GetBinContent(binx,biny)+hist->GetBinContent(binx,biny));
      	avgTmMap->SetBinContent(binx,biny,avgTmMap->GetBinContent(binx,biny)+Thist->GetBinContent(binx,biny));
      }
    }
    

  }


  ////////////////////////
  // MC truth collection
  if ( !m_isData ) {
    Mono::MonoTruthSnoop snoopy(iEvent,iSetup);
    const HepMC::GenParticle *mono = snoopy.mono(Mono::monopole);
    const HepMC::GenParticle *amon = snoopy.mono(Mono::anti_monopole);

    Mono::MonoGenTrackExtrapolator extrap;

    if ( mono ) {
      m_mono_eta = mono->momentum().eta();
      m_mono_phi = mono->momentum().phi();
      m_mono_pT = mono->momentum().perp();

      extrap.setMonopole(*mono);
      m_mono_extrapEta = extrap.etaVr(1.29);
    }

    if ( amon ) {
      m_amon_eta = amon->momentum().eta();
      m_amon_phi = amon->momentum().phi();
      m_amon_pT = amon->momentum().perp();

      extrap.setMonopole(*amon);
      m_amon_extrapEta = extrap.etaVr(1.29);
    }
  }



  // get BasicCluster Collection
  Handle<BasicClusterCollection> bClusters;
  edm::InputTag bcClusterTag("hybridSuperClusters","hybridBarrelBasicClusters"); 
  iEvent.getByLabel(bcClusterTag,bClusters);
  const unsigned nbClusters = bClusters->size();

  tagger.clearTags();
  if ( !m_isData && nbClusters ) tagger.tag(nbClusters,&(*bClusters)[0]);

  for ( unsigned i=0; i != nbClusters; i++ ) {
    m_egClust_E.push_back( (*bClusters)[i].energy() );
    m_egClust_size.push_back( (*bClusters)[i].size() );
    m_egClust_eta.push_back( (*bClusters)[i].eta() );
    m_egClust_phi.push_back( (*bClusters)[i].phi() );
    if ( !m_isData ) {
      m_egClust_matchDR.push_back(tagger.matchDR()[i]);
      m_egClust_tagged.push_back(tagger.tagResult()[i]);
      m_egClust_matchPID.push_back(tagger.matchPID()[i]);
    }
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

    m_ehit_kWeird.push_back( (*itHit).checkFlag(EcalRecHit::kWeird) );
    m_ehit_kDiWeird.push_back( (*itHit).checkFlag(EcalRecHit::kDiWeird) );
    m_ehit_flag.push_back( (*itHit).recoFlag() );

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

  m_histDir = new TFileDirectory( m_fs->mkdir("clusterMaps"));
  m_avgDir = new TFileDirectory( m_fs->mkdir("avgClusterMaps"));

  m_tree = m_fs->make<TTree>("tree","tree");

  m_tree->Branch("run",&m_run,"run/i");
  m_tree->Branch("lumiBlock",&m_lumi,"lumiBlock/i");
  m_tree->Branch("event",&m_event,"evnet/i");

  m_tree->Branch("NPV",&m_NPV,"NPV/i");

  m_tree->Branch("betas_E",&m_betas);

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

  m_tree->Branch("mono_eta",&m_mono_eta,"mono_eta/D");
  m_tree->Branch("mono_phi",&m_mono_phi,"mono_phi/D");
  m_tree->Branch("mono_pT",&m_mono_pT,"mono_pT/D");
  m_tree->Branch("mono_extrapEta",&m_mono_extrapEta,"mono_extrapEta/D");
  
  m_tree->Branch("amon_eta",&m_amon_eta,"amon_eta/D");
  m_tree->Branch("amon_phi",&m_amon_phi,"amon_phi/D");
  m_tree->Branch("amon_pT",&m_amon_pT,"amon_pT/D");
  m_tree->Branch("amon_extrapEta",&m_amon_extrapEta,"amon_extrapEta/D");

  m_tree->Branch("egClust_N",&m_nClusterEgamma,"egClust_N/i");
  m_tree->Branch("egClust_E",&m_egClust_E);
  m_tree->Branch("egClust_size",&m_egClust_size);
  m_tree->Branch("egClust_eta",&m_egClust_eta);
  m_tree->Branch("egClust_phi",&m_egClust_phi);
  m_tree->Branch("egClust_matchDR",&m_egClust_matchDR);
  m_tree->Branch("egClust_matchPID",&m_egClust_matchPID);
  m_tree->Branch("egClust_tagged",&m_egClust_tagged);

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
  m_tree->Branch("ehit_E",&m_ehit_energy);
  m_tree->Branch("ehit_kWeird",&m_ehit_kWeird);
  m_tree->Branch("ehit_kDiWeird",&m_ehit_kDiWeird);
  m_tree->Branch("ehit_flag",&m_ehit_flag);


}



void MonoAnalysis::clear()
{ 


     m_run = 0;
     m_lumi = 0;
     m_event = 0;
  
    m_NPV = 0;

    m_betas.clear();

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

    m_mono_eta = 0.;
    m_mono_phi = 0.;
    m_mono_pT = 0.;
    m_mono_extrapEta = 0.;
    m_amon_eta = 0.;
    m_amon_phi = 0.;
    m_amon_pT = 0.;
    m_amon_extrapEta = 0.;


    m_nClusterEgamma = 0;
    m_egClust_E.clear();
    m_egClust_size.clear();
    m_egClust_eta.clear();
    m_egClust_phi.clear();
    m_egClust_matchDR.clear();
    m_egClust_matchPID.clear();
    m_egClust_tagged.clear();


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
    m_ehit_kWeird.clear();
    m_ehit_kDiWeird.clear();
    m_ehit_flag.clear();
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
