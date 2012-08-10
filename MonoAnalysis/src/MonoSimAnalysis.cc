// -*- C++ -*-
//
// Package:    MonoSimAnalysis
// Class:      MonoSimAnalysis
// 
/**\class MonoSimAnalysis MonoSimAnalysis.cc Monopoles/MonoAnalysis/src/MonoSimAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christopher Cowden
//         Created:  Tue Feb  7 16:21:08 CST 2012
// $Id: MonoSimAnalysis.cc,v 1.2 2012/06/14 17:08:14 cowden Exp $
//
//


// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


//Data Formats
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"


// Monopole includes
#include "Monopoles/MonoAlgorithms/interface/MonoDefs.h"
#include "Monopoles/MonoAlgorithms/interface/MonoSimTracker.h"
#include "Monopoles/MonoAlgorithms/interface/MonoTruthSnooper.h"


// ROOT include files
#include "TLorentzVector.h"
#include "TH2D.h"
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TLatex.h"

//
// class declaration
//

class MonoSimAnalysis : public edm::EDAnalyzer {
   public:
      explicit MonoSimAnalysis(const edm::ParameterSet&);
      ~MonoSimAnalysis();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

    // clear TBranche variables
    void clear();

    // do analysis of monopole Ecal strands
    void strandAnalysis();



      // ----------member data ---------------------------

    // TFileService
    edm::Service<TFileService> m_fs;

    // edm InputTags to retrieve simulation collections
    edm::InputTag m_TagSimTracks;
    edm::InputTag m_TagSimHits;
    edm::InputTag m_TagCaloHits;


    // TTree
    TTree * m_tree;

    // TTree TBranches
    unsigned m_mono_Ecal_N;
    std::vector<double> m_mono_Ecal_x;
    std::vector<double> m_mono_Ecal_y;
    std::vector<double> m_mono_Ecal_z;
    std::vector<double> m_mono_Ecal_rho;
    std::vector<double> m_mono_Ecal_eta;
    std::vector<double> m_mono_Ecal_phi;
    std::vector<double> m_mono_Ecal_time;
    std::vector<double> m_mono_Ecal_energy;
    std::vector<unsigned> m_mono_Ecal_id;

    unsigned m_amon_Ecal_N;
    std::vector<double> m_amon_Ecal_x;
    std::vector<double> m_amon_Ecal_y;
    std::vector<double> m_amon_Ecal_z;
    std::vector<double> m_amon_Ecal_rho;
    std::vector<double> m_amon_Ecal_eta;
    std::vector<double> m_amon_Ecal_phi;
    std::vector<double> m_amon_Ecal_time;
    std::vector<double> m_amon_Ecal_energy;


    unsigned m_mono_EcalSum_Nids;
    std::vector<unsigned> m_mono_EcalSum_id;
    std::vector<unsigned> m_mono_EcalSum_N;
    std::vector<double>   m_mono_EcalSum_energy;
    std::vector<double>   m_mono_EcalSum_eta;
    std::vector<double>   m_mono_EcalSum_phi;

    unsigned m_amon_EcalSum_Nids;
    std::vector<unsigned> m_amon_EcalSum_id;
    std::vector<unsigned> m_amon_EcalSum_N;
    std::vector<double>   m_amon_EcalSum_energy;
    std::vector<double>   m_amon_EcalSum_eta;
    std::vector<double>   m_amon_EcalSum_phi;
    


    unsigned m_mono_Pix_N;
    std::vector<double> m_mono_Pix_x;
    std::vector<double> m_mono_Pix_y;
    std::vector<double> m_mono_Pix_z;
    std::vector<double> m_mono_Pix_rho;
    std::vector<double> m_mono_Pix_eta;
    std::vector<double> m_mono_Pix_phi;
    std::vector<double> m_mono_Pix_tof;
    std::vector<double> m_mono_Pix_energy;
    std::vector<double> m_mono_Pix_length;

    unsigned m_amon_Pix_N;
    std::vector<double> m_amon_Pix_x;
    std::vector<double> m_amon_Pix_y;
    std::vector<double> m_amon_Pix_z;
    std::vector<double> m_amon_Pix_rho;
    std::vector<double> m_amon_Pix_eta;
    std::vector<double> m_amon_Pix_phi;
    std::vector<double> m_amon_Pix_tof;
    std::vector<double> m_amon_Pix_energy;
    std::vector<double> m_amon_Pix_length;


    unsigned m_mono_PixSum_Nids;
    std::vector<unsigned> m_mono_PixSum_id;
    std::vector<unsigned> m_mono_PixSum_N;
    std::vector<double>   m_mono_PixSum_energy;
    std::vector<double>   m_mono_PixSum_eta;
    std::vector<double>   m_mono_PixSum_phi;

    unsigned m_amon_PixSum_Nids;
    std::vector<unsigned> m_amon_PixSum_id;
    std::vector<unsigned> m_amon_PixSum_N;
    std::vector<double>   m_amon_PixSum_energy;
    std::vector<double>   m_amon_PixSum_eta;
    std::vector<double>   m_amon_PixSum_phi;



    // Generator level branches
    double m_mono_p;
    double m_mono_eta;
    double m_mono_phi;
    double m_mono_m;

    double m_amon_p;
    double m_amon_eta;
    double m_amon_phi;
    double m_amon_m;



};

//
// constants, enums and typedefs
//

namespace cow {

double mag ( double x, double y) {
  return sqrt( x*x + y*y );
}

double mag ( double x, double y, double z){
  return sqrt( x*x + y*y + z*z );
}


}


//
// static data member definitions
//

//
// constructors and destructor
//
MonoSimAnalysis::MonoSimAnalysis(const edm::ParameterSet& iConfig)
  :m_TagSimTracks(iConfig.getParameter<edm::InputTag>("SimTracksTag") )
  ,m_TagSimHits(iConfig.getParameter<edm::InputTag>("SimHitsTag") )
  ,m_TagCaloHits(iConfig.getParameter<edm::InputTag>("SimCaloTag") )
{
   //now do what ever initialization is needed

}


MonoSimAnalysis::~MonoSimAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MonoSimAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  clear();


  Mono::MonoSimTracker<PSimHit,TrackerGeometry> monoPixST(iEvent,iSetup,Mono::PixelEBLowTof);
  Mono::MonoSimTracker<PCaloHit,CaloGeometry> monoEcalST(iEvent,iSetup,Mono::EcalEB);

  Mono::MonoEnum m = Mono::monopole;
  Mono::MonoEnum a = Mono::anti_monopole;

  // fill monopole Pixel SimHit vectors
  for ( unsigned i=0; i != monoPixST.size(Mono::monopole); i++ ) {

    const double x = monoPixST.x(m,i);
    const double y = monoPixST.y(m,i);
    const double z = monoPixST.z(m,i);

    const double rho = sqrt(x*x + y*y);

    m_mono_Pix_x.push_back( x );
    m_mono_Pix_y.push_back( y );
    m_mono_Pix_z.push_back( z );
    m_mono_Pix_rho.push_back( rho );
    m_mono_Pix_eta.push_back( monoPixST.eta(m,i) );
    m_mono_Pix_phi.push_back( monoPixST.phi(m,i) );

    const double tof = monoPixST.hit(m,i) ? monoPixST.hit(m,i)->timeOfFlight() : -999;
    const double energy = monoPixST.hit(m,i) ? monoPixST.hit(m,i)->energyLoss() : -999;

    m_mono_Pix_tof.push_back( tof );
    m_mono_Pix_energy.push_back( energy );

    // calculate path length in tracker element
    if ( monoPixST.hit(m,i) ) {

      LocalPoint entry = monoPixST.hit(m,i)->entryPoint();
      LocalPoint exit =  monoPixST.hit(m,i)->exitPoint();

      LocalVector path( exit - entry );

      m_mono_Pix_length.push_back( path.mag() );

    }

    
   
  }
  m_mono_Pix_N = monoPixST.size(Mono::monopole);


  // aggregate summed information
  const std::map<unsigned, Mono::SumStruct> idSumMap = monoPixST.idSumMap(m);
  std::map<unsigned, Mono::SumStruct>::const_iterator it = idSumMap.begin();
  for ( ; it != idSumMap.end(); it++ ) {
    m_mono_PixSum_id.push_back( (*it).first );
    m_mono_PixSum_N.push_back( (*it).second.N );
    m_mono_PixSum_energy.push_back( (*it).second.sum );
    m_mono_PixSum_Nids++;
    m_mono_PixSum_eta.push_back( monoPixST.eta( (*it).first ) );
    m_mono_PixSum_phi.push_back( monoPixST.phi( (*it).first ) );
  }





  // fill anti-monopole Pixel SimHit vectors
  for ( unsigned i=0; i != monoPixST.size(Mono::anti_monopole); i++ ) {

    const double x = monoPixST.x(a,i);
    const double y = monoPixST.y(a,i);
    const double z = monoPixST.z(a,i);

    const double rho = sqrt(x*x + y*y);

    m_amon_Pix_x.push_back( x );
    m_amon_Pix_y.push_back( y );
    m_amon_Pix_z.push_back( z );
    m_amon_Pix_rho.push_back( rho );
    m_amon_Pix_eta.push_back( monoPixST.eta(a,i) );
    m_amon_Pix_phi.push_back( monoPixST.phi(a,i) );

    const double tof = monoPixST.hit(a,i) ? monoPixST.hit(a,i)->timeOfFlight() : -999;
    const double energy = monoPixST.hit(a,i) ? monoPixST.hit(a,i)->energyLoss() : -999;

    m_amon_Pix_tof.push_back( tof );
    m_amon_Pix_energy.push_back( energy );


    // calculate path length in tracker element
    if ( monoPixST.hit(a,i) ) {

      LocalPoint entry = monoPixST.hit(a,i)->entryPoint();
      LocalPoint exit =  monoPixST.hit(a,i)->exitPoint();

      LocalVector path( exit - entry );

      m_amon_Pix_length.push_back( path.mag() );

    }

   




  }
  m_amon_Pix_N = monoPixST.size(Mono::anti_monopole);



  // aggregate summed information
  const std::map<unsigned,Mono::SumStruct> amonPixIdSumMap = monoPixST.idSumMap(a);
  it = amonPixIdSumMap.begin();
  for ( ; it != amonPixIdSumMap.end(); it++ ) {
    m_amon_PixSum_id.push_back( (*it).first );
    m_amon_PixSum_N.push_back( (*it).second.N );
    m_amon_PixSum_energy.push_back( (*it).second.sum );
    m_amon_PixSum_Nids++;
    m_amon_PixSum_eta.push_back( monoPixST.eta( (*it).first ) );
    m_amon_PixSum_phi.push_back( monoPixST.phi( (*it).first ) );
  }





  // fill monopole Ecal SimHit vectors
  for ( unsigned i=0; i != monoEcalST.size(Mono::monopole); i++ ) {

    const double x = monoEcalST.x(m,i);
    const double y = monoEcalST.y(m,i);
    const double z = monoEcalST.z(m,i);
    
    const double rho = sqrt(x*x + y*y);

    m_mono_Ecal_x.push_back( x );
    m_mono_Ecal_y.push_back( y );
    m_mono_Ecal_z.push_back( z );
    m_mono_Ecal_rho.push_back( rho );
    m_mono_Ecal_eta.push_back( monoEcalST.eta(m,i) );
    m_mono_Ecal_phi.push_back( monoEcalST.phi(m,i) );

    const double time = monoEcalST.hit(m,i) ? monoEcalST.hit(m,i)->time() : -999;
    const double energy = monoEcalST.hit(m,i) ? monoEcalST.hit(m,i)->energy() : -999;

    m_mono_Ecal_time.push_back( time );
    m_mono_Ecal_energy.push_back( energy );





  }
  m_mono_Ecal_N = monoEcalST.size(Mono::monopole);
  

 
  // get monopole aggregate data 
  const std::map<unsigned, Mono::SumStruct> monoEcalIdSumMap = monoEcalST.idSumMap(m);
  it = monoEcalIdSumMap.begin();  
  for ( ; it != monoEcalIdSumMap.end(); it++ ) {
    m_mono_EcalSum_id.push_back( (*it).first );
    m_mono_EcalSum_N.push_back( (*it).second.N );
    m_mono_EcalSum_energy.push_back( (*it).second.sum );
    m_mono_EcalSum_Nids++;
    m_mono_EcalSum_eta.push_back( monoEcalST.eta( (*it).first ) );
    m_mono_EcalSum_phi.push_back( monoEcalST.phi( (*it).first) );
  }



  // fill anti-monopole Ecal SimHit vectors
  for ( unsigned i=0; i != monoEcalST.size(Mono::anti_monopole); i++ ) {

    const double x = monoEcalST.x(a,i);
    const double y = monoEcalST.y(a,i);
    const double z = monoEcalST.z(a,i);
    
    const double rho = sqrt(x*x + y*y);

    m_amon_Ecal_x.push_back( x );
    m_amon_Ecal_y.push_back( y );
    m_amon_Ecal_z.push_back( z );
    m_amon_Ecal_rho.push_back( rho );
    m_amon_Ecal_eta.push_back( monoEcalST.eta(a,i) );
    m_amon_Ecal_phi.push_back( monoEcalST.phi(a,i) );

    const double time = monoEcalST.hit(a,i) ? monoEcalST.hit(a,i)->time() : -999;
    const double energy = monoEcalST.hit(a,i) ? monoEcalST.hit(a,i)->energy() : -999;

    m_amon_Ecal_time.push_back( time );
    m_amon_Ecal_energy.push_back( energy );



  }
  m_amon_Ecal_N = monoEcalST.size(Mono::anti_monopole);


  // get anti-monopole aggregate data
  const std::map<unsigned,Mono::SumStruct> amonEcalIdSumMap = monoEcalST.idSumMap(a);
  it = amonEcalIdSumMap.begin();
  for ( ; it != amonEcalIdSumMap.end(); it++ ) {
    m_amon_EcalSum_id.push_back( (*it).first );
    m_amon_EcalSum_N.push_back( (*it).second.N );
    m_amon_EcalSum_energy.push_back( (*it).second.sum );
    m_amon_EcalSum_Nids++;
    m_amon_EcalSum_eta.push_back( monoEcalST.eta( (*it).first ) );
    m_amon_EcalSum_phi.push_back( monoEcalST.phi( (*it).first) );
  }


  // find generator level information
  Mono::MonoTruthSnoop snoopy(iEvent,iSetup);
  const HepMC::GenParticle *mono = snoopy.mono(Mono::monopole);
  const HepMC::GenParticle *amon = snoopy.mono(Mono::anti_monopole);

  if ( mono ) {
    m_mono_p = cow::mag( mono->momentum().px(), mono->momentum().py(), mono->momentum().pz() );
    m_mono_eta = mono->momentum().eta();
    m_mono_phi = mono->momentum().phi();
    m_mono_m = mono->momentum().m();
  }
 
  if ( amon ) { 
    m_amon_p = cow::mag( amon->momentum().px(), amon->momentum().py(), amon->momentum().pz() );
    m_amon_eta = amon->momentum().eta();
    m_amon_phi = amon->momentum().phi();
    m_amon_m = amon->momentum().m(); 
  }

  
  

  m_tree->Fill();

}


// ----------- strandAnalysis -----------------------------
void MonoSimAnalysis::strandAnalysis()
{

}


// ------------ method called once each job just before starting event loop  ------------
void 
MonoSimAnalysis::beginJob()
{
 
  m_tree = m_fs->make<TTree>("SimTree","SimTree");


  m_tree->Branch("mono_Ecal_N", &m_mono_Ecal_N, "mono_Ecal_N/i");
  m_tree->Branch("mono_Ecal_x", &m_mono_Ecal_x);
  m_tree->Branch("mono_Ecal_y", &m_mono_Ecal_y);
  m_tree->Branch("mono_Ecal_z", &m_mono_Ecal_z);
  m_tree->Branch("mono_Ecal_rho", &m_mono_Ecal_rho);
  m_tree->Branch("mono_Ecal_eta", &m_mono_Ecal_eta);
  m_tree->Branch("mono_Ecal_phi", &m_mono_Ecal_phi);
  m_tree->Branch("mono_Ecal_time", &m_mono_Ecal_time);
  m_tree->Branch("mono_Ecal_energy", &m_mono_Ecal_energy);
  m_tree->Branch("mono_Ecal_id", &m_mono_Ecal_id);

  m_tree->Branch("amon_Ecal_N", &m_amon_Ecal_N, "amon_Ecal_N/i");
  m_tree->Branch("amon_Ecal_x", &m_amon_Ecal_x);
  m_tree->Branch("amon_Ecal_y", &m_amon_Ecal_y);
  m_tree->Branch("amon_Ecal_z", &m_amon_Ecal_z);
  m_tree->Branch("amon_Ecal_rho", &m_amon_Ecal_rho);
  m_tree->Branch("amon_Ecal_eta", &m_amon_Ecal_eta);
  m_tree->Branch("amon_Ecal_phi", &m_amon_Ecal_phi);
  m_tree->Branch("amon_Ecal_time", &m_amon_Ecal_time);
  m_tree->Branch("amon_Ecal_energy", &m_amon_Ecal_energy);

  m_tree->Branch("mono_EcalSum_Nids", &m_mono_EcalSum_Nids, "mono_EcalSum_Nids/i");
  m_tree->Branch("mono_EcalSum_id", &m_mono_EcalSum_id);
  m_tree->Branch("mono_EcalSum_N", &m_mono_EcalSum_N);
  m_tree->Branch("mono_EcalSum_energy", &m_mono_EcalSum_energy);
  m_tree->Branch("mono_EcalSum_eta", &m_mono_EcalSum_eta);
  m_tree->Branch("mono_EcalSum_phi", &m_mono_EcalSum_phi);

  m_tree->Branch("amon_EcalSum_Nids", &m_amon_EcalSum_Nids, "amon_EcalSum_Nids/i");
  m_tree->Branch("amon_EcalSum_id", &m_amon_EcalSum_id);
  m_tree->Branch("amon_EcalSum_N", &m_amon_EcalSum_N);
  m_tree->Branch("amon_EcalSum_energy", &m_amon_EcalSum_energy);
  m_tree->Branch("amon_EcalSum_eta", &m_amon_EcalSum_eta);
  m_tree->Branch("amon_EcalSum_phi", &m_amon_EcalSum_phi);


  m_tree->Branch("mono_Pix_N", &m_mono_Pix_N, "mono_Pix_N/i");
  m_tree->Branch("mono_Pix_x", &m_mono_Pix_x);
  m_tree->Branch("mono_Pix_y", &m_mono_Pix_y);
  m_tree->Branch("mono_Pix_z", &m_mono_Pix_z);
  m_tree->Branch("mono_Pix_rho", &m_mono_Pix_rho);
  m_tree->Branch("mono_Pix_eta", &m_mono_Pix_eta);
  m_tree->Branch("mono_Pix_phi", &m_mono_Pix_phi);
  m_tree->Branch("mono_Pix_tof", &m_mono_Pix_tof);
  m_tree->Branch("mono_Pix_energy", &m_mono_Pix_energy);
  m_tree->Branch("mono_Pix_length", &m_mono_Pix_length);

  m_tree->Branch("amon_Pix_N", &m_amon_Pix_N, "amon_Pix_N/i");
  m_tree->Branch("amon_Pix_x", &m_amon_Pix_x);
  m_tree->Branch("amon_Pix_y", &m_amon_Pix_y);
  m_tree->Branch("amon_Pix_z", &m_amon_Pix_z);
  m_tree->Branch("amon_Pix_rho", &m_amon_Pix_rho);
  m_tree->Branch("amon_Pix_eta", &m_amon_Pix_eta);
  m_tree->Branch("amon_Pix_phi", &m_amon_Pix_phi);
  m_tree->Branch("amon_Pix_tof", &m_amon_Pix_tof);
  m_tree->Branch("amon_Pix_energy", &m_amon_Pix_energy);
  m_tree->Branch("amon_Pix_length", &m_amon_Pix_length);


  m_tree->Branch("mono_PixSum_Nids", &m_mono_PixSum_Nids, "mono_PixSum_Nids/i");
  m_tree->Branch("mono_PixSum_id", &m_mono_PixSum_id);
  m_tree->Branch("mono_PixSum_N", &m_mono_PixSum_N);
  m_tree->Branch("mono_PixSum_energy", &m_mono_PixSum_energy);
  m_tree->Branch("mono_PixSum_eta", &m_mono_PixSum_eta);
  m_tree->Branch("mono_PixSum_phi", &m_mono_PixSum_phi);

  m_tree->Branch("amon_PixSum_Nids", &m_amon_PixSum_Nids, "amon_PixSum_Nids/i");
  m_tree->Branch("amon_PixSum_id", &m_amon_PixSum_id);
  m_tree->Branch("amon_PixSum_N", &m_amon_PixSum_N);
  m_tree->Branch("amon_PixSum_energy", &m_amon_PixSum_energy);
  m_tree->Branch("amon_PixSum_eta", &m_amon_PixSum_eta);
  m_tree->Branch("amon_PixSum_phi", &m_amon_PixSum_phi);


  m_tree->Branch("mono_p", &m_mono_p, "mono_p/D");
  m_tree->Branch("mono_eta", &m_mono_eta, "mono_eta/D");
  m_tree->Branch("mono_phi", &m_mono_phi, "mono_phi/D");
  m_tree->Branch("mono_m", &m_mono_m, "mono_m/D");

  m_tree->Branch("amon_p", &m_amon_p, "amon_p/D");
  m_tree->Branch("amon_eta", &m_amon_eta, "amon_eta/D");
  m_tree->Branch("amon_phi", &m_amon_phi, "amon_phi/D");
  m_tree->Branch("amon_m", &m_amon_m, "amon_m/D");


}

// ------------ method called once each job just after ending the event loop  ------------
void 
MonoSimAnalysis::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MonoSimAnalysis::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MonoSimAnalysis::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MonoSimAnalysis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MonoSimAnalysis::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MonoSimAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


// clear
void
MonoSimAnalysis::clear()
{

  m_mono_Ecal_N = 0;
  m_mono_Ecal_x.clear();
  m_mono_Ecal_y.clear();
  m_mono_Ecal_z.clear();
  m_mono_Ecal_rho.clear();
  m_mono_Ecal_eta.clear();
  m_mono_Ecal_phi.clear();
  m_mono_Ecal_time.clear();
  m_mono_Ecal_energy.clear();
  m_mono_Ecal_id.clear();

  m_amon_Ecal_N = 0;
  m_amon_Ecal_x.clear();
  m_amon_Ecal_y.clear();
  m_amon_Ecal_z.clear();
  m_amon_Ecal_rho.clear();
  m_amon_Ecal_eta.clear();
  m_amon_Ecal_phi.clear();
  m_amon_Ecal_time.clear();
  m_amon_Ecal_energy.clear();

  m_mono_EcalSum_Nids = 0;
  m_mono_EcalSum_id.clear();
  m_mono_EcalSum_N.clear();
  m_mono_EcalSum_energy.clear();
  m_mono_EcalSum_eta.clear();
  m_mono_EcalSum_phi.clear();


  m_amon_EcalSum_Nids = 0;
  m_amon_EcalSum_id.clear();
  m_amon_EcalSum_N.clear();
  m_amon_EcalSum_energy.clear();
  m_amon_EcalSum_eta.clear();
  m_amon_EcalSum_phi.clear();


  m_mono_Pix_N = 0;
  m_mono_Pix_x.clear();
  m_mono_Pix_y.clear();
  m_mono_Pix_z.clear();
  m_mono_Pix_rho.clear();
  m_mono_Pix_eta.clear();
  m_mono_Pix_phi.clear();
  m_mono_Pix_tof.clear();
  m_mono_Pix_energy.clear();
  m_mono_Pix_length.clear();

  m_amon_Pix_N = 0;
  m_amon_Pix_x.clear();
  m_amon_Pix_y.clear();
  m_amon_Pix_z.clear();
  m_amon_Pix_rho.clear();
  m_amon_Pix_eta.clear();
  m_amon_Pix_phi.clear();
  m_amon_Pix_tof.clear();
  m_amon_Pix_energy.clear();
  m_amon_Pix_length.clear();



  m_mono_PixSum_Nids = 0;
  m_mono_PixSum_id.clear();
  m_mono_PixSum_N.clear();
  m_mono_PixSum_energy.clear();
  m_mono_PixSum_eta.clear();
  m_mono_PixSum_phi.clear();

  m_amon_PixSum_Nids = 0;
  m_amon_PixSum_id.clear();
  m_amon_PixSum_N.clear();
  m_amon_PixSum_energy.clear();
  m_amon_PixSum_eta.clear();
  m_amon_PixSum_phi.clear();



  m_mono_p = 0.;
  m_mono_eta = 0.;
  m_mono_phi = 0.;
  m_mono_m = 0.;

  m_amon_p = 0.;
  m_amon_eta = 0.;
  m_amon_phi = 0.;
  m_amon_m = 0.;


}

//define this as a plug-in
DEFINE_FWK_MODULE(MonoSimAnalysis);
