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
// $Id: MonoSimAnalysis.cc,v 1.7 2013/03/17 12:44:31 cowden Exp $
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
#include "DataFormats/Math/interface/deltaR.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"


// Monopole includes
#include "Monopoles/MonoAlgorithms/interface/MonoDefs.h"
#include "Monopoles/MonoAlgorithms/interface/MonoSimTracker.h"
#include "Monopoles/MonoAlgorithms/interface/MonoTruthSnooper.h"
#include "Monopoles/MonoAlgorithms/interface/MonoGenTrackExtrapolator.h"


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
    edm::InputTag m_TagEEHits;


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
    
    unsigned m_mono_EE_N;
    std::vector<double> m_mono_EE_x;
    std::vector<double> m_mono_EE_y;
    std::vector<double> m_mono_EE_z;
    std::vector<double> m_mono_EE_rho;
    std::vector<double> m_mono_EE_eta;
    std::vector<double> m_mono_EE_phi;
    std::vector<double> m_mono_EE_time;
    std::vector<double> m_mono_EE_energy;
    std::vector<unsigned> m_mono_EE_id;

    unsigned m_amon_EE_N;
    std::vector<double> m_amon_EE_x;
    std::vector<double> m_amon_EE_y;
    std::vector<double> m_amon_EE_z;
    std::vector<double> m_amon_EE_rho;
    std::vector<double> m_amon_EE_eta;
    std::vector<double> m_amon_EE_phi;
    std::vector<double> m_amon_EE_time;
    std::vector<double> m_amon_EE_energy;


    unsigned m_mono_EESum_Nids;
    std::vector<unsigned> m_mono_EESum_id;
    std::vector<unsigned> m_mono_EESum_N;
    std::vector<double>   m_mono_EESum_energy;
    std::vector<double>   m_mono_EESum_eta;
    std::vector<double>   m_mono_EESum_phi;

    unsigned m_amon_EESum_Nids;
    std::vector<unsigned> m_amon_EESum_id;
    std::vector<unsigned> m_amon_EESum_N;
    std::vector<double>   m_amon_EESum_energy;
    std::vector<double>   m_amon_EESum_eta;
    std::vector<double>   m_amon_EESum_phi;


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
    double m_mono_px;
    double m_mono_py;
    double m_mono_pz;
    double m_mono_x;
    double m_mono_y;
    double m_mono_z;

    double m_monoExp_eta;
    double m_monoExp_phi;

    double m_amon_p;
    double m_amon_eta;
    double m_amon_phi;
    double m_amon_m;
    double m_amon_px;
    double m_amon_py;
    double m_amon_pz;
    double m_amon_x;
    double m_amon_y;
    double m_amon_z;

    double m_amonExp_eta;
    double m_amonExp_phi;

    double m_monoVirt_m;
    double m_monoVirt_eta;
    double m_monoVirt_phi;
    double m_monoVirt_p;
    double m_monoVirt_pt;
    double m_monoVirt_px; 
    double m_monoVirt_py;
    double m_monoVirt_pz;
    double m_monoVirt_E;
    double m_monoVirt_dRmono;
    
    double m_amonVirt_m;
    double m_amonVirt_eta;
    double m_amonVirt_phi;
    double m_amonVirt_p;
    double m_amonVirt_pt;
    double m_amonVirt_px; 
    double m_amonVirt_py;
    double m_amonVirt_pz;
    double m_amonVirt_E;
    double m_amonVirt_dRamon;

    unsigned m_md_N;
    unsigned m_ad_N;
    std::vector<double> m_md_eta;
    std::vector<double> m_md_phi;
    std::vector<double> m_md_pdgId;
    std::vector<double> m_md_pt;
    std::vector<double> m_ad_eta;
    std::vector<double> m_ad_phi;
    std::vector<double> m_ad_pdgId;
    std::vector<double> m_ad_pt;


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
  Mono::MonoSimTracker<PCaloHit,CaloGeometry> monoEEST(iEvent,iSetup,Mono::EcalEE);

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


  {
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
  std::map<unsigned, Mono::SumStruct>::const_iterator it = monoEcalIdSumMap.begin();  
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


  // fill monopole EE SimHit vectors
  for ( unsigned i=0; i != monoEEST.size(Mono::monopole); i++ ) {

    const double x = monoEEST.x(m,i);
    const double y = monoEEST.y(m,i);
    const double z = monoEEST.z(m,i);
    
    const double rho = sqrt(x*x + y*y);

    m_mono_EE_x.push_back( x );
    m_mono_EE_y.push_back( y );
    m_mono_EE_z.push_back( z );
    m_mono_EE_rho.push_back( rho );
    m_mono_EE_eta.push_back( monoEEST.eta(m,i) );
    m_mono_EE_phi.push_back( monoEEST.phi(m,i) );

    const double time = monoEEST.hit(m,i) ? monoEEST.hit(m,i)->time() : -999;
    const double energy = monoEEST.hit(m,i) ? monoEEST.hit(m,i)->energy() : -999;

    m_mono_EE_time.push_back( time );
    m_mono_EE_energy.push_back( energy );





  }
  m_mono_EE_N = monoEEST.size(Mono::monopole);
  

 
  // get monopole aggregate data 
  const std::map<unsigned, Mono::SumStruct> monoEEIdSumMap = monoEEST.idSumMap(m);
  it = monoEEIdSumMap.begin();  
  for ( ; it != monoEEIdSumMap.end(); it++ ) {
    m_mono_EESum_id.push_back( (*it).first );
    m_mono_EESum_N.push_back( (*it).second.N );
    m_mono_EESum_energy.push_back( (*it).second.sum );
    m_mono_EESum_Nids++;
    m_mono_EESum_eta.push_back( monoEEST.eta( (*it).first ) );
    m_mono_EESum_phi.push_back( monoEEST.phi( (*it).first) );
  }



  // fill anti-monopole EE SimHit vectors
  for ( unsigned i=0; i != monoEEST.size(Mono::anti_monopole); i++ ) {

    const double x = monoEEST.x(a,i);
    const double y = monoEEST.y(a,i);
    const double z = monoEEST.z(a,i);
    
    const double rho = sqrt(x*x + y*y);

    m_amon_EE_x.push_back( x );
    m_amon_EE_y.push_back( y );
    m_amon_EE_z.push_back( z );
    m_amon_EE_rho.push_back( rho );
    m_amon_EE_eta.push_back( monoEEST.eta(a,i) );
    m_amon_EE_phi.push_back( monoEEST.phi(a,i) );

    const double time = monoEEST.hit(a,i) ? monoEEST.hit(a,i)->time() : -999;
    const double energy = monoEEST.hit(a,i) ? monoEEST.hit(a,i)->energy() : -999;

    m_amon_EE_time.push_back( time );
    m_amon_EE_energy.push_back( energy );



  }
  m_amon_EE_N = monoEEST.size(Mono::anti_monopole);


  // get anti-monopole aggregate data
  const std::map<unsigned,Mono::SumStruct> amonEEIdSumMap = monoEEST.idSumMap(a);
  it = amonEEIdSumMap.begin();
  for ( ; it != amonEEIdSumMap.end(); it++ ) {
    m_amon_EESum_id.push_back( (*it).first );
    m_amon_EESum_N.push_back( (*it).second.N );
    m_amon_EESum_energy.push_back( (*it).second.sum );
    m_amon_EESum_Nids++;
    m_amon_EESum_eta.push_back( monoEEST.eta( (*it).first ) );
    m_amon_EESum_phi.push_back( monoEEST.phi( (*it).first) );
  }

    
  // find generator level information
  Mono::MonoTruthSnoop snoopy(iEvent,iSetup);
  const HepMC::GenParticle *mono = snoopy.mono(Mono::monopole);
  const HepMC::GenParticle *amon = snoopy.mono(Mono::anti_monopole);

  Mono::MonoGenTrackExtrapolator extrap;

  if ( mono ) {
    m_mono_p = cow::mag( mono->momentum().px(), mono->momentum().py(), mono->momentum().pz() );
    m_mono_eta = mono->momentum().eta();
    m_mono_phi = mono->momentum().phi();
    m_mono_m = mono->momentum().m();
    m_mono_px = mono->momentum().px();
    m_mono_py = mono->momentum().py();
    m_mono_pz = mono->momentum().pz();
    m_mono_x = mono->momentum().x();
    m_mono_y = mono->momentum().y();
    m_mono_z = mono->momentum().z();

    extrap.setMonopole(*mono);
    m_monoExp_eta = extrap.etaVr(1.29);
    m_monoExp_phi = extrap.phi();
  }
 
  if ( amon ) { 
    m_amon_p = cow::mag( amon->momentum().px(), amon->momentum().py(), amon->momentum().pz() );
    m_amon_eta = amon->momentum().eta();
    m_amon_phi = amon->momentum().phi();
    m_amon_m = amon->momentum().m(); 
    m_amon_px = amon->momentum().px();
    m_amon_py = amon->momentum().py();
    m_amon_pz = amon->momentum().pz();
    m_amon_x = amon->momentum().x();
    m_amon_y = amon->momentum().y();
    m_amon_z = amon->momentum().z();

    extrap.setMonopole(*amon);
    m_amonExp_eta = extrap.etaVr(1.29);
    m_amonExp_phi = extrap.phi();
  }

  std::vector<const HepMC::GenParticle *> d1 = snoopy.monoDaughter(Mono::monopole);
  std::vector<const HepMC::GenParticle *> d2 = snoopy.monoDaughter(Mono::anti_monopole);

  for ( unsigned i=0; i < d1.size(); i++ ) {
    std::cout << "Found monopole daughter!" << std::endl;
    m_md_N++;
    m_md_eta.push_back( d1[i]->momentum().eta() );
    m_md_phi.push_back( d1[i]->momentum().phi() );
    m_md_pdgId.push_back( d1[i]->pdg_id() );
    m_md_pt.push_back( d1[i]->momentum().perp() );
  }

  for ( unsigned i=0; i < d2.size(); i++ ) {
    std::cout << "Found monopole daughter!" << std::endl;
    m_ad_N++;
    m_ad_eta.push_back( d2[i]->momentum().eta() );
    m_ad_phi.push_back( d2[i]->momentum().phi() );
    m_ad_pdgId.push_back( d2[i]->pdg_id() );
    m_ad_pt.push_back( d2[i]->momentum().perp() );
  }
 
  const HepMC::GenParticle * monoVirt = snoopy.monoGen(Mono::monopole);
  const HepMC::GenParticle * amonVirt = snoopy.monoGen(Mono::anti_monopole);

  if ( monoVirt ) {
    const HepMC::FourVector & mom = monoVirt->momentum();
    m_monoVirt_m = mom.m();
    m_monoVirt_eta = mom.eta();
    m_monoVirt_phi = mom.phi();
    m_monoVirt_p = cow::mag(mom.px(),mom.py(),mom.pz());
    m_monoVirt_pt = mom.perp();
    m_monoVirt_px = mom.px();
    m_monoVirt_py = mom.py();
    m_monoVirt_pz = mom.pz();
    m_monoVirt_E  = mom.e();
    if ( mono ) m_monoVirt_dRmono = reco::deltaR(m_monoVirt_eta,m_monoVirt_phi,m_mono_eta,m_mono_phi);
  }

  if ( amonVirt ) {
    const HepMC::FourVector & mom = amonVirt->momentum();
    m_amonVirt_m = mom.m();
    m_amonVirt_eta = mom.eta();
    m_amonVirt_phi = mom.phi();
    m_amonVirt_p = cow::mag(mom.px(),mom.py(),mom.pz());
    m_amonVirt_pt = mom.perp();
    m_amonVirt_px = mom.px();
    m_amonVirt_py = mom.py();
    m_amonVirt_pz = mom.pz();
    m_amonVirt_E  = mom.e();
    if ( amon ) m_amonVirt_dRamon = reco::deltaR(m_amonVirt_eta,m_amonVirt_phi,m_amon_eta,m_amon_phi);
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

  m_tree->Branch("mono_EE_N", &m_mono_EE_N, "mono_EE_N/i");
  m_tree->Branch("mono_EE_x", &m_mono_EE_x);
  m_tree->Branch("mono_EE_y", &m_mono_EE_y);
  m_tree->Branch("mono_EE_z", &m_mono_EE_z);
  m_tree->Branch("mono_EE_rho", &m_mono_EE_rho);
  m_tree->Branch("mono_EE_eta", &m_mono_EE_eta);
  m_tree->Branch("mono_EE_phi", &m_mono_EE_phi);
  m_tree->Branch("mono_EE_time", &m_mono_EE_time);
  m_tree->Branch("mono_EE_energy", &m_mono_EE_energy);
  m_tree->Branch("mono_EE_id", &m_mono_EE_id);

  m_tree->Branch("amon_EE_N", &m_amon_EE_N, "amon_EE_N/i");
  m_tree->Branch("amon_EE_x", &m_amon_EE_x);
  m_tree->Branch("amon_EE_y", &m_amon_EE_y);
  m_tree->Branch("amon_EE_z", &m_amon_EE_z);
  m_tree->Branch("amon_EE_rho", &m_amon_EE_rho);
  m_tree->Branch("amon_EE_eta", &m_amon_EE_eta);
  m_tree->Branch("amon_EE_phi", &m_amon_EE_phi);
  m_tree->Branch("amon_EE_time", &m_amon_EE_time);
  m_tree->Branch("amon_EE_energy", &m_amon_EE_energy);

  m_tree->Branch("mono_EESum_Nids", &m_mono_EESum_Nids, "mono_EESum_Nids/i");
  m_tree->Branch("mono_EESum_id", &m_mono_EESum_id);
  m_tree->Branch("mono_EESum_N", &m_mono_EESum_N);
  m_tree->Branch("mono_EESum_energy", &m_mono_EESum_energy);
  m_tree->Branch("mono_EESum_eta", &m_mono_EESum_eta);
  m_tree->Branch("mono_EESum_phi", &m_mono_EESum_phi);

  m_tree->Branch("amon_EESum_Nids", &m_amon_EESum_Nids, "amon_EESum_Nids/i");
  m_tree->Branch("amon_EESum_id", &m_amon_EESum_id);
  m_tree->Branch("amon_EESum_N", &m_amon_EESum_N);
  m_tree->Branch("amon_EESum_energy", &m_amon_EESum_energy);
  m_tree->Branch("amon_EESum_eta", &m_amon_EESum_eta);
  m_tree->Branch("amon_EESum_phi", &m_amon_EESum_phi);


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
  m_tree->Branch("mono_px",&m_mono_px, "mono_px/D");
  m_tree->Branch("mono_py",&m_mono_py, "mono_py/D");
  m_tree->Branch("mono_pz",&m_mono_pz, "mono_pz/D");
  m_tree->Branch("mono_x",&m_mono_x, "mono_x/D");
  m_tree->Branch("mono_y",&m_mono_y, "mono_y/D");
  m_tree->Branch("mono_z",&m_mono_z, "mono_z/D");
  
  m_tree->Branch("monoExp_eta",&m_monoExp_eta, "monoExp_eta/D");
  m_tree->Branch("monoExp_phi",&m_monoExp_phi, "monoExp_phi/D");

  m_tree->Branch("amon_p", &m_amon_p, "amon_p/D");
  m_tree->Branch("amon_eta", &m_amon_eta, "amon_eta/D");
  m_tree->Branch("amon_phi", &m_amon_phi, "amon_phi/D");
  m_tree->Branch("amon_m", &m_amon_m, "amon_m/D");
  m_tree->Branch("amon_px",&m_amon_px, "amon_px/D");
  m_tree->Branch("amon_py",&m_amon_py, "amon_py/D");
  m_tree->Branch("amon_pz",&m_amon_pz, "amon_pz/D");
  m_tree->Branch("amon_x",&m_amon_x, "amon_x/D");
  m_tree->Branch("amon_y",&m_amon_y, "amon_y/D");
  m_tree->Branch("amon_z",&m_amon_z, "amon_z/D");

  m_tree->Branch("amonExp_eta",&m_amonExp_eta, "amonExp_eta/D");
  m_tree->Branch("amonExp_phi",&m_amonExp_phi, "amonExp_phi/D");

  m_tree->Branch("monoVirt_m",&m_monoVirt_m,"monoVirt_m");
  m_tree->Branch("monoVirt_eta",&m_monoVirt_eta,"monoVirt_eta/D");
  m_tree->Branch("monoVirt_phi",&m_monoVirt_phi,"monoVirt_phi/D");
  m_tree->Branch("monoVirt_p",&m_monoVirt_p,"monoVirt_p/D");
  m_tree->Branch("monoVirt_pt",&m_monoVirt_pt,"monoVirt_pt/D");
  m_tree->Branch("monoVirt_px",&m_monoVirt_px,"monoVirt_px/D"); 
  m_tree->Branch("monoVirt_py",&m_monoVirt_py,"monoVirt_py/D");
  m_tree->Branch("monoVirt_pz",&m_monoVirt_pz,"monoVirt_pz/D");
  m_tree->Branch("monoVirt_E",&m_monoVirt_E,"monoVirt_E/D");
  m_tree->Branch("monoVirt_dRmono",&m_monoVirt_dRmono,"monoVirt_dRmono/D");

  m_tree->Branch("amonVirt_m",&m_amonVirt_m,"amonVirt_m");
  m_tree->Branch("amonVirt_eta",&m_amonVirt_eta,"amonVirt_eta/D");
  m_tree->Branch("amonVirt_phi",&m_amonVirt_phi,"amonVirt_phi/D");
  m_tree->Branch("amonVirt_p",&m_amonVirt_p,"amonVirt_p/D");
  m_tree->Branch("amonVirt_pt",&m_amonVirt_pt,"amonVirt_pt/D");
  m_tree->Branch("amonVirt_px",&m_amonVirt_px,"amonVirt_px/D"); 
  m_tree->Branch("amonVirt_py",&m_amonVirt_py,"amonVirt_py/D");
  m_tree->Branch("amonVirt_pz",&m_amonVirt_pz,"amonVirt_pz/D");
  m_tree->Branch("amonVirt_E",&m_amonVirt_E,"amonVirt_E/D");
  m_tree->Branch("amonVirt_dRamon",&m_amonVirt_dRamon,"amonVirt_dRamon/D");

  m_tree->Branch("md_N",&m_md_N,"md_N/i");
  m_tree->Branch("md_eta",&m_md_eta);
  m_tree->Branch("md_phi",&m_md_phi);
  m_tree->Branch("md_pdgId",&m_md_pdgId);
  m_tree->Branch("md_pt",&m_md_pt);
 
  m_tree->Branch("ad_N",&m_ad_N,"ad_N/i"); 
  m_tree->Branch("ad_eta",&m_ad_eta);
  m_tree->Branch("ad_phi",&m_ad_phi);
  m_tree->Branch("ad_pdgId",&m_ad_pdgId);
  m_tree->Branch("ad_pt",&m_ad_pt);

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
  m_mono_px = 0.;
  m_mono_py = 0.;
  m_mono_pz = 0.;
  m_mono_x = 0.;
  m_mono_y = 0.;
  m_mono_z = 0.;

  m_monoExp_eta = 0.;
  m_monoExp_phi = 0.;

  m_amon_p = 0.;
  m_amon_eta = 0.;
  m_amon_phi = 0.;
  m_amon_m = 0.;
  m_amon_px = 0.;
  m_amon_py = 0.;
  m_amon_pz = 0.;
  m_amon_x = 0.;
  m_amon_y = 0.;
  m_amon_z = 0.;

  m_amonExp_eta = 0.;
  m_amonExp_phi = 0.;

  m_monoVirt_m = 0.;
  m_monoVirt_eta = 0.;
  m_monoVirt_phi = 0.;
  m_monoVirt_pt = 0.;
  m_monoVirt_px = 0.;
  m_monoVirt_py = 0.;
  m_monoVirt_pz = 0.;
  m_monoVirt_E = 0.;
  m_monoVirt_dRmono = 0.;

  m_amonVirt_m = 0.;
  m_amonVirt_eta = 0.;
  m_amonVirt_phi = 0.;
  m_amonVirt_pt = 0.;
  m_amonVirt_px = 0.;
  m_amonVirt_py = 0.;
  m_amonVirt_pz = 0.;
  m_amonVirt_E = 0.;
  m_amonVirt_dRamon = 0.;

  m_md_N = 0U;
  m_md_eta.clear();
  m_md_phi.clear();
  m_md_pdgId.clear();
  m_md_pt.clear();

  m_ad_N = 0U;
  m_ad_eta.clear();
  m_ad_phi.clear();
  m_ad_pdgId.clear();
  m_ad_pt.clear();

}

//define this as a plug-in
DEFINE_FWK_MODULE(MonoSimAnalysis);
