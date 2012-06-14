// -*- C++ -*-
//
// Package:    MonoGen
// Class:      MonoGen
// 
/**\class MonoGen MonoGen.cc Monopoles/MonoGen/src/MonoGen.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christopher Cowden
//         Created:  Mon Jan 30 13:13:11 CST 2012
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>
#include <cstdio>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "pythia6.h"

// ROOT includes
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"


// reflex stuff
#ifdef __GCCXML__
template class std::vector<TLorentzVector>;
#endif

//
// class declaration
//

class MonoGen : public edm::EDAnalyzer {
   public:
      explicit MonoGen(const edm::ParameterSet&);
      ~MonoGen();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

    void clear();  // method to clear tree branch variables 

      // ----------member data ---------------------------

    /// TFileService
    edm::Service<TFileService> m_fs;


    TTree * m_kin;



    // kinematic tree branch variables
    // put these last in the class so constructor will not create member
    // data out of order.
    int    m_m1_pdgId;
    TLorentzVector m_m1;


    int    m_m2_pdgId;
    TLorentzVector m_m2;


    int    m_g1_pdgId;
    TLorentzVector m_g1;
    
    int    m_g2_pdgId;
    TLorentzVector m_g2;

    std::vector<double> m_d1_pT;
    std::vector<double> m_d1_pX;
    std::vector<double> m_d1_pY;
    std::vector<double> m_d1_pZ;
    std::vector<double> m_d1_E;
    std::vector<double> m_d1_eta;
    std::vector<double> m_d1_phi;
    std::vector<int> m_d1_pdgId;
    std::vector<TLorentzVector> m_d1;

    std::vector<double> m_d2_pT;
    std::vector<double> m_d2_pX;
    std::vector<double> m_d2_pY;
    std::vector<double> m_d2_pZ;
    std::vector<double> m_d2_E;
    std::vector<double> m_d2_eta;
    std::vector<double> m_d2_phi;
    std::vector<int> m_d2_pdgId;
    std::vector<TLorentzVector> m_d2;


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
MonoGen::MonoGen(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


MonoGen::~MonoGen()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MonoGen::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  clear();

    // Get the MC event
    edm::Handle<edm::HepMCProduct> mcproduct;
    iEvent.getByLabel(edm::InputTag("generator"), mcproduct);
    const HepMC::GenEvent* mc = mcproduct->GetEvent();
    assert(mc);

    std::cout << "\n**** HepMC Printout ****" << std::endl;



    // Cycle over MC particles
    const HepMC::GenEvent::particle_const_iterator end = mc->particles_end();
    for (HepMC::GenEvent::particle_const_iterator p = mc->particles_begin();
         p != end; ++p)
    {
        const HepMC::GenParticle* particle = *p;
        const reco::Candidate::LorentzVector p4( particle->momentum());

        if ( abs(particle->pdg_id()) != 4110000 ) continue;

     /*   printf("%4u %4d  %4u  %10.5f   %10.5f   %10.5f   %10.5f   %10.5f %10.5f\n",
               counter, particle->pdg_id(),  particle->status(), p4.Px(), p4.Py(), p4.Pz(), p4.E(),
               p4.M(),p4.Eta());*/


      // if monopole is part of the hard process
      if ( particle->status() == 3 ) {
        std::cout << "\\ " << particle->pdg_id() << " " << particle->status();
	std::cout << " : " << p4.Pt() << " " << p4.Eta() << std::endl;


	// loop over decay products
	const HepMC::GenVertex * vert = particle->end_vertex();
   	HepMC::GenVertex::particles_out_const_iterator daughter = vert->particles_out_const_begin();
        for ( ; daughter != vert->particles_out_const_end(); daughter++ ) {

	  const reco::Candidate::LorentzVector d4( (*daughter)->momentum() );

	  std::cout << "\t-> " << (*daughter)->pdg_id() << " " << (*daughter)->status();
	  std::cout << " : " << d4.Pt() << " " << d4.Eta() << std::endl;

	  if ( abs((*daughter)->pdg_id()) != 4110000 && particle->pdg_id() > 0 )  {
	    m_d1_pT.push_back( d4.Pt() );
	    m_d1_pX.push_back( d4.Px() );
	    m_d1_pY.push_back( d4.Py() );
	    m_d1_pZ.push_back( d4.Pz() );
	    m_d1_E.push_back( d4.E() );
	    m_d1_eta.push_back( d4.Eta() );
	    m_d1_phi.push_back( d4.Phi() );
	    m_d1_pdgId.push_back( (*daughter)->pdg_id() );
	    m_d1.push_back( TLorentzVector( d4.Px(), d4.Py(), d4.Pz(), d4.E() ) );
	  }
	  else if ( abs((*daughter)->pdg_id()) != 4110000 && particle->pdg_id() < 0 ) {
	    m_d2_pT.push_back( d4.Pt() );
	    m_d2_pX.push_back( d4.Px() );
	    m_d2_pY.push_back( d4.Py() );
	    m_d2_pZ.push_back( d4.Pz() );
	    m_d2_E.push_back( d4.E() );
	    m_d2_eta.push_back( d4.Eta() );
	    m_d2_phi.push_back( d4.Phi() );
	    m_d2_pdgId.push_back( (*daughter)->pdg_id() );
	    m_d2.push_back( TLorentzVector( d4.Px(), d4.Py(), d4.Pz(), d4.E() ) );
	  }

	}

	if ( particle->pdg_id() > 0 ) {
	  m_g1_pdgId = particle->pdg_id();
	  m_g1 = TLorentzVector(p4.Px(),p4.Py(),p4.Pz(),p4.E());

	} else if ( particle->pdg_id() < 0 ) {
	  m_g2_pdgId = particle->pdg_id();
	  m_g2 = TLorentzVector(p4.Px(),p4.Py(),p4.Pz(),p4.E());

	}


      } else if ( particle->status() == 1 ) {


      // write out the two stable monopole's kinematic information
      if ( particle->pdg_id() > 0 ) {
	m_m1_pdgId = particle->pdg_id();
	m_m1 = TLorentzVector(p4.Px(),p4.Py(),p4.Pz(),p4.E());
      } else if ( particle->pdg_id() < 0 ) {
	m_m2_pdgId = particle->pdg_id();
	m_m2 = TLorentzVector(p4.Px(),p4.Py(),p4.Pz(),p4.E());
      }

    } 


    }
    std::cout << std::endl;


  /*std::cout.flush();
  std::cerr.flush();
  fflush(stdout);
  fflush(stderr);

  pylist(2);

  fflush(stdout);
  fflush(stderr);*/


  m_kin->Fill();

}



void MonoGen::clear()
{
  m_m1_pdgId = 0;
  m_m1 = TLorentzVector(0.,0.,0.,0.); 

  m_m2_pdgId = 0;
  m_m2 = TLorentzVector(0.,0.,0.,0.); 


  m_g1_pdgId = 0;
  m_g1 = TLorentzVector(0.,0.,0.,0.); 

  m_g2_pdgId = 0;
  m_g2 = TLorentzVector(0.,0.,0.,0.); 

  m_d1_pT.clear();
  m_d1_pX.clear();
  m_d1_pY.clear();
  m_d1_pZ.clear();
  m_d1_E.clear();
  m_d1_eta.clear();
  m_d1_phi.clear();
  m_d1_pdgId.clear();
  m_d1.clear();

  m_d2_pT.clear();
  m_d2_pX.clear();
  m_d2_pY.clear();
  m_d2_pZ.clear();
  m_d2_E.clear();
  m_d2_eta.clear();
  m_d2_phi.clear();
  m_d2_pdgId.clear();
  m_d2.clear();

}


// ------------ method called once each job just before starting event loop  ------------
void 
MonoGen::beginJob()
{

  m_kin = m_fs->make<TTree>("kinematics","kinematics");


  // set kinematic tree branches
  m_kin->Branch("m1_pdgId",&m_m1_pdgId,"m1_pdgId/I");
  m_kin->Branch("m1",&m_m1);
  
  m_kin->Branch("m2_pdgId",&m_m2_pdgId,"m2_pdgId/I");
  m_kin->Branch("m2",&m_m2);


  m_kin->Branch("g1_pdgId",&m_g1_pdgId,"g1_pdgId/I");
  m_kin->Branch("g1",&m_g1);

  m_kin->Branch("g2_pdgId",&m_g2_pdgId,"g2_pdgId/I");
  m_kin->Branch("g2",&m_g2);


  m_kin->Branch("d1_pT",&m_d1_pT);
  m_kin->Branch("d1_pX",&m_d1_pX);
  m_kin->Branch("d1_pY",&m_d1_pY);
  m_kin->Branch("d1_pZ",&m_d1_pZ);
  m_kin->Branch("d1_E",&m_d1_E);
  m_kin->Branch("d1_eta",&m_d1_eta);
  m_kin->Branch("d1_phi",&m_d1_phi);
  m_kin->Branch("d1_pdgId",&m_d1_pdgId);
  //m_kin->Branch("d1",&m_d1);

  m_kin->Branch("d2_pT",&m_d2_pT);
  m_kin->Branch("d2_pX",&m_d2_pX);
  m_kin->Branch("d2_pY",&m_d2_pY);
  m_kin->Branch("d2_pZ",&m_d2_pZ);
  m_kin->Branch("d2_E",&m_d2_E);
  m_kin->Branch("d2_eta",&m_d2_eta);
  m_kin->Branch("d2_phi",&m_d2_phi);
  m_kin->Branch("d2_pdgId",&m_d2_pdgId);
  //m_kin->Branch("d2",&m_d2);


}


// ------------ method called once each job just after ending the event loop  ------------
void 
MonoGen::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MonoGen::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MonoGen::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MonoGen::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MonoGen::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MonoGen::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MonoGen);
