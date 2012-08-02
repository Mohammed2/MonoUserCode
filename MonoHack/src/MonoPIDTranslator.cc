// -*- C++ -*-
//
// Package:    MonoPIDTranslator
// Class:      MonoPIDTranslator
// 
/**\class MonoPIDTranslator MonoPIDTranslator.cc Monopoles/MonoHack/src/MonoPIDTranslator.cc

 Description: Translate 4th generation leptoin PID to monopole PID based on generated mass.

 Implementation:
The MonoPIDTranslator class translates the PID of the 4th generation lepton generated by
Pythia (PID 17) to the monopole PID.  The analyzer accomplishes this by taking the monopole
mass from the python configuration script and looking for the correct PID from the 
particle data table.  Each monopole mass has its own unique PID starting at 4110000.
If the translator does not find the given mass, it throws a cms exception, thus terminating
the run.
*/
//
// Original Author:  Christopher Cowden
//         Created:  Wed Aug  1 14:25:42 CST 2012
// $Id: MonoPIDTranslator.cc,v 1.1 2012/06/14 18:23:18 cowden Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"


//data formats
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

// HEP particle data table
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

//
// class declaration
//

class MonoPIDTranslator : public edm::EDAnalyzer {
   public:
      explicit MonoPIDTranslator(const edm::ParameterSet&);
      ~MonoPIDTranslator();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
  const edm::InputTag m_mc_label;

  double m_genMass; // mass of the generated particle (given by config file)

  int m_PID;  // pid of monopole for given mass (m_genMass)

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
MonoPIDTranslator::MonoPIDTranslator(const edm::ParameterSet& iConfig)
  :m_mc_label(iConfig.getParameter<edm::InputTag>("MC_Label"))
  ,m_genMass(iConfig.getParameter<double>("Mass"))
  ,m_PID(0)
{
   //now do what ever initialization is needed


}


MonoPIDTranslator::~MonoPIDTranslator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MonoPIDTranslator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


  if ( m_PID == 0 ) {
    ESHandle<ParticleDataTable> pdt;
    iSetup.getData( pdt );
  
    HepPDT::ParticleDataTable::const_iterator p=pdt->begin();
    for ( ; p != pdt->end(); ++p ) {
  
      HepPDT::ParticleData particle = (p->second);
      std::string particleName = (particle.name()).substr(0,8);
      if ( particleName.find("Monopole") != std::string::npos && particle.mass() == m_genMass )
        m_PID = abs(particle.pid());
  
    }
  
    if ( m_PID == 0 ) {
      throw cms::Exception("Undefined monopole mass") << "The monopole mass " << m_genMass << " is not defined in the particle data table!!";
    }
  }




  edm::Handle<edm::HepMCProduct> mcproduct;
  iEvent.getByLabel(m_mc_label, mcproduct);
  HepMC::GenEvent *mc = const_cast<HepMC::GenEvent *>(mcproduct->GetEvent());
  assert(mc);


  HepMC::GenEvent::particle_iterator p = mc->particles_begin();
  for ( ; p != mc->particles_end(); ++p ) {

    if ( abs( (*p)->pdg_id()) == 17 ) {
      //std::cout << "Monopole: " << (*p)->pdg_id() << " " ;
      (*p)->set_pdg_id( m_PID*(*p)->pdg_id()/abs((*p)->pdg_id()) );
      //std::cout << (*p)->pdg_id() << std::endl;
    }

  } 
  
  

   
}


// ------------ method called once each job just before starting event loop  ------------
void 
MonoPIDTranslator::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MonoPIDTranslator::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MonoPIDTranslator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MonoPIDTranslator::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MonoPIDTranslator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MonoPIDTranslator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MonoPIDTranslator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MonoPIDTranslator);