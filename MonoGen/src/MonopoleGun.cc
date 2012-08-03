

#include "Monopoles/MonoGen/interface/MonopoleGun.h"


#include <iostream>
#include <string>


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"





using namespace Mono;

MonopoleGun::MonopoleGun(const edm::ParameterSet& pset )
 :m_pid(pset.getParameter<int>("PID"))
 ,m_mass(0.)
 ,m_constructdir(pset.getParameter<std::vector<double> >("momentum"))
 ,m_constructpos(pset.getParameter<std::vector<double> >("location"))
{ 



  // find global momentum vector
  double tmp_dir[3] = {0.};

  for ( unsigned i=0; i != m_constructdir.size(); i++ ) 
    tmp_dir[i] = m_constructdir[i];

  
  m_dir = GlobalVector(tmp_dir[0],tmp_dir[1],tmp_dir[2]);

 
  // find global position vector 
  double tmp_pos[3] = {0.};

  for ( unsigned i=0; i != m_constructpos.size(); i++ )
    tmp_pos[i] = m_constructpos[i];

  m_pos = GlobalPoint(tmp_pos[0],tmp_pos[1],tmp_pos[2]);



  produces<edm::HepMCProduct>();
}



MonopoleGun::~MonopoleGun()
{ }





void MonopoleGun::beginJob()
{ }


void MonopoleGun::endJob()
{ }



void MonopoleGun::beginRun( edm::Run &r, edm::EventSetup const &es )
{ }


void MonopoleGun::endRun( edm::Run &r, edm::EventSetup const &es )
{ }




void MonopoleGun::produce( edm::Event &evt, edm::EventSetup const &es )
{ 


  //find mass
  edm::ESHandle<ParticleDataTable> pdt;
  es.getData( pdt );

  HepPDT::ParticleDataTable::const_iterator p=pdt->begin();
  for ( ; p != pdt->end(); ++p ) {

    HepPDT::ParticleData particle = (p->second);
    std::string particleName = (particle.name()).substr(0,8);
    if ( particleName.find("Monopole") != std::string::npos )
      m_mass = particle.mass();

  }

  // create GenEvent
  m_evt = new HepMC::GenEvent();


  // generate the event
  generateEvent();

  m_evt->set_beam_particles(0,0);
  m_evt->set_event_number(evt.id().event());
  m_evt->set_signal_process_id( 0 );


  // load the event into the event record
  loadEvent(evt);

  //m_evt->print();


}




/// load HepMC::GenEvent into the event record
void MonopoleGun::loadEvent( edm::Event &evt )
{
 
 std::auto_ptr<edm::HepMCProduct> bare_product(new edm::HepMCProduct());  

 if(m_evt)  bare_product->addHepMCData( m_evt );

 evt.put(bare_product);

}



/// generate the event
void MonopoleGun::generateEvent()
{ 

  // the primary gun vertex, pass the location vector as argument to constructor
  HepMC::GenVertex *vtx = new HepMC::GenVertex( HepMC::FourVector(m_pos.x(),m_pos.y(),m_pos.z() ) );

  double px=m_dir.x(), py=m_dir.y(), pz=m_dir.z(), E=0.;


  // find E
  E = sqrt( m_mass*m_mass+m_dir.mag2() );


  HepMC::FourVector p(px,py,pz,E);
  HepMC::GenParticle *part = new HepMC::GenParticle(p,m_pid,1);
  part->suggest_barcode( 1 );

  vtx->add_particle_out( part );  

  m_evt->add_vertex( vtx );



}



// define this as a plug-in
DEFINE_FWK_MODULE(MonopoleGun);
