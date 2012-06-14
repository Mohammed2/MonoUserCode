#ifndef MONO_MONOPOLEGUN_H
#define MONO_MONOPOLEGUN_H

////////////////////////////////////////////////////////////////
// C S Cowde 					12 March 2012
// Monopole gun EDProducer class to shoot monopoles at 
// particular detector elements.
//
// Pass it a Point and a vector for position and direction 
// respectively.
// The particle momentum must be given to determine kinematics.
////////////////////////////////////////////////////////////////


#include <vector>
#include <string>

#include "HepMC/GenEvent.h"

#include "FWCore/Framework/interface/EDProducer.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"


namespace Mono {


class MonopoleGun : public edm::EDProducer {


  public:
    MonopoleGun( const edm::ParameterSet& );
    virtual ~MonopoleGun();
    void beginJob() ;
    void endJob();
    void beginRun( edm::Run &, edm::EventSetup const& );
    void endRun( edm::Run &, edm::EventSetup const& );
    void produce( edm::Event&, const edm::EventSetup& ) ;

  private:

    void loadEvent( edm::Event& );
    void generateEvent();



    HepMC::GenEvent * m_evt;


    int m_pid;
    double m_mass;

    std::vector<double> m_constructdir;
    std::vector<double> m_constructpos;
 
    // momentum vector
    GlobalVector m_dir;

    // position of gun
    GlobalPoint  m_pos;



};  // end MonopoleGun class declaration



} // end Mono namespace




#endif
