#ifndef MONOALGORITHMS_MONOECALOBS0_H
#define MONOALGORITHMS_MONOECALOBS0_H
/////////////////////////////////////////////////////////
// C S Cowden                       20 June 2012
// First monopole physics observable for CMS ECAL
/////////////////////////////////////////////////////////

#include <vector>
#include <cassert>
#include <cmath>
#include <cfloat>

#include "Monopoles/MonoAlgorithms/interface/MonoEcalSeed.h"
#include "Monopoles/MonoAlgorithms/interface/MonoEcalCluster.h"
#include "Monopoles/MonoAlgorithms/interface/ClustCategorizer.h"
#include "Monopoles/MonoAlgorithms/interface/EnergyFlowFunctor.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

// forward declarations
namespace edm {
  class Event;
  class EventSetup;
}

namespace Mono {


// ---------------------------------------------------------------------
// Ecal barrel map class
class EBmap {

public:
  inline EBmap() { constructMap(); }

  virtual inline ~EBmap() { }

  // construct the map
  inline void constructMap() 
    {
      const EBDetId did;
      m_nEta = 2*did.MAX_IETA;
      m_nPhi = did.MAX_IPHI;
      m_nCells = m_nPhi*m_nEta;

      m_minEta = -did.MAX_IETA*did.crystalUnitToEta;
      m_maxEta = did.MAX_IETA*did.crystalUnitToEta;
      m_minPhi = -M_PI;
      m_maxPhi = M_PI;

      m_etaWidth = did.crystalUnitToEta;
      m_phiWidth = 2.*M_PI/m_nPhi;

      m_ecalMap.resize(m_nCells);
      m_ecalTMap.resize(m_nCells);
    }

  // fill map
  void fillMap(const edm::Event &ev);

  // clear map
  void clear();

  // construct Ecal geometry stuff
  void constructGeo(const edm::EventSetup &);


  // accessor methods
  const inline unsigned nCells() const { return m_nCells; }
  const inline unsigned nEta() const { return m_nEta; }
  const inline unsigned nPhi() const { return m_nPhi; }

  // return the energy of the given bin
  const double inline operator[](const unsigned bin) const
   { 
     assert(m_ecalMap.size());
     assert(bin < m_nCells );
     return m_ecalMap[bin];
   }

  // return the time in the given bin
  const double inline time(const unsigned bin ) const
    {
      assert(m_ecalTMap.size());
      assert(bin < m_nCells);
      return m_ecalTMap[bin];
    }


  // return the ecalMap bin of an EBdetId
  const unsigned findBin(const EBDetId & ) const;
  const unsigned findBinEtaPhi(double,double ) const;


  // return eta and phi given eta/phi bin (not combined bin)
  const inline  double eta(unsigned i) const
    {
      assert(i < m_nEta );
      return i*m_etaWidth+m_minEta;
    }
  const inline  double phi(unsigned i) const
    {
      assert( i < m_nPhi );
      return i*m_phiWidth+m_minPhi;
    }


private:
  
  // ecalMap
  unsigned m_nCells;
  unsigned m_nEta;
  unsigned m_nPhi;
  double m_minEta;
  double m_maxEta;
  double m_minPhi;
  double m_maxPhi;
  double m_etaWidth;
  double m_phiWidth;

  // map of energy
  std::vector<double> m_ecalMap;
  // map of time
  std::vector<double> m_ecalTMap;


  // calorimetry geometry
  const CaloSubdetectorGeometry *m_geom;

};


// ----------------------------------------------------------------------
// Seed finder class based on strips in eta
class StripSeedFinder {

public:
  inline StripSeedFinder()
    :m_seedLength(3U)
    ,m_threshold(20.)
    ,m_nSeeds(0U) 
    { }

  inline StripSeedFinder(const unsigned seedLength,const double threshold, const unsigned cells)
    :m_seedLength(seedLength)
    ,m_threshold(threshold)
    ,m_nSeeds(0U)
    ,m_maxSeeds(cells)
    { }


  inline virtual ~StripSeedFinder()
   { }

  // initialize the memory for seed array
  inline void initialize() { m_seeds.resize(m_maxSeeds); m_stage.resize(m_maxSeeds); }

  // construct the seed array and obtain geometry related stuff
  void constructGeo(const edm::EventSetup &);

  // run the finder
  // return false of an error occurred
  // This function expects the following arguments:
  // edm::Event, EBmap reference
  bool find(const edm::Event &, const EBmap &);

  // clear seeds
  inline void clear() { m_nSeeds = 0U; }


  // accessor methods
  inline const unsigned nSeeds() const { return m_nSeeds; }
  inline const MonoEcalSeed * seeds() const { return &m_seeds[0]; }
  inline const unsigned seedLength() const { return m_seedLength; }

  inline const bool adjacentInPhi(const unsigned i,const unsigned j) const
    {
      assert ( i < m_nSeeds );
      assert ( j < m_nSeeds );

      const MonoEcalSeed & iSeed = m_seeds[i];
      const MonoEcalSeed & jSeed = m_seeds[j];

      const unsigned iEta = iSeed.ieta();
      const unsigned jEta = jSeed.ieta();
      const unsigned iLength = iSeed.seedLength();
      const unsigned jLength = jSeed.seedLength();

      if ( iEta - jEta > jLength ) return false;
      else if ( jEta - iEta > iLength ) return false;
  
      const unsigned iPhi = m_seeds[i].iphi(); 
      const unsigned jPhi = m_seeds[j].iphi();
      return abs(iPhi-jPhi) == 1;
    }

  const bool adjacentInEta(const unsigned i,const unsigned j) const
    {
      assert ( i < m_nSeeds );
      assert ( j < m_nSeeds );

      const MonoEcalSeed & iSeed = m_seeds[i];
      const MonoEcalSeed & jSeed = m_seeds[j];

      if ( iSeed.iphi() != jSeed.iphi() ) return false;

      const unsigned iEta = iSeed.ieta();
      const unsigned jEta = jSeed.ieta();
      const unsigned iLength = iSeed.seedLength();
      const unsigned jLength = jSeed.seedLength();

      bool res = iEta+iLength+1U == jEta;
      res = res && (jEta+jLength+1U == iEta);
      return res;
    }

  const bool overlapInEta(const unsigned i,const unsigned j) const
    {
      assert ( i < m_nSeeds );
      assert ( j < m_nSeeds );

      const MonoEcalSeed & iSeed = m_seeds[i];
      const MonoEcalSeed & jSeed = m_seeds[j];

      if ( iSeed.iphi() != jSeed.iphi() ) return false;

      const unsigned iEta = iSeed.ieta();
      const unsigned jEta = jSeed.ieta();
      const unsigned iLength = iSeed.seedLength();
      const unsigned jLength = jSeed.seedLength();

      bool res = iEta+iLength >= jEta;
      res = res && ( jEta+jLength >= iEta );
      return res;
    }

private:

  // add seed to seed list
  inline void addSeed(const unsigned ieta, const unsigned iphi, const double E )
   {
      m_seeds[m_nSeeds++] = MonoEcalSeed(m_seedLength,ieta,iphi,E);
   }

  // merge seeds
  void mergeSeeds(const EBmap &);

  // center hot spot of seeds in strip
  void centerSeeds(const EBmap &);

  // length of the seeds to find (number of cells in eta)
  unsigned m_seedLength;

  // energy threshold of seed
  double m_threshold;
 
  // number of seeds found 
  unsigned m_nSeeds;

  // seed array
  unsigned m_maxSeeds;
  std::vector<MonoEcalSeed> m_seeds;
  unsigned m_stageSize;
  std::vector<unsigned> m_stage;

  // calorimetry geometry
  const CaloSubdetectorGeometry *m_geom;

};

// ----------------------------------------------------------
// Simply extends seeds in eta strips
class SimplePathFinder {

};


// ---------------------------------------------------------
// Build clusters from seeds
class ClusterBuilder {

public:
  inline ClusterBuilder(): m_nClusters(0U) { m_clusters.resize(50U); }

  inline virtual ~ClusterBuilder () { }

  // build clusters around seeds
  void buildClusters(unsigned,const MonoEcalSeed *, const EBmap &);

  // accessor methods
  inline const unsigned nClusters() const { return m_nClusters; }
  inline const MonoEcalCluster * clusters() const { return &m_clusters[0]; }

private:

  // number of clusters
  unsigned m_nClusters;

  std::vector<MonoEcalCluster> m_clusters;

};



// ---------------------------------------------------------
// MonoEcalObs0 
class MonoEcalObs0 {

public:

  inline MonoEcalObs0(const edm::ParameterSet &ps) 
    :m_seedLength(ps.getParameter<unsigned>("StripSeedLength") )
    ,m_threshold(ps.getParameter<double>("SeedThreshold") )
    {
      m_seedFinder = StripSeedFinder(m_seedLength,m_threshold,m_ecalMap.nCells());
      m_seedFinder.initialize();
    } 

  inline virtual ~MonoEcalObs0()
    {
     // if ( m_ecalMap ) delete m_ecalMap;
    }


  // calculate observable method
  double calculate(const edm::EventSetup &,const edm::Event &);

  // accessor methods
  inline const StripSeedFinder & finder() const { return m_seedFinder; }
  inline const ClusterBuilder & clusterBuilder() const { return m_clusterBuilder; }
  inline const EBmap & ecalMap() const { return m_ecalMap; }

private:

  // -- private member functions
  
  // return the expected energy in bin i 
  double eBarI(unsigned i);

  // find Beta ij
  double betaij(unsigned i, unsigned j);

  // calculate M ij
  double mij(unsigned i, unsigned j);

  // load H matrix lookup tables
  void loadHMatTables();


  // -- private member data
  // seed length
  unsigned m_seedLength;
  // seed threshold
  double m_threshold;
  

  // ecalMap
  EBmap m_ecalMap;

  // the seed finder
  StripSeedFinder m_seedFinder;

  // the cluster builder
  ClusterBuilder m_clusterBuilder;

  // H matrix look up table map
  std::map<ClustCategorizer,std::vector<double> > m_hMatMap;

};  


//
// Calibrator class for the Monopole Ecal observable
class MonoEcalObs0Calibrator {

public:

  inline MonoEcalObs0Calibrator(const edm::ParameterSet &ps)
    :m_seedLength(ps.getParameter<unsigned>("StripSeedLength"))
    ,m_threshold(ps.getParameter<double>("SeedThreshold"))
    ,m_calibName(ps.getParameter<std::string>("CalibrationName")) 
    ,m_wsSize(50U)
    {
      m_seedFinder = StripSeedFinder(m_seedLength,m_threshold,m_ecalMap.nCells());
      m_seedFinder.initialize();

      m_workspace.resize(m_wsSize);

      double pars[3] = {1.,0.4,0.4};
      m_functor.setParameters(3,pars);
    }

  inline virtual ~MonoEcalObs0Calibrator() { }

  // calculate M_ij^n called for every event
  void calculateMijn(const edm::EventSetup &, const edm::Event &);

  // computed H_ij at end of run
  void calculateHij();

  // dump calibration
  void dumpCalibration();

  // set the functor parameters
  inline void setClusterParameters(const unsigned N, const double * pars)
    {
      m_functor.setParameters(N,pars);
    }

  // set the functor
  inline void setFunctor(const EnergyFlowFunctor &functor)
    {
      m_functor = EnergyFlowFunctor(functor);
    }

private:

  // -- private member functions
  // compute M_ij at end of run in order to find H_ij
  void computeMij();

  // -- private member data
  std::map<ClustCategorizer,std::vector<double> > m_hij;
  std::map<ClustCategorizer,std::vector<std::vector<double> > > m_Mijn;

  // seed length
  unsigned m_seedLength;
  // seed threshold
  double m_threshold;
  // output calibration file name
  std::string m_calibName;

  EBmap m_ecalMap;
  StripSeedFinder m_seedFinder;
  ClusterBuilder m_clusterBuilder;

  const unsigned m_wsSize;
  std::vector<double> m_workspace;

  // energy flow 
  EnergyFlowFunctor m_functor;


};  // end calibrator class


}  // end Mono namespace

#endif
