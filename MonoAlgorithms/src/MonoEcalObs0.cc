
#include "Monopoles/MonoAlgorithms/interface/MonoEcalObs0.h"

#include <iostream>

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"


namespace Mono {

// --------------------------- EBmap member functions ---------------------------------

void EBmap::fillMap(const edm::Event &ev)
{
  // clear the maps
  clear();

  // fill Ecal Map with Rec Hit energy
  edm::InputTag tagEcalEB_RecHits("ecalRecHit","EcalRecHitsEB");
  edm::Handle<EBRecHitCollection> ecalRecHits;
  ev.getByLabel(tagEcalEB_RecHits,ecalRecHits);

  double maxE = 0.;

  const unsigned nHits = ecalRecHits->size();
  for ( unsigned i=0; i != nHits; i++ ) {
    EBDetId detId( (*ecalRecHits)[i].id() );
    const unsigned loc = findBin(detId);
    const double energy = (*ecalRecHits)[i].energy();
    if ( energy > maxE ) maxE = energy;
    m_ecalMap[loc] = energy;
    m_ecalTMap[loc] = (*ecalRecHits)[i].time();
  }  

}

void EBmap::clear()
{
  for ( unsigned i=0; i != m_nCells; i++ ) {
    m_ecalMap[i] = 0.;
    m_ecalTMap[i] = 0.;
  }
}

void EBmap::constructGeo(const edm::EventSetup & es)
{
  // get a handle on the Ecal Barrel subdetector geometry
  edm::ESHandle<CaloGeometry> calo;
  es.get<CaloGeometryRecord>().get(calo);
  const CaloGeometry * caloGeo = (const CaloGeometry*)calo.product();
  m_geom = caloGeo->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
}

const unsigned EBmap::findBin( const EBDetId &did ) const
{
  
  const CaloCellGeometry *cell = m_geom->getGeometry( did );
  const double eta = cell->getPosition().eta();
  const double phi = cell->getPosition().phi();
  return findBinEtaPhi(eta,phi);

}

const unsigned EBmap::findBinEtaPhi(const double eta, const double phi ) const
{

  if ( eta > m_maxEta || eta < m_minEta || phi > m_maxPhi || phi < m_minPhi )
    throw cms::Exception("eta,phi out of range") << "The eta and/or phi given to EBmap::findBinEtaPhi is out of range";

  assert(phi <= M_PI && phi >= -M_PI );
 
  unsigned phiBin = (phi - m_minPhi)/m_phiWidth;
  unsigned etaBin = (eta-m_minEta)/m_etaWidth;

  const unsigned bin = phiBin*m_nEta+etaBin;
  assert( bin < m_nCells );

  return bin;
 
} 


// --------------------------- StripSeedFinder member functions ------------------------

void StripSeedFinder::constructGeo(const edm::EventSetup &es)
{

  // get a handle on the Ecal Barrel subdetector geometry
  edm::ESHandle<CaloGeometry> calo;
  es.get<CaloGeometryRecord>().get(calo);
  const CaloGeometry * caloGeo = (const CaloGeometry*)calo.product();
  m_geom = caloGeo->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);


}


bool StripSeedFinder::find(const edm::Event &ev, const EBmap &ecalMap) 
{

  const unsigned nCells = ecalMap.nCells();
  assert( m_seeds.size() == nCells );

  // clear old seeds if any
  clear();

  const unsigned nEta = ecalMap.nEta();
  const unsigned nPhi = ecalMap.nPhi();


  // start looking for seeds
  const unsigned etaSearch = nEta - m_seedLength + 1U;
  for ( unsigned i=0; i != nPhi; i++ ) {
    for ( unsigned j=0; j != etaSearch; j++ ) {
      const unsigned loc = i*nEta+j;
      double energySum = ecalMap[loc];
      for ( unsigned s=1; s != m_seedLength; s++ ) 
	energySum += ecalMap[loc+s];

      if ( energySum > m_threshold ) addSeed(j,i,energySum);
    }
  } 

  // search in phi direction (cross-check of signature)
  /*const unsigned phiSearch = nPhi - m_seedLength + 1U;
  for ( unsigned i=0; i != phiSearch; i++ ) {
    for ( unsigned j=0; j != nEta; j++ ) {
      const unsigned loc = i*nEta+j;
      double energySum = ecalMap[loc];
      for ( unsigned s=1; s != m_seedLength; s++ )
	energySum += ecalMap[(i+s)*nEta+j];

      if ( energySum > m_threshold ) addSeed(j,i,energySum);
    }
  }*/

  if ( m_nSeeds > 1 ) mergeSeeds(ecalMap);
  centerSeeds(ecalMap);
  if ( m_nSeeds > 1 ) mergeSeeds(ecalMap);

  return 1;
}



void StripSeedFinder::mergeSeeds(const EBmap &map)
{

  const unsigned nEta = map.nEta();

  // stage of indeces into m_seeds
  unsigned m_stageSize = 0U;

  // cycle over seeds
  // if no overlap of seed in stage, add to stage
  // else merge into seed in stage 
  m_stage[m_stageSize++] = 0U;
  for ( unsigned s = 1U; s != m_nSeeds; s++ ) {
    const MonoEcalSeed & seed = m_seeds[s];
    bool didMerge = false;

    for ( unsigned t = 0; t != m_stageSize; t++ ) {
      const unsigned seedIndex = m_stage[t];

      if ( overlapInEta(s,seedIndex) ) {
       	const unsigned sEta = seed.ieta();
     	const unsigned sPhi = seed.iphi();
        const MonoEcalSeed & teed = m_seeds[seedIndex];
    	const unsigned tEta = teed.ieta();
	const unsigned tPhi = teed.iphi();
      	const unsigned sLength = seed.seedLength();
	const unsigned tLength = seed.seedLength();

	if ( sEta > tEta ) {
	  const unsigned newLength = sLength + tLength - (tEta+tLength-sEta);
	  double newEnergy = teed.energy();
	  for ( unsigned i=tLength; i != newLength; i++ )
	    newEnergy += map[tPhi*nEta+tEta+i];
	  m_seeds[seedIndex] = MonoEcalSeed(newLength,tEta,tPhi,newEnergy);
	} else if ( tEta > sEta ) {
	  const unsigned newLength = tLength + sLength - (sEta+sLength-tEta);	
	  double newEnergy = seed.energy();
	  for ( unsigned i=sLength; i != newLength; i++ )
	    newEnergy += map[sPhi*nEta+sEta+i];
	  m_seeds[seedIndex] = MonoEcalSeed(newLength,sEta,sPhi,newEnergy);
	}

      	didMerge = true;	
      }
    }

    if ( !didMerge ) m_stage[m_stageSize++] = s;

  }

  // copy stage into m_seeds
  for ( unsigned s=0; s != m_stageSize; s++ ) 
    m_seeds[s] = m_seeds[m_stage[s]];
  m_nSeeds = m_stageSize;

  m_stageSize = 0U;

}


void StripSeedFinder::centerSeeds(const EBmap &map)
{

  const unsigned nEta = map.nEta();

  // cycle over seeds
  for ( unsigned s=0; s != m_nSeeds; s++ ) {
    const MonoEcalSeed & seed = m_seeds[s];
    const unsigned length = seed.seedLength();
    const unsigned loc = seed.iphi()*nEta+seed.ieta();
    unsigned hotSpot = 0U;
    double hotPeak = 0.;
    for ( unsigned i=0U; i != length; i++ ) {
      const double E = map[loc+i];
      if ( E > hotPeak ) {
	hotPeak = E;
	hotSpot = i;
      }
    }

    // is the hotspot centered?
    // if not...buffer the short side 
    const int middle = length/2;
    const int dist = hotSpot - middle;
    unsigned newLength = length;
    unsigned newEta = seed.ieta();
    if ( length % 2 == 0 ) {
      if ( dist > 1 ) {
	 newLength += dist;
      } else if ( dist < -1 ) {
	newLength -= dist;
	(unsigned)abs(dist) < newEta ? newEta -= (unsigned)abs(dist) : newEta = 0U;
      }
      else continue; // no need to do anything to the seed since the hot spot is centered
    } else if ( dist != 0 ) { 
      newLength += abs(dist);
      if ( dist < 0 ) (unsigned)abs(dist) < newEta ? newEta -= (unsigned)abs(dist) : newEta = 0U;
    } else continue; // no need to do anything to the seed since the hot spot is centered

    if ( newEta+newLength >= nEta ) newLength = nEta-newEta-1; 

    assert( newEta < nEta );
    assert( newLength + newEta < nEta );

    // find the new energy
    const unsigned newLoc = seed.iphi()*nEta+newEta;
    double newE = 0.;
    for ( unsigned i=0; i != newLength; i++ )
      newE += map[loc+i];
   
    // update seed list   
    m_seeds[s] = MonoEcalSeed(newLength,newEta,seed.iphi(),newLength);
    
  }

}

// --------------------------- ClusterBuilder member functions ---------------------

void ClusterBuilder::buildClusters(const unsigned nSeeds, const MonoEcalSeed *seeds, const EBmap &map)
{

  assert( seeds );
  m_nClusters = 0U;
  if ( nSeeds > m_clusters.size() ) m_clusters.resize(nSeeds);

  const unsigned nEta = map.nEta();
  const unsigned nPhi = map.nPhi();

  // add N eta strips to each side of the seed
  const unsigned N=2;

  for ( unsigned s=0; s != nSeeds; s++ ) {
    const MonoEcalSeed & seed = seeds[s];
    const unsigned sEta = seed.ieta();
    const unsigned sPhi = seed.iphi();
    const double sEnergy = seed.energy();
    const unsigned length = seed.seedLength();

    double energy = sEnergy;
    // add below in phi
    for ( unsigned i=N; i != 0; i-- ) {
      int newPhiS = (int)sPhi-i;
      unsigned newPhi = 0U;
      if ( newPhiS < 0 ) newPhi = nPhi+newPhiS;
      else newPhi = newPhiS;
      assert( newPhi < nPhi );
      for ( unsigned i=0; i != length; i++ ) {
	energy += map[newPhi*nEta+sEta+i];
      } 
    }
    // add above in phi
    for ( unsigned i=1; i != N+1; i++ ) {
      unsigned newPhi = sPhi+i;
      if ( newPhi >= nPhi ) newPhi = newPhi-nPhi;
      assert ( newPhi < nPhi );
      for ( unsigned i=0; i != length; i++ ) {
	energy += map[newPhi*nEta+sEta+i];
      }
    }
 
    m_clusters[m_nClusters++] = MonoEcalCluster(length,2*N+1U,sEta,sPhi,energy,seed); 
  }
  
}


// --------------------------- MonoEcalObs0 member functions ------------------------

double MonoEcalObs0::calculate(const edm::EventSetup &es, const edm::Event &ev)
{
  m_ecalMap.constructGeo(es);
  m_seedFinder.constructGeo(es);

  // fill ecal map
  m_ecalMap.fillMap(ev);

  // run seed finder
  m_seedFinder.find(ev,m_ecalMap);

  // run the cluster builders
  m_clusterBuilder.buildClusters(m_seedFinder.nSeeds(),m_seedFinder.seeds(),m_ecalMap);

  return 0.;
}

double MonoEcalObs0::eBarI(const unsigned i)
{
  return 0.;
}

double MonoEcalObs0::betaij(const unsigned i, const unsigned j)
{
  return 0;
}

double MonoEcalObs0::mij(const unsigned i,const unsigned j)
{
  return 0.;
}

void loadHMatTables()
{ }

// ------------------------- MonoEcalOb0Calibrator ---------------------------------

void MonoEcalObs0Calibrator::calculateMijn(const edm::EventSetup &es, const edm::Event &ev)
{
  m_ecalMap.constructGeo(es);
  m_seedFinder.constructGeo(es);

  // fill ecal map
  m_ecalMap.fillMap(ev);

  // run seed finder
  m_seedFinder.find(ev,m_ecalMap);

  // run the cluster builders
  m_clusterBuilder.buildClusters(m_seedFinder.nSeeds(),m_seedFinder.seeds(),m_ecalMap);


  const unsigned nEta = m_ecalMap.nEta();

  // cycle over clusters and fill appropriate map
  const unsigned nClusters = m_clusterBuilder.nClusters();
  const MonoEcalCluster * clusters = m_clusterBuilder.clusters();
  for ( unsigned i=0; i != nClusters; i++ ) {
    const unsigned length = clusters[i].clusterLength();
    const unsigned width = clusters[i].clusterWidth();

    const ClustCategorizer cat(length,width);

    // check if catergory exists in map
    std::map<ClustCategorizer,std::vector<std::vector<double> > >::iterator iter = m_Mijn.find(cat);
    if ( iter == m_Mijn.end() ) {
      std::vector<std::vector<double> > & vec = m_Mijn[cat];
      vec.resize( width*width*length*length );
    }

    std::vector<std::vector<double> > & vec = m_Mijn[cat];

    // clear workspace
    for ( unsigned j=0; j != m_wsSize; j++ ) m_workspace[j]=0.;


    const unsigned ieta = clusters[i].ieta();
    const unsigned iphi = clusters[i].iphi();
    const unsigned size = width*length;
    const double eTot = clusters[i].clusterEnergy();
    assert( size <= m_wsSize );
    // fill workspace
    for ( unsigned k=0; k != width; k++ ) {
      unsigned newPhi = 0U;  
      unsigned newPhiDiff = 0U;
      if ( k < width/2 )  newPhiDiff = width/2-k;
      else newPhiDiff = k-width/2;
      
      for ( unsigned j=0; j != length; j++ ) {
	const unsigned loc = (iphi+k-width/2)*nEta+ieta+j; 
     	const unsigned num = k*length+j;
	m_workspace[num] = m_ecalMap[loc]/eTot-m_functor(j-length/2,k-width/2);
      }
    }
  

    // calculate M_ij^n
    assert( size*size == vec.size() );
    for ( unsigned j=0; j != size; j++ ) {
      const double E1 = m_workspace[j];
      for ( unsigned k=0; k != size; k++ ) {
	vec[j*size+k].push_back(E1*m_workspace[k]);
      }
    }
	 
  }

}


void MonoEcalObs0Calibrator::calculateHij()
{ }

void MonoEcalObs0Calibrator::dumpCalibration()
{ }

void MonoEcalObs0Calibrator::computeMij()
{ }


} //end Mono namespace
