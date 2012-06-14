#ifndef Monopoles_MonoAlgorithms_NPVHelper_h
#define Monopoles_MonoAlgorithms_NPVHelper_h


#include "FWCore/Framework/interface/Event.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"



namespace Mono {

inline unsigned int getNPV(const edm::Event & iEvent, const edm::EventSetup &iSetup )
{
  // PrimaryVertex analysis
  edm::InputTag m_PVTag = edm::InputTag("offlinePrimaryVertices","");

  edm::Handle<reco::VertexCollection> handlePV;
  iEvent.getByLabel(m_PVTag,handlePV);

  int totalNPV = 0.;

  reco::VertexCollection::const_iterator pv = handlePV->begin();
  for ( ; pv != handlePV->end(); pv++ ) {
    if ( !pv->isFake() && pv->ndof() > 4.0 ) {
      ++totalNPV;
    }
  }

  return totalNPV;
}


} // end Mono namespace


#endif
