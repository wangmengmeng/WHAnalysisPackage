#ifndef RecoTauTag_RecoTau_RecoTauVertexAssociator_h
#define RecoTauTag_RecoTau_RecoTauVertexAssociator_h

/* RecoTauVertexAssociator
 *
 * Authors: Evan K. Friis, Christian Veelken, UC Davis
 *          Michalis Bachtis, UW Madison
 *
 * The associatedVertex member function retrieves the vertex from the event
 * associated to a given tau.  This class is configured using a cms.PSet.
 *
 * The required arguments are:
 *  o primaryVertexSrc - InputTag with the vertex collection
 *  o useClosestPV - Use the "closest to lead track in z" to find the vertex.
 *
 * The setEvent method must be called at least once per event.
 *
 */

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/Handle.h"

// Forward declarations
namespace edm {
  class ParameterSet;
  class Event;
}

namespace reco {
  class PFTau;
  class PFJet;
}

namespace reco { namespace tau {

class RecoTauVertexAssociator {
  public:
    enum Algorithm {
      kHighestPtInEvent,
      kClosestDeltaZ,
      kHighestWeigtForLeadTrack
    };

    RecoTauVertexAssociator (const edm::ParameterSet& pset);
    virtual ~RecoTauVertexAssociator (){}
    /// Get the primary vertex associated to a given jet. Returns a null Ref if
    /// no vertex is found.
    reco::VertexRef associatedVertex(const PFJet& tau) const;
    /// Convenience function to get the PV associated to the jet that
    /// seeded this tau.
    reco::VertexRef associatedVertex(const PFTau& tau) const;
    /// Load the vertices from the event.
    void setEvent(const edm::Event& evt);
  private:
    std::vector<reco::VertexRef> vertices_;
    edm::InputTag vertexTag_;
    Algorithm algo_;
};

} /* tau */ } /* reco */

#endif /* end of include guard: RecoTauTag_RecoTau_RecoTauVertexAssociator_h */
