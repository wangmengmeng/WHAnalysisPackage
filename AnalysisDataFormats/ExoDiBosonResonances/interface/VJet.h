#ifndef _AnalysisDataFormats_ExoDiBosonResonances_VJet_H_
#define _AnalysisDataFormats_ExoDiBosonResonances_VJet_H_

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "CMGTools/External/interface/PileupJetIdentifier.h"
#include "AnalysisDataFormats/CMGTools/interface/TriBool.h"
#include "AnalysisDataFormats/CMGTools/interface/UnSet.h"

#include "AnalysisDataFormats/CMGTools/interface/BaseJet.h"
#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"
#include "AnalysisDataFormats/CMGTools/interface/PFJetComponent.h"

#include <vector>

namespace cmg {

  //forward def needed
  class VJet;

  /** Class representing V-tagged jets in the CMG framework.
      
  In addition to the attributes from the BaseJet and PFJet mother classes, 
  VJets contains...

  */


  class VJet : public PFJet {
  public:

    ///number of subjets in pruned jet
    //   static const unsigned NSUBJETS;
    //empty constructor
    VJet() :  qjet_(-1.0), tau1_(0.0), tau2_(-1.0), tau3_(99.0),tau4_(99.0), tau5_(99.0),mdrop_(-99.0), prunedMass_(-111.0), nfatjetL_(-1.), nfatjetM_(-1.),nfatjetT_(-1.), nsubjetL_(-1.), nsubjetM_(-1.),nsubjetT_(-1.),subjet1_phi_(10.),subjet1_eta_(10.),subjet2_phi_(10.),subjet2_eta_(10.), dR_subjet_(10.0), subjet_flavor_(30){};

    //constructor
      VJet(const value& j): PFJet(j), qjet_(-1.0), tau1_(0.0), tau2_(-1.0), tau3_(99.0),tau4_(99.0), tau5_(99.0),mdrop_(-99.0), prunedMass_(-111.0),nfatjetL_(-1.),     nfatjetM_(-1.),nfatjetT_(-1.), nsubjetL_(-1.), nsubjetM_(-1.),nsubjetT_(-1.),subjet1_phi_(10.),subjet1_eta_(10.),subjet2_phi_(10.),subjet2_eta_(10.),dR_subjet_(10.0),subjet_flavor_(30){};

    virtual ~VJet(){}

    float qjet() const {return qjet_;}
    float tau1() const {return tau1_;}
    float tau2() const {return tau2_;}
    float tau3() const {return tau3_;}
    float ntau21() const {return tau2_/tau1_;}
    float ntau32() const {return tau3_/tau2_;}
    float ntau31() const {return tau3_/tau1_;}
// add tau4 tau5 tat41 etc
    float tau4() const {return tau4_;}
    float tau5() const {return tau5_;}
    float ntau41() const {return tau4_/tau1_;}
    float ntau42() const {return tau4_/tau2_;}
    float ntau43() const {return tau4_/tau3_;}
    float ntau51() const {return tau5_/tau1_;}
    float ntau52() const {return tau5_/tau2_;}
    float ntau53() const {return tau5_/tau3_;}
    float ntau54() const {return tau5_/tau4_;}

    float mdrop() const {return mdrop_;}
    float prunedMass() const {return prunedMass_;}
// add fat & sub b jet here for WH
//    float subjet_flavor() const {return subjet_flavor_;}
    int subjet_flavor() const {return subjet_flavor_;}
    float subjet1_phi() const {return subjet1_phi_;}
    float subjet1_eta() const {return subjet1_eta_;}
    float subjet2_phi() const {return subjet2_phi_;}
    float subjet2_eta() const {return subjet2_eta_;}
    float dR_subjet() const {return dR_subjet_;}
    float nsubjetL() const {return nsubjetL_;}
    float nsubjetM() const {return nsubjetM_;}
    float nsubjetT() const {return nsubjetT_;}
    float nfatjetL() const {return nfatjetL_;}
    float nfatjetM() const {return nfatjetM_;}
    float nfatjetT() const {return nfatjetT_;}

    //return pointer to matched pruned jet
    pat::JetPtr const* prunedJetPtr() const{
      return &prunedJetPtr_;
    }

    void setPrunedJetPtr(const pat::JetPtr& ptr){
      prunedJetPtr_=ptr;
    }

    friend class VJetFactory;

  private:

    float qjet_; //qjet volatility
    float tau1_, tau2_,tau3_, tau4_, tau5_;
    float mdrop_, prunedMass_;
// add fat & sub b jet here for WH
//    float subjet_flavor_;
    int subjet_flavor_;
    float subjet1_phi_, subjet1_eta_, subjet2_phi_,subjet2_eta_;
    float dR_subjet_;
    float nsubjetL_, nsubjetM_, nsubjetT_;
    float nfatjetL_, nfatjetM_, nfatjetT_;

    pat::JetPtr prunedJetPtr_;

  };
}


#endif
