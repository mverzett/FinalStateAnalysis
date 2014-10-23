#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "FinalStateAnalysis/PatTools/interface/METhelpers.h"

class PATMETUserCandEmbedder : public edm::EDProducer {
  public:
    typedef reco::LeafCandidate ShiftedCand;
    typedef std::vector<ShiftedCand> ShiftedCandCollection;
    typedef reco::Candidate::LorentzVector LorentzVector;
    typedef std::vector<edm::InputTag> VInputTag;
    typedef std::vector<edm::ParameterSet> VPSet;

    PATMETUserCandEmbedder(const edm::ParameterSet& pset);
    virtual ~PATMETUserCandEmbedder(){}
    void produce(edm::Event& evt, const edm::EventSetup& es);
    void embedShift(edm::Event& evt, pat::MET& met, 
		const std::string& embedName, const std::string& branchname, 
		const reco::Candidate::LorentzVector p4);
  private:
    edm::InputTag metSrc_;
    VPSet userCands_;
};

PATMETUserCandEmbedder::PATMETUserCandEmbedder(
    const edm::ParameterSet& pset):
  metSrc_(pset.getParameter<edm::InputTag>("src")),
  userCands_(pset.getParameter<VPSet>("userCands")){
  produces<pat::METCollection>();
  for(VPSet::const_iterator pset = userCands_.begin(); pset != userCands_.end(); ++pset){
    std::string branch_name = pset->getParameter<std::string>("branchName");
    produces<ShiftedCandCollection>(branch_name);      
  }
}

void PATMETUserCandEmbedder::embedShift(edm::Event& evt, pat::MET& met, 
		const std::string& embedName, const std::string& branchname, 
		const reco::Candidate::LorentzVector p4){
  typedef reco::LeafCandidate ShiftedCand;
  typedef std::vector<ShiftedCand> ShiftedCandCollection;
  typedef edm::OrphanHandle<ShiftedCandCollection> PutHandle;
  typedef reco::CandidatePtr CandidatePtr;

  std::auto_ptr<ShiftedCandCollection> output(new ShiftedCandCollection);
  ShiftedCand newCand(met);
  newCand.setP4(p4);
  output->push_back(newCand);
  PutHandle outputH = evt.put(output, branchname);
  met.addUserCand(embedName, CandidatePtr(outputH, 0));
}


void PATMETUserCandEmbedder::produce(edm::Event& evt, const edm::EventSetup& es) {


  edm::Handle<pat::METCollection> mets;
  evt.getByLabel(metSrc_, mets);

  assert(mets->size() == 1);

  const pat::MET& inputMET = mets->at(0);
  pat::MET outputMET = inputMET;

  for(VPSet::const_iterator pset = userCands_.begin(); pset != userCands_.end(); ++pset){
    std::string branch_name = pset->getParameter<std::string>("branchName");
    std::string embed_name  = pset->getParameter<std::string>("embedName");
    edm::InputTag cand_src  = pset->getParameter<edm::InputTag>("src");

    edm::Handle< edm::View<reco::MET> > candidate;
    evt.getByLabel(cand_src, candidate);
    
    embedShift(evt, outputMET, embed_name, branch_name, candidate->at(0).p4());
  }

  std::auto_ptr<pat::METCollection> outputColl(new pat::METCollection);
  outputColl->push_back(outputMET);

  evt.put(outputColl);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATMETUserCandEmbedder);
