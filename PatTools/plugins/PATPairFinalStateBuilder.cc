#include "FinalStateAnalysis/PatTools/plugins/PATPairFinalStateBuilderT.h"
#include "FinalStateAnalysis/DataFormats/interface/PATDiLeptonFinalStates.h"

typedef PATPairFinalStateBuilderT<PATElecElecFinalState> PATElecElecFinalStateProducer;
typedef PATPairFinalStateBuilderT<PATElecMuFinalState> PATElecMuFinalStateProducer;
typedef PATPairFinalStateBuilderT<PATElecTauFinalState> PATElecTauFinalStateProducer;
typedef PATPairFinalStateBuilderT<PATElecPhoFinalState> PATElecPhoFinalStateProducer;
typedef PATPairFinalStateBuilderT<PATMuMuFinalState> PATMuMuFinalStateProducer;
typedef PATPairFinalStateBuilderT<PATMuTauFinalState> PATMuTauFinalStateProducer;
typedef PATPairFinalStateBuilderT<PATMuPhoFinalState> PATMuPhoFinalStateProducer;
typedef PATPairFinalStateBuilderT<PATTauTauFinalState> PATTauTauFinalStateProducer;
typedef PATPairFinalStateBuilderT<PATTauPhoFinalState> PATTauPhoFinalStateProducer;
typedef PATPairFinalStateBuilderT<PATPhoPhoFinalState> PATPhoPhoFinalStateProducer;
typedef PATPairFinalStateBuilderT<PATMuJetFinalState> PATMuJetFinalStateProducer;
typedef PATPairFinalStateBuilderT<PATElecJetFinalState> PATElecJetFinalStateProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATElecElecFinalStateProducer);
DEFINE_FWK_MODULE(PATElecMuFinalStateProducer);
DEFINE_FWK_MODULE(PATElecTauFinalStateProducer);
DEFINE_FWK_MODULE(PATElecPhoFinalStateProducer);
DEFINE_FWK_MODULE(PATMuMuFinalStateProducer);
DEFINE_FWK_MODULE(PATMuTauFinalStateProducer);
DEFINE_FWK_MODULE(PATMuPhoFinalStateProducer);
DEFINE_FWK_MODULE(PATTauTauFinalStateProducer);
DEFINE_FWK_MODULE(PATTauPhoFinalStateProducer);
DEFINE_FWK_MODULE(PATPhoPhoFinalStateProducer);
DEFINE_FWK_MODULE(PATMuJetFinalStateProducer);
DEFINE_FWK_MODULE(PATElecJetFinalStateProducer);

