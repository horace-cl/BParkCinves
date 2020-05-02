#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "CommonTools/CandUtils/interface/Booster.h"
#include "Math/GenVector/Boost.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"


#include <vector>
#include <memory> 
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"
#include <Math/VectorUtil.h>
 
class PrimaryVertexSelector : public edm::global::EDProducer<> {

public:

  explicit PrimaryVertexSelector(const edm::ParameterSet &cfg):
    bMesons_(consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("bMesons") )),
    trgMuonToken_(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("trgMuon"))),
    vertexSrc_( consumes<reco::VertexCollection> ( cfg.getParameter<edm::InputTag>( "vertexCollection" ) ) ),
    dzCut_(cfg.getParameter<double>("dzCut")){
		produces<nanoaod::FlatTable>();
    }

  ~PrimaryVertexSelector() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  //virtual void produce(edm::Event&, const edm::EventSetup&);

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}

  
private:
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> bMesons_;
  const edm::EDGetTokenT<pat::MuonCollection> trgMuonToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_; 
  const double dzCut_;
};

//void PrimaryVertexSelector::produce(edm::Event& evt, const edm::EventSetup& iSetup) {
void PrimaryVertexSelector::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

	//input
	edm::Handle<pat::CompositeCandidateCollection> bMesons;
	evt.getByToken(bMesons_, bMesons);

	edm::Handle<reco::VertexCollection> vertexHandle;
	evt.getByToken(vertexSrc_, vertexHandle);

	edm::Handle<pat::MuonCollection> trgMuons;
	evt.getByToken(trgMuonToken_, trgMuons);


	std::vector<float> vx,vy,vz, dzTrgMu;// ,dlenSig,pAngle;

	dzTrgMu.push_back(fabs(vertexHandle->front().position().z()-trgMuons->front().vz()));
	vx.push_back(vertexHandle->front().position().x());
	vy.push_back(vertexHandle->front().position().y());
	vz.push_back(vertexHandle->front().position().z());

	int nv = 0;
	for(size_t i=1;i < (*vertexHandle).size(); i++){
		for (const pat::Muon & trgmu : *trgMuons){
			if(nv==4) break;
			if(fabs((*vertexHandle)[i].position().z()-trgmu.vz())<dzCut_) {
				dzTrgMu.push_back(fabs((*vertexHandle)[i].position().z()-trgmu.vz()));
				vx.push_back((*vertexHandle)[i].position().x());
				vy.push_back((*vertexHandle)[i].position().y());
				vz.push_back((*vertexHandle)[i].position().z());
				nv++;
			}
		}
	}



	// output
    auto pvTable = std::make_unique<nanoaod::FlatTable>(dzTrgMu.size(),"PV",false);

    pvTable->addColumn<float>("dzTrgMu",dzTrgMu,"abs(vz-trgMuon.vz)",nanoaod::FlatTable::FloatColumn,10);
    pvTable->addColumn<float>("vx",vx,"vx in cm",nanoaod::FlatTable::FloatColumn,10);
    pvTable->addColumn<float>("vy",vy,"vy in cm",nanoaod::FlatTable::FloatColumn,10);
    pvTable->addColumn<float>("vz",vz,"vz in cm",nanoaod::FlatTable::FloatColumn,10);


	// std::cout << "B Mesons:  "<< bMesons->size() << "\n";
	// std::cout << "Vertices:  "<< vertexHandle->size() << "\n";
	// std::cout << "Trg muons: "<< trgMuons->size() << "\n";


	evt.put(std::move(pvTable));

}


DEFINE_FWK_MODULE(PrimaryVertexSelector);