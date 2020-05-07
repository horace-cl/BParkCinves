//Crearemos el Builder desde cero


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <TLorentzVector.h>
#include "helper.h"


//Esta clase debe Heredar de un EDProducer
class BToKmm : public edm::EDProducer{
	//Solo definimos los elementos
	//Elementos publicos de la clase
	public:


		//Constructor 
		//Debe tener el mismo nombre que la clase y recibir ese argumento
		explicit BToKmm(const edm::ParameterSet &iConfig);

		//Destructor
		~MuonTriggerSelector() override {};


	//Elementos privados de la clase
	private:

		//Metodo produce, el cual se ejecuta durante el "RunTime"
		virtual void produce(edm::Event&, const edm::EventSetup&);

		//Objetos que recibe la clase a travez del configurador
		//TODO
	    edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;
	    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
	    edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
	    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
	    edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;

	    
	    // Parametros para la seleccion de Muones
	    const double maxdR_; //for trigger match
	    const double dzTrg_cleaning_; // selects primary vertex
	    const double ptMin_;          // min pT in all muons for B candidates
	    const double absEtaMax_;      //max eta ""
	    const bool softMuonsOnly_;    //cuts muons without soft ID

}



