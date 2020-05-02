// class to produce 2 pat::MuonCollections
// one matched to the Park triggers
// another fitered wrt the Park triggers


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

using namespace std;

constexpr bool debug = false;

class MuonTriggerSelector : public edm::EDProducer {
    
public:
    
    explicit MuonTriggerSelector(const edm::ParameterSet &iConfig);
    
    ~MuonTriggerSelector() override {};
    
    
private:

    virtual void produce(edm::Event&, const edm::EventSetup&);

    edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
    edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;

    //for trigger match
    const double maxdR_;

    //for filter wrt trigger
    const double dzTrg_cleaning_; // selects primary vertex

    const double ptMin_;          // min pT in all muons for B candidates
    const double absEtaMax_;      //max eta ""
    const bool softMuonsOnly_;    //cuts muons without soft ID
}; 


MuonTriggerSelector::MuonTriggerSelector(const edm::ParameterSet &iConfig):
  muonSrc_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  vertexSrc_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) ), 
  maxdR_(iConfig.getParameter<double>("maxdR_matching")),
  dzTrg_cleaning_(iConfig.getParameter<double>("dzForCleaning_wrtTrgMuon")),
  ptMin_(iConfig.getParameter<double>("ptMin")),
  absEtaMax_(iConfig.getParameter<double>("absEtaMax")),
  softMuonsOnly_(iConfig.getParameter<bool>("softMuonsOnly"))
{
  // produce 2 collections: trgMuons (tags) and SelectedMuons (probes & tags if survive preselection cuts)
    produces<pat::MuonCollection>("trgMuons"); 
    produces<pat::MuonCollection>("SelectedMuons");
    produces<TransientTrackCollection>("SelectedTransientMuons");  
}
 
 

void MuonTriggerSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
    edm::ESHandle<MagneticField> bFieldHandle;
    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

    edm::Handle<reco::VertexCollection> vertexHandle;
    iEvent.getByToken(vertexSrc_, vertexHandle);
    const reco::Vertex & PV = vertexHandle->front();

    if(debug) std::cout << " MuonTriggerSelector::produce " << std::endl;

    // Not so clear
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMPathsAndTriggerBits
    edm::Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken(triggerBits_, triggerBits);
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

    std::vector<pat::TriggerObjectStandAlone> triggeringMuons;

    //taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#Trigger
    edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
    iEvent.getByToken(triggerObjects_, triggerObjects);
    if(debug) std::cout << "\n TRIGGER OBJECTS " << std::endl;


    int int_obj = 0;
    // std::cout << "\n\n\n------------>>>>>>NEW RECORD NEW RECORD NEW RECORD NEW RECORD"<<"\n";
    
 
    //Itera sobre los objetos 
    //corrobora que al menos haya un filterID = 83
    //corrobora que al menos uno tenga 'hltL3' y 'Park'
    // el objeto que satisfaga estas dos condiciones se agrega
    // al vector triggeringMuons, definido previamente
    // std::cout << "\n\n ------------>>>>>>>>TriggerObjectStandAlone TriggerObjectStandAlone TriggerObjectStandAlone TriggerObjectStandAlone"<<"\n";
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
      
      int_obj++; 
       // std::cout << "---->>obj number = " <<int_obj << "\n";
      obj.unpackFilterLabels(iEvent, *triggerBits);
      obj.unpackPathNames(names);
      bool isTriggerMuon = false;

      // checa que al menos un elemento sea un muon (ID=83)
      // que pasa con los demas?
      // std::cout << "\tfilterIds size:   " << obj.filterIds().size()<< "\n";
      for (unsigned h = 0; h < obj.filterIds().size(); ++h)
      if(obj.filterIds()[h] == 83){ 
        isTriggerMuon = true; 
         // std::cout << "\t\tType IDs:   " << 83 <<"\n";  //83 = muon
        //break;
      } else {
         // std::cout << "\t\tXXXXXXXXXX   Not Muon:   " << obj.filterIds()[h] <<"\n";
         isTriggerMuon = false;
      }

      if(!isTriggerMuon) continue;
      // Ahora checa que dentro de los filterlabes en al menos uno
      // exista hltL3 y Park
      isTriggerMuon = false;
       // std::cout << "\tfilterLabels size:  " << obj.filterLabels().size()<< "\n";
      for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
        
        std::string filterName = obj.filterLabels()[h];
         // std::cout << "\t\tfilterlabes:  " << h << filterName << "\n";

        if(filterName.find("hltL3") != std::string::npos  && filterName.find("Park") != std::string::npos){
            isTriggerMuon = true;
             // std::cout << "\t\tVVVVVVVV  Filter:   " << filterName<<"\n"; 
        }
	      //else{ isTriggerMuon = false; }
      }
      //std::cout << "\n\n\n";

      if(!isTriggerMuon) continue;
      triggeringMuons.push_back(obj);


      if(debug){ 
         std::cout << "\t\t\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      	 //Print trigger object collection and type
	      std::cout << "\t\t\tCollection: " << obj.collection() << std::endl;
       }
    }//trigger objects

    if(debug){
      
      // std::cout << "\n Total n of triggering muons = " << triggeringMuons.size() << std::endl;
      for(auto ij : triggeringMuons){
	       std::cout << " \t>>> components (pt, eta, phi) = (" << ij.pt() << ", " << ij.eta() << ", " << ij.phi() << ")\n";
      }
    }
     // std::cout << "\n\nTriggerObjectStandAlone TriggerObjectStandAlone TriggerObjectStandAlone TriggerObjectStandAlone <<<<<<------------"<<"\n";
    






    // Los objetos de salida
    std::unique_ptr<pat::MuonCollection>      trgmuons_out   ( new pat::MuonCollection );
    std::unique_ptr<pat::MuonCollection>      muons_out      ( new pat::MuonCollection );
    std::unique_ptr<TransientTrackCollection> trans_muons_out( new TransientTrackCollection );

    //now check for reco muons matched to triggering muons
    edm::Handle<std::vector<pat::Muon>> muons;
    iEvent.getByToken(muonSrc_, muons);

    // Se crea un vector de 0 con tamaño de muones ??? 
    std::vector<int> muonIsTrigger(muons->size(), 0);
    // Se itera sobre los elementos de muonCollection
     // std::cout << "--------->> muonCollection\n";
     // std::cout <<"Muons size:  " << muons->size() << "\n\n";
    
    //Se itera sobre todos los elementos de muons
    for(const pat::Muon & muon : *muons){

      //this is for triggering muon not really need to be configurable
      unsigned int iMuo(&muon - &(muons->at(0)) );
       // std::cout << "\tiMuo " << iMuo<< "\n";


       // Se verifica primero que sea LooseMuon y SoftMuon
       //                
       //  C1      C2     and    !and
       //
       //   T       T      T      F
       //   T       F      F      T
       //   F       T      F      T
       //   F       F      F      T

       // isLooseMuon    --  Particle identified as a muon by the Particle-Flow 
       //                    event reconstruction, and that is also reconstructed 
       //                    either as a global-muon or as an arbitrated tracker-muon
       //

       // isSoftMuon(PV) --  This selection requires the candidate to be a Tracker Muon, 
       //                     with the additional requirement that a muon segment is 
       //                     matched in both x and y coordinates with the extrapolated tracker track, 
       //                     such that the pull for local x and y is less than 3
       //                     This selection is used in quarkonia and B-physics analyses in CMS


      if(!(muon.isLooseMuon() && muon.isSoftMuon(PV))){
        // if (!(muon.isLooseMuon())){
        //   std::cout << "\n\t\tIs Not Loose Muon";
        // }
        // if (!muon.isSoftMuon(PV)){
        //   std::cout << "\n\t\tIs Not Soft Muon";
        // }
        continue;
      }

      float dRMuonMatching = -1.;
      int recoMuonMatching_index = -1;
      int trgMuonMatching_index = -1;
      // float dR_H = 0;

      // Se itera sobre todos los Triggering Muons, del primer for
      // y se calcula el dr = √(eta1-eta2)^2+(phi1-phi2)^2
      // Se considera el Triggering Muon que tenga el menor dr respecto 
      // a los slimmedMuon, y que ademas sea menor al maxdR impuesto por el usuario
      for(unsigned int iTrg=0; iTrg<triggeringMuons.size(); ++iTrg){
        float dR = reco::deltaR(triggeringMuons[iTrg], muon);
        // std::cout << "\n\t\tdR = " << dR << "\n";
	      if((dR < dRMuonMatching || dRMuonMatching == -1) && dR < maxdR_){
          dRMuonMatching = dR;
          recoMuonMatching_index = iMuo;
          trgMuonMatching_index = iTrg;
          float eta = muon.eta() - triggeringMuons[iTrg].eta();
          float phi = muon.phi() - triggeringMuons[iTrg].phi();
          // dR_H = std::sqrt(eta*eta+phi*phi);

          // std::cout << "\n\t\t dR_H"<< iTrg <<" = " << dR_H
          //   << "\n\t\treco (pt, eta, phi) = " << muon.pt() << " " << muon.eta() << " " << muon.phi() << " " 
          //   << "\n\t\tHLT (pt, eta, phi)  = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " " << triggeringMuons[iTrg].phi()
          //   << std::endl;
	      }
      }

      // Si algun Slimmed Muon satisfiso la condicion en dR con algun TrackMuon
      // entra en esta condicion, se guarda el muon "reco" a la coleccion trgmuons_out
      // se le agrega un entero "trgMuonIndex" el entero correspondiente al Trigger Muon
      // de la coleccion obtenida en el primer for
      // Se modifica la entrada del vector de ceros correspondiente a la
      // entrada del muon que se agrega
      
      //save reco muon 
      if(recoMuonMatching_index != -1){
        pat::Muon recoTriggerMuonCand (muon);
        recoTriggerMuonCand.addUserInt("trgMuonIndex", trgMuonMatching_index);
        trgmuons_out->emplace_back(recoTriggerMuonCand);
        //keep track of original muon index for SelectedMuons collection
        muonIsTrigger[iMuo] = 1;
      }
    }
    // std::cout << "\n<<<--------- muonCollection\n\n\n";








 

    // now produce output for analysis (code simplified loop of trg inside)
    // trigger muon + all compatible in dz with any tag
    for(unsigned int muIdx=0; muIdx<muons->size(); ++muIdx) {
      const pat::Muon& mu = (*muons)[muIdx];
      //selection cuts
      if (mu.pt() < ptMin_) continue;
      if (fabs(mu.eta()) > absEtaMax_) continue;
      //following ID is needed for trigger muons not here
      // anyway it is off in the configuration
      if (softMuonsOnly_ && !mu.isSoftMuon(PV)) continue;

      // same PV as the tag muon, both tag and probe only dz selection
      bool SkipMuon=true;
      for (const pat::Muon & trgmu : *trgmuons_out) {
        if( fabs(mu.vz()-trgmu.vz()) > dzTrg_cleaning_ && dzTrg_cleaning_ >0 )
        continue;
        SkipMuon=false;
      } 
      // needs decission: what about events without trg muon? now we SKIP them
      if (SkipMuon)  continue;


      // build transient track
      const reco::TransientTrack muonTT((*(mu.bestTrack())), &(*bFieldHandle)); //sara: check, why not using inner track for muons? 
      if (!muonTT.isValid()) continue;

      // // Momento del muon
      // math::XYZTLorentzVector momento = mu.p4();
      // mu.addUserFloat("p_X", momento.x());
      // mu.addUserFloat("p_Y", momento.y());
      // mu.addUserFloat("p_Y", momento.z());
      // mu.addUserFloat("p_T", momento.t());

      muons_out->emplace_back(mu);
      muons_out->back().addUserInt("isTriggering", muonIsTrigger[muIdx]);

      trans_muons_out->emplace_back(muonTT);



      }



    iEvent.put(std::move(trgmuons_out),    "trgMuons");
    iEvent.put(std::move(muons_out),       "SelectedMuons");
    iEvent.put(std::move(trans_muons_out), "SelectedTransientMuons");
}



DEFINE_FWK_MODULE(MuonTriggerSelector);
