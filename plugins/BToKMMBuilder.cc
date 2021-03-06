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

#include "TVector3.h"
#include "TMatrixD.h"

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "DataFormats/Math/interface/deltaR.h"
 
class BToKMMBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit BToKMMBuilder(const edm::ParameterSet &cfg):
    k_selection_{cfg.getParameter<std::string>("kaonSelection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    //HCL
    dzCut_{cfg.getParameter<double>("dzCut")},
    trgMuonToken_{consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("trgMuon"))},
    vertexSrc_{consumes<reco::VertexCollection> ( cfg.getParameter<edm::InputTag>( "vertexCollection" ) ) },
    //HCL
    dileptons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dileptons") )},
    leptons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("leptonTransientTracks") )},
    kaons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("kaons") )},
    kaons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("kaonsTransientTracks") )},
    isotracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    isolostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
    isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} {
      produces<pat::CompositeCandidateCollection>("Bcollection");
      //HCL
      produces<nanoaod::FlatTable>("VertexTable");
      //HCL
    }

  ~BToKMMBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::CompositeCandidate> k_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit

  //HCL
  const double dzCut_;
  const edm::EDGetTokenT<pat::MuonCollection> trgMuonToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_; 
  //HCL

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dileptons_;
  const edm::EDGetTokenT<TransientTrackCollection> leptons_ttracks_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> kaons_;
  const edm::EDGetTokenT<TransientTrackCollection> kaons_ttracks_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;
  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; 

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

void BToKMMBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  evt.getByToken(dileptons_, dileptons);
  
  edm::Handle<TransientTrackCollection> leptons_ttracks;
  evt.getByToken(leptons_ttracks_, leptons_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> kaons;
  evt.getByToken(kaons_, kaons);
  
  edm::Handle<TransientTrackCollection> kaons_ttracks;
  evt.getByToken(kaons_ttracks_, kaons_ttracks);  

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  





  //HCL PRIMERO OBTENGAMOS LOS VERTICES QUE NOS INTERESAN
  edm::Handle<reco::VertexCollection> vertexHandle;
  evt.getByToken(vertexSrc_, vertexHandle);

  edm::Handle<pat::MuonCollection> trgMuons;
  evt.getByToken(trgMuonToken_, trgMuons);

  std::vector<float> vx,vy,vz, dzTrgMu;// ,dlenSig,pAngle;

  const reco::Vertex & vertex = vertexHandle->front(); 
  const pat::Muon & trgmu = trgMuons->front();
  dzTrgMu.push_back(fabs(vertex.position().z()-trgmu.vz()));
  vx.push_back(vertex.position().x());
  vy.push_back(vertex.position().y());
  vz.push_back(vertex.position().z());  
  
  int nv = 0;
  for (const reco::Vertex & vertex : *vertexHandle){
    for (const pat::Muon & trgmu : *trgMuons){
      if(nv==3) break;
      if(fabs(vertex.position().z()-trgmu.vz())<dzCut_) {
        dzTrgMu.push_back(fabs(vertex.position().z()-trgmu.vz()));
        vx.push_back(vertex.position().x());
        vy.push_back(vertex.position().y());
        vz.push_back(vertex.position().z());
        nv++;
      }
    }
  }
  // if (dzTrgMu.size()==0){
  //      const reco::Vertex & vertex = vertexHandle->front(); 
  //      const pat::Muon & trgmu = trgMuons->front();
  //       dzTrgMu.push_back(fabs(vertex.position().z()-trgmu.vz()));
  //       vx.push_back(vertex.position().x());
  //       vy.push_back(vertex.position().y());
  //       vz.push_back(vertex.position().z());	
  // }
  // output
  auto pvTable = std::make_unique<nanoaod::FlatTable>(dzTrgMu.size(),"PV",false);

  pvTable->addColumn<float>("dzTrgMu",dzTrgMu,"abs(vz-trgMuon.vz)",nanoaod::FlatTable::FloatColumn,10);
  pvTable->addColumn<float>("vx",vx,"vx in cm",nanoaod::FlatTable::FloatColumn,10);
  pvTable->addColumn<float>("vy",vy,"vy in cm",nanoaod::FlatTable::FloatColumn,10);
  pvTable->addColumn<float>("vz",vz,"vz in cm",nanoaod::FlatTable::FloatColumn,10);
  //std::cout << "Pasamos el almacenamiento delos PVs\n\n" ;
  //HCL







  //for isolation
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);
  edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  evt.getByToken(isolostTracksToken_, iso_lostTracks);
  unsigned int nTracks     = iso_tracks->size();
  unsigned int totalTracks = nTracks + iso_lostTracks->size();

  std::vector<int> used_lep1_id, used_lep2_id, used_trk_id;

  //std::cout << "Kaons Size:  "<< kaons->size() << std::endl;
  //std::cout << "Dileptons Size:  "<< dileptons->size() << std::endl;
  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());

  for(size_t k_idx = 0; k_idx < kaons->size(); ++k_idx) {
    edm::Ptr<pat::CompositeCandidate> k_ptr(kaons, k_idx);
    
   // std::cout << "K DCASig : " <<k_ptr->userFloat("DCASig") << std::endl;
   // std::cout << "K nValidHits : " <<k_ptr->userInt("nValidHits") << std::endl;
    
   // try {
    if( !k_selection_(*k_ptr) ) continue;
    
    math::PtEtaPhiMLorentzVector k_p4(
      k_ptr->pt(), 
      k_ptr->eta(),
      k_ptr->phi(),
      K_MASS
      );


    for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
      
     // std::cout << "Dimu id : " << ll_idx << std::endl;

      edm::Ptr<pat::CompositeCandidate> ll_prt(dileptons, ll_idx);
      edm::Ptr<reco::Candidate> l1_ptr = ll_prt->userCand("l1");
      edm::Ptr<reco::Candidate> l2_ptr = ll_prt->userCand("l2");
      int l1_idx = ll_prt->userInt("l1_idx");
      int l2_idx = ll_prt->userInt("l2_idx");
    
      pat::CompositeCandidate cand;
      cand.setP4(ll_prt->p4() + k_p4);
      cand.setCharge(ll_prt->charge() + k_ptr->charge());
      // Use UserCands as they should not use memory but keep the Ptr itself
      // Put the lepton passing the corresponding selection
      cand.addUserCand("l1", l1_ptr);
      cand.addUserCand("l2", l2_ptr);
      cand.addUserCand("K", k_ptr);
      cand.addUserCand("dilepton", ll_prt);

      cand.addUserInt("l1_idx", l1_idx);
      cand.addUserInt("l2_idx", l2_idx);
      cand.addUserInt("k_idx", k_idx);
     
      auto dr_info = min_max_dr({l1_ptr, l2_ptr, k_ptr});
      cand.addUserFloat("min_dr", dr_info.first);
      cand.addUserFloat("max_dr", dr_info.second);

      float dr1 = reco::deltaR(*l1_ptr, *k_ptr);
      float dr2 = reco::deltaR(*l2_ptr, *k_ptr);
      cand.addUserFloat("dr_l1", dr1);
      cand.addUserFloat("dr_l2", dr2);
      cand.addUserFloat("dz_l1", fabs( l1_ptr->vz()-k_ptr->vz() ) );
      cand.addUserFloat("dz_l2", fabs( l2_ptr->vz()-k_ptr->vz() ) );

      // TODO add meaningful variables
      // Variables Pre fitting
      // muon1
      // cand.addUserFloat("l1_vx", l1_ptr->vx());
      // cand.addUserFloat("l1_vy", l1_ptr->vy());
      // cand.addUserFloat("l1_vz", l1_ptr->vz());
      // // muon2
      // cand.addUserFloat("l2_vx", l2_ptr->vx());
      // cand.addUserFloat("l2_vy", l2_ptr->vy());
      // cand.addUserFloat("l2_vz", l2_ptr->vz());
      // //kaon 
      // cand.addUserFloat("k_vx", k_ptr->vx());
      // cand.addUserFloat("k_vy", k_ptr->vy());
      // cand.addUserFloat("k_vz", k_ptr->vz());
      // //dimuon
      // cand.addUserFloat("dmuon_vx", ll_prt->vx());
      // cand.addUserFloat("dmuon_vy", ll_prt->vy());
      // cand.addUserFloat("dmuon_vz", ll_prt->vz());
      

      if( !pre_vtx_selection_(cand) ) continue;
    
      //HCL
      //std::cout << "BBB PreVertex Fitting    ";
    
      KinVtxFitter fitter(
        {leptons_ttracks->at(l1_idx), leptons_ttracks->at(l2_idx), kaons_ttracks->at(k_idx)},
        {l1_ptr->mass(), l2_ptr->mass(), K_MASS},
        {LEP_SIGMA, LEP_SIGMA, K_SIGMA} //some small sigma for the lepton mass
        );
      //HCL Se puede extraer el Vertex???
      //std::cout << "BBB PostVertex Fitting    \n";
      // RefCountedKinematicVertex myVertex = fitter.fitted_vtx_();
      // std::cout << myVertex->vertexIsValid << "\n\n";


      if(!fitter.success()) continue; // hardcoded, but do we need otherwise?
      cand.setVertex( 
        reco::Candidate::Point( 
          fitter.fitted_vtx().x(),
          fitter.fitted_vtx().y(),
          fitter.fitted_vtx().z()
          )  
        );
      
     // std::cout << "----->  SUCCESS!!!\n";
      used_lep1_id.emplace_back(l1_idx);
      used_lep2_id.emplace_back(l2_idx);
      used_trk_id.emplace_back(k_idx);
      cand.addUserInt("sv_OK" , fitter.success()); 
      cand.addUserFloat("sv_chi2", fitter.chi2());
      cand.addUserFloat("sv_ndof", fitter.dof()); // float??
      cand.addUserFloat("sv_prob", fitter.prob());
 
      cand.addUserFloat("fitted_mll" , (fitter.daughter_p4(0) + fitter.daughter_p4(1)).mass());
      cand.addUserFloat("fitted_pt_ll" , (fitter.daughter_p4(0) + fitter.daughter_p4(1)).pt());
      cand.addUserFloat("fitted_eta_ll" , (fitter.daughter_p4(0) + fitter.daughter_p4(1)).eta());
      cand.addUserFloat("fitted_phi_ll" , (fitter.daughter_p4(0) + fitter.daughter_p4(1)).phi());
      
      auto fit_p4 = fitter.fitted_p4();
      cand.addUserFloat("fitted_pt"  , fit_p4.pt()); 
      // cand.addUserFloat("fitted_px"  , fitter.fitted_candidate().globalMomentum().x()); 
      // cand.addUserFloat("fitted_py"  , fitter.fitted_candidate().globalMomentum().y()); 
      // cand.addUserFloat("fitted_pz"  , fitter.fitted_candidate().globalMomentum().z()); 
      cand.addUserFloat("fitted_eta" , fit_p4.eta());
      cand.addUserFloat("fitted_phi" , fit_p4.phi());
      cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass());      
      cand.addUserFloat("fitted_massErr", sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)));      
      cand.addUserFloat(
        "cos_theta_2D", 
        cos_theta_2D(fitter, *beamspot, cand.p4())
        );
      cand.addUserFloat(
        "fitted_cos_theta_2D", 
        cos_theta_2D(fitter, *beamspot, fit_p4)
        );
      auto lxy = l_xy(fitter, *beamspot);
      cand.addUserFloat("l_xy", lxy.value());
      cand.addUserFloat("l_xy_unc", lxy.error());

      // //Almacenemos la informacion del beamspot 
      // reco::BeamSpot bs = *beamspot;
      // cand.addUserFloat("beamSpot_x", bs.x(cand.vz()));
      // cand.addUserFloat("beamSpot_y", bs.y(cand.vz()));

      cand.addUserFloat("vtx_x", cand.vx());
      cand.addUserFloat("vtx_y", cand.vy());
      cand.addUserFloat("vtx_z", cand.vz());
      cand.addUserFloat("vtx_ex" , sqrt(fitter.fitted_vtx_uncertainty().cxx()));
      cand.addUserFloat("vtx_ey" , sqrt(fitter.fitted_vtx_uncertainty().cyy()));
      cand.addUserFloat("vtx_ez" , sqrt(fitter.fitted_vtx_uncertainty().czz()));
      try{
        cand.addUserFloat("vtx_eyx", fitter.fitted_vtx_uncertainty().cyx());
        cand.addUserFloat("vtx_ezx", fitter.fitted_vtx_uncertainty().czx());
        cand.addUserFloat("vtx_ezy", fitter.fitted_vtx_uncertainty().czy());
      }
      catch(...){
        cand.addUserFloat("vtx_eyx", -1);
        cand.addUserFloat("vtx_ezx", -1);
        cand.addUserFloat("vtx_ezy", -1);
      }
      cand.addUserFloat("fitted_l1_pt" , fitter.daughter_p4(0).pt()); 
      cand.addUserFloat("fitted_l1_eta", fitter.daughter_p4(0).eta());
      cand.addUserFloat("fitted_l1_phi", fitter.daughter_p4(0).phi());
      cand.addUserFloat("l1_charge", l1_ptr->charge());

      cand.addUserFloat("fitted_l2_pt" , fitter.daughter_p4(1).pt()); 
      cand.addUserFloat("fitted_l2_eta", fitter.daughter_p4(1).eta());
      cand.addUserFloat("fitted_l2_phi", fitter.daughter_p4(1).phi());
      cand.addUserFloat("l2_charge", l2_ptr->charge());

      cand.addUserFloat("fitted_k_pt"  , fitter.daughter_p4(2).pt()); 
      cand.addUserFloat("fitted_k_eta" , fitter.daughter_p4(2).eta());
      cand.addUserFloat("fitted_k_phi" , fitter.daughter_p4(2).phi());
      cand.addUserFloat("k_charge", k_ptr->charge());
      

      if( !post_vtx_selection_(cand) ) continue;        

      //compute isolation
      float l1_iso03 = 0;
      float l1_iso04 = 0;
      float l2_iso03 = 0;
      float l2_iso04 = 0;
      float k_iso03  = 0;
      float k_iso04  = 0;
      float b_iso03  = 0;
      float b_iso04  = 0;

      for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
      
        const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];
        // define selections for iso tracks (pT, eta, ...)
        if( !isotrk_selection_(trk) ) continue;
        // check if the track is the kaon
        if (k_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        // check if the track is one of the two leptons
        if (track_to_lepton_match(l1_ptr, iso_tracks.id(), iTrk) || 
            track_to_lepton_match(l2_ptr, iso_tracks.id(), iTrk) ) continue;

        // add to final particle iso if dR < cone
        float dr_to_l1 = deltaR(cand.userFloat("fitted_l1_eta"), cand.userFloat("fitted_l1_phi"), trk.eta(), trk.phi());
        float dr_to_l2 = deltaR(cand.userFloat("fitted_l2_eta"), cand.userFloat("fitted_l2_phi"), trk.eta(), trk.phi());
        float dr_to_k  = deltaR(cand.userFloat("fitted_k_eta") , cand.userFloat("fitted_k_phi") , trk.eta(), trk.phi());
        float dr_to_b  = deltaR(cand.userFloat("fitted_eta")   , cand.userFloat("fitted_phi") , trk.eta(), trk.phi());

        if (dr_to_l1 < 0.4){
          l1_iso04 += trk.pt();
          if ( dr_to_l1 < 0.3) l1_iso03 += trk.pt();
        }
        if (dr_to_l2 < 0.4){
          l2_iso04 += trk.pt();
          if (dr_to_l2 < 0.3)  l2_iso03 += trk.pt();
        }
        if (dr_to_k < 0.4){
          k_iso04 += trk.pt();
          if (dr_to_k < 0.3) k_iso03 += trk.pt();
        }
        if (dr_to_b < 0.4){
          b_iso04 += trk.pt();
          if (dr_to_b < 0.3) b_iso03 += trk.pt();
        }
      } 


      cand.addUserFloat("l1_iso03", l1_iso03);
      cand.addUserFloat("l1_iso04", l1_iso04);
      cand.addUserFloat("l2_iso03", l2_iso03);
      cand.addUserFloat("l2_iso04", l2_iso04);
      cand.addUserFloat("k_iso03" , k_iso03 );
      cand.addUserFloat("k_iso04" , k_iso04 );
      cand.addUserFloat("b_iso03" , b_iso03 );
      cand.addUserFloat("b_iso04" , b_iso04 );






      // Aqui creemos el boost al CM del dilepton
      math::XYZTLorentzVector dilep = ll_prt->p4();
      ROOT::Math::Boost cmboost(dilep.BoostToCM());

      math::XYZTLorentzVector kaonCM(  cmboost( k_ptr->p4() )  );
      math::XYZTLorentzVector muonCM1, muonCM2;

      //where the thetal is the angle between the
      //direction of the m-(m+) lepton and the K+(K-)
      if (l1_ptr->charge()==k_ptr->charge()){
        muonCM1 = cmboost(fitter.daughter_p4(1)) ;
        muonCM2 = cmboost(fitter.daughter_p4(0)) ;
      }
      else {
        muonCM1 = cmboost(fitter.daughter_p4(0)) ;
        muonCM2 = cmboost(fitter.daughter_p4(1)) ;
              } 

      float costhetaL = ( muonCM1.x()*muonCM2.x() 
                         + muonCM1.y()*muonCM2.y() 
                         + muonCM1.z()*muonCM2.z() ) / (muonCM1.P()*muonCM2.P() );

      float costhetaKL = ( muonCM1.x()*kaonCM.x()
                         + muonCM1.y()*kaonCM.y()
                         + muonCM1.z()*kaonCM.z() ) / (muonCM1.P()*kaonCM.P() );

      cand.addUserFloat("cosTheta_mm", costhetaL);
      cand.addUserFloat("cosTheta_km", costhetaKL);


      std::vector<float> cosAlpha(4, -2);
      std::vector<float> lxy_pv(4,-1);
      std::vector<float> errP(4,-1);

      // CHEQUEMOS QUE EL CANDIDATO A "B" SATIFAGA LOS CORTES EN COS(ALPHA) Y SIGNIFICANCIA
      for( unsigned int iPV=0; iPV<dzTrgMu.size(); ++iPV ){ 
        //Debemos calcular el vector que une los vertices primario y el ajustado
        math::XYZVector vPS(cand.vx()-vx.at(iPV), cand.vy()-vy.at(iPV), cand.vz()-vz.at(iPV));
        //Momento espacial del candidato
        math::XYZVector Bp(fitter.fitted_candidate().globalMomentum().x(), fitter.fitted_candidate().globalMomentum().y(), fitter.fitted_candidate().globalMomentum().z());
        //CosALPHA
        cosAlpha[iPV] = vPS.Dot(Bp)/(vPS.R()*Bp.R());
        
        //Para significancia:  
        GlobalError err = fitter.fitted_vtx_uncertainty();
        GlobalPoint delta(cand.vx()-vx.at(iPV), cand.vy()-vy.at(iPV), 0.);  

        lxy_pv[iPV] = delta.perp();
        errP[iPV] = sqrt(err.rerr(delta));

      }

      cand.addUserFloat("cosAlpha0", cosAlpha[0]);
      cand.addUserFloat("cosAlpha1", cosAlpha[1]);
      cand.addUserFloat("cosAlpha2", cosAlpha[2]);
      cand.addUserFloat("cosAlpha3", cosAlpha[3]);

      cand.addUserFloat("lxy_pv0", lxy_pv[0]);
      cand.addUserFloat("lxy_pv1", lxy_pv[1]);
      cand.addUserFloat("lxy_pv2", lxy_pv[2]);
      cand.addUserFloat("lxy_pv3", lxy_pv[3]);

      cand.addUserFloat("significance0", lxy_pv[0]/errP[0]);
      cand.addUserFloat("significance1", lxy_pv[1]/errP[1]);
      cand.addUserFloat("significance2", lxy_pv[2]/errP[2]);
      cand.addUserFloat("significance3", lxy_pv[3]/errP[3]);
      

      TVector3 pv, sv, pT;
      pT.SetXYZ(fit_p4.px(),fit_p4.py(),0.0);
      pv.SetXYZ(vx.at(0),vy.at(0),vz.at(0));
      sv.SetXYZ(cand.vx(),cand.vy(),cand.vz());


      TMatrix ESV(3,3);
      TMatrix EPV(3,3);

      ESV(0,0) = fitter.fitted_vtx_uncertainty().cxx();
      ESV(1,1) = fitter.fitted_vtx_uncertainty().cyy();
      ESV(2,2) = fitter.fitted_vtx_uncertainty().czz();
      ESV(0,1) = fitter.fitted_vtx_uncertainty().cyx();
      ESV(0,2) = fitter.fitted_vtx_uncertainty().czx();
      ESV(1,2) = fitter.fitted_vtx_uncertainty().czy();
       
      EPV(0,0) = vertexHandle->front().covariance(0,0);
      EPV(1,1) = vertexHandle->front().covariance(1,1);
      EPV(2,2) = vertexHandle->front().covariance(2,2);
      EPV(0,1) = vertexHandle->front().covariance(0,1);
      EPV(0,2) = vertexHandle->front().covariance(0,2);
      EPV(1,2) = vertexHandle->front().covariance(1,2);



      double ct, ect;
      V0_Lifetime(pv,sv,EPV,ESV, 5.27932, pT, ct, ect);

      
      cand.addUserFloat("PDL", ct);
      cand.addUserFloat("ePDL", ect);



      // VARIABLES DE TABLA DE KAONES 
      cand.addUserFloat("k_DCASig", k_ptr->userFloat("DCASig") );
      cand.addUserFloat("k_nValidHits", k_ptr->userInt("nValidHits") );
      cand.addUserInt("k_isMatchedToMuon", k_ptr->userInt("isMatchedToMuon") );
      cand.addUserInt("k_isMatchedToLooseMuon", k_ptr->userInt("isMatchedToLooseMuon") );
      cand.addUserInt("k_isMatchedToSoftMuon", k_ptr->userInt("isMatchedToSoftMuon") );
      cand.addUserInt("k_isMatchedToMediumMuon", k_ptr->userInt("isMatchedToMediumMuon") );
      cand.addUserInt("k_isMatchedToEle", k_ptr->userInt("isMatchedToEle") ); 

      cand.addUserInt("k_HighPurity", k_ptr->userInt("trackHighPurity") ); 
      cand.addUserInt("k_numberOfHits", k_ptr->userInt("numberOfHits") ); 
      cand.addUserInt("k_numberOfPixelHits", k_ptr->userInt("numberOfPixelHits") ); 
      cand.addUserInt("k_lostInnerHits", k_ptr->userInt("lostInnerHits") ); 

      // std::cout<<"\nBToKMMBuilder\n-------------\ntrackHighPurity  " << k_ptr->userInt("trackHighPurity") << "\n";
      // std::cout<<"numberOfHits  " << k_ptr->userInt("numberOfHits") << "\n";
      // std::cout<<"numberOfPixelHits  " << k_ptr->userInt("numberOfPixelHits") << "\n";
      // std::cout<<"lostInnerHits  " << k_ptr->userInt("lostInnerHits") << "\n";
      // std::cout << "------------\n\n";

 
      // VARIABLES DE TABLA DE MUONES
      // cand.addUserInt("l1_isPFMuon", l1_ptr->userInt("isPFMuon"));
      // cand.addUserInt("l1_isGlobalMuon", l1_ptr->userInt("isGlobalMuon"));
      // cand.addUserInt("l1_isTrackerMuon", l1_ptr->userInt("isTrackerMuon"));
      // cand.addUserInt("l1_isTriggering", l1_ptr->userInt("isTriggering"));
      // cand.addUserInt("l1_isSoft", l1_ptr->userInt("isSoft"));

      // cand.addUserInt("l2_isPFMuon", l2_ptr->userInt("isPFMuon"));
      // cand.addUserInt("l2_isGlobalMuon", l2_ptr->userInt("isGlobalMuon"));
      // cand.addUserInt("l2_isTrackerMuon", l2_ptr->userInt("isTrackerMuon"));
      // cand.addUserInt("l2_isTriggering", l2_ptr->userInt("isTriggering"));
      // cand.addUserInt("l2_isSoft", l2_ptr->userInt("isSoft"));



      ret_val->push_back(cand);




    } // for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
  
 //   }
   // catch(...){
     // std::cout << "Error!!!\n  "; 
    //}  

  } // for(size_t k_idx = 0; k_idx < kaons->size(); ++k_idx)



  for (auto & cand: *ret_val){
    cand.addUserInt("n_k_used", std::count(used_trk_id.begin(),used_trk_id.end(),cand.userInt("k_idx")));
    cand.addUserInt("n_l1_used", std::count(used_lep1_id.begin(),used_lep1_id.end(),cand.userInt("l1_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l1_idx")));
    cand.addUserInt("n_l2_used", std::count(used_lep1_id.begin(),used_lep1_id.end(),cand.userInt("l2_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l2_idx")));
  }


  evt.put(std::move(pvTable), "VertexTable");
  evt.put(std::move(ret_val), "Bcollection");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BToKMMBuilder);
