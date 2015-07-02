// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include <DataFormats/MuonReco/interface/Muon.h>

#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"


//
// class declaration
//

class MuonDxyDzPV : public edm::EDProducer {
public:
  explicit MuonDxyDzPV(const edm::ParameterSet&);
  ~MuonDxyDzPV();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  const edm::InputTag probes_;    
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
MuonDxyDzPV::MuonDxyDzPV(const edm::ParameterSet& iConfig):
probes_(iConfig.getParameter<edm::InputTag>("probes"))

{
  produces<edm::ValueMap<float> >("dxyBS");
  produces<edm::ValueMap<float> >("dzBS");
  produces<edm::ValueMap<float> >("dxyPV");
  produces<edm::ValueMap<float> >("dzPV");
}


MuonDxyDzPV::~MuonDxyDzPV()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuonDxyDzPV::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // read input
  Handle<View<reco::Muon> > probes;
  iEvent.getByLabel(probes_,  probes);
  
  edm::Handle<std::vector<reco::Vertex> > primaryVerticesHandle;
  iEvent.getByLabel("offlinePrimaryVertices", primaryVerticesHandle);
      
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot" ,recoBeamSpotHandle);
  reco::BeamSpot bs = *recoBeamSpotHandle;
  math::XYZPoint BSPosition;
  BSPosition = bs.position();

  // prepare vector for output    
  std::vector<double> muon_dxyBS;
  std::vector<double> muon_dzBS;
  std::vector<double> muon_dxyPV;
  std::vector<double> muon_dzPV;

  // fill
  View<reco::Muon>::const_iterator probe, endprobes = probes->end();

  // loop on PROBES
  for (probe = probes->begin(); probe != endprobes; ++probe) {
    
    Double_t dxyBS = -99999.;
    Double_t dzBS = -99999.;
    Double_t dxyPV = -99999.;
    Double_t dzPV = -99999.;

    if(probe->muonBestTrack().isNonnull())
    {  
      edm::ESHandle<TransientTrackBuilder> ttBuilder;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);
      
      reco::TrackRef muonTrack = probe->muonBestTrack();
      reco::TransientTrack tTrack = ttBuilder->build(*muonTrack);

      dxyBS = muonTrack->dxy(BSPosition);
      dzBS = muonTrack->dz(BSPosition);
      dxyPV = muonTrack->dxy(primaryVerticesHandle->at(0).position());
      dzPV = muonTrack->dz(primaryVerticesHandle->at(0).position());

      muon_dxyBS.push_back(dxyBS);
      muon_dzBS.push_back(dzBS);
      muon_dxyPV.push_back(dxyPV);
      muon_dzPV.push_back(dzPV);
    }

  }// end loop on probes

  // convert into ValueMap and store
  std::auto_ptr<ValueMap<float> > dxyBS(new ValueMap<float>());
  std::auto_ptr<ValueMap<float> > dzBS(new ValueMap<float>());
  std::auto_ptr<ValueMap<float> > dxyPV(new ValueMap<float>());
  std::auto_ptr<ValueMap<float> > dzPV(new ValueMap<float>());

  ValueMap<float>::Filler filler0(*dxyBS);
  ValueMap<float>::Filler filler1(*dzBS);
  ValueMap<float>::Filler filler2(*dxyPV);
  ValueMap<float>::Filler filler3(*dzPV);

  filler0.insert(probes, muon_dxyBS.begin(), muon_dxyBS.end());
  filler1.insert(probes, muon_dzBS.begin(), muon_dzBS.end());
  filler2.insert(probes, muon_dxyPV.begin(), muon_dxyPV.end());
  filler3.insert(probes, muon_dzPV.begin(), muon_dzPV.end());
  
  filler0.fill();
  filler1.fill();
  filler2.fill();
  filler3.fill();

  iEvent.put(dxyBS, "dxyBS");
  iEvent.put(dzBS, "dzBS");
  iEvent.put(dxyPV, "dxyPV");
  iEvent.put(dzPV, "dzPV");

}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonDxyDzPV::beginJob()
{


}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonDxyDzPV::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonDxyDzPV);
