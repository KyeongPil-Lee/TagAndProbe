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
// #include <DataFormats/MuonReco/interface/Muon.h>
#include "DataFormats/PatCandidates/interface/Muon.h"

// #include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
// #include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
// #include "DataFormats/TrackReco/interface/TrackBase.h"
// #include "DataFormats/BeamSpot/interface/BeamSpot.h"

// #include "TrackingTools/IPTools/interface/IPTools.h"
// #include "TrackingTools/Records/interface/TransientTrackRecord.h"
// #include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"


//
// class declaration
//

class ChangeValueMapDoubleToFloat : public edm::EDProducer {
public:
  explicit ChangeValueMapDoubleToFloat(const edm::ParameterSet&);
  ~ChangeValueMapDoubleToFloat();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  const edm::InputTag muons_;
  const edm::InputTag IsoDouble_;     
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
ChangeValueMapDoubleToFloat::ChangeValueMapDoubleToFloat(const edm::ParameterSet& iConfig):
muons_(iConfig.getParameter<edm::InputTag>("muons")),
IsoDouble_(iConfig.getParameter<edm::InputTag>("IsoDouble"))
{
  produces<edm::ValueMap<float> >();
}


ChangeValueMapDoubleToFloat::~ChangeValueMapDoubleToFloat()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ChangeValueMapDoubleToFloat::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace reco;
  using namespace edm;
  using namespace pat;

  // read input
  Handle< pat::MuonCollection > muons;
  iEvent.getByLabel(muons_,  muons);

  Handle< edm::ValueMap<double> > IsoDouble;
  iEvent.getByLabel(IsoDouble_, IsoDouble);

  // prepare vector for output    
  std::vector<float> Muon_IsoFloat;

  pat::MuonCollection::const_iterator amuon, endmuon = muons->end();
  Int_t qq = 0;
  // loop on PROBES
  for (amuon = muons->begin(); amuon != endmuon; ++amuon, ++qq)
  {
    edm::Ref<pat::MuonCollection> muRef(muons, qq);
    Muon_IsoFloat.push_back( (*IsoDouble)[muRef] );
  }// end loop on probes

  // convert into ValueMap and store
  std::auto_ptr<ValueMap<float> > IsoFloat(new ValueMap<float>());

  ValueMap<float>::Filler filler(*IsoFloat);

  filler.insert(muons, Muon_IsoFloat.begin(), Muon_IsoFloat.end());

  filler.fill();

  iEvent.put(IsoFloat);
}

// ------------ method called once each job just before starting event loop  ------------
void 
ChangeValueMapDoubleToFloat::beginJob()
{


}

// ------------ method called once each job just after ending the event loop  ------------
void 
ChangeValueMapDoubleToFloat::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ChangeValueMapDoubleToFloat);
