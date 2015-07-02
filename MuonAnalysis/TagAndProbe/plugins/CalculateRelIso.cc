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

class CalculateRelIso : public edm::EDProducer {
public:
  explicit CalculateRelIso(const edm::ParameterSet&);
  ~CalculateRelIso();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  const edm::InputTag muons_;
  const edm::InputTag SumCH_;
  const edm::InputTag SumNH_;
  const edm::InputTag SumPh_;     
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
CalculateRelIso::CalculateRelIso(const edm::ParameterSet& iConfig):
muons_(iConfig.getParameter<edm::InputTag>("muons")),
SumCH_(iConfig.getParameter<edm::InputTag>("SumCH")),
SumNH_(iConfig.getParameter<edm::InputTag>("SumNH")),
SumPh_(iConfig.getParameter<edm::InputTag>("SumPh"))
{
  produces<edm::ValueMap<float> >();
}


CalculateRelIso::~CalculateRelIso()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
CalculateRelIso::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace reco;
  using namespace edm;
  using namespace pat;

  // read input
  Handle< pat::MuonCollection > muons;
  iEvent.getByLabel(muons_,  muons);

  Handle< edm::ValueMap<double> > SumCH;
  iEvent.getByLabel(SumCH_, SumCH);
  Handle< edm::ValueMap<double> > SumNH;
  iEvent.getByLabel(SumNH_, SumNH);
  Handle< edm::ValueMap<double> > SumPh;
  iEvent.getByLabel(SumPh_, SumPh);

  // prepare vector for output    
  std::vector<float> Muon_RelIso;

  pat::MuonCollection::const_iterator amuon, endmuon = muons->end();
  Int_t qq = 0;
  // loop on PROBES
  for (amuon = muons->begin(); amuon != endmuon; ++amuon, ++qq)
  {
    edm::Ref<pat::MuonCollection> muRef(muons, qq);
    double SumE = (*SumCH)[muRef] + (*SumNH)[muRef] + (*SumPh)[muRef];
    double Pt = amuon->pt();

    Muon_RelIso.push_back( SumE / Pt );

  }// end loop on probes

  // convert into ValueMap and store
  std::auto_ptr<ValueMap<float> > RelIso(new ValueMap<float>());

  ValueMap<float>::Filler filler(*RelIso);

  filler.insert(muons, Muon_RelIso.begin(), Muon_RelIso.end());

  filler.fill();

  iEvent.put(RelIso);
}

// ------------ method called once each job just before starting event loop  ------------
void 
CalculateRelIso::beginJob()
{


}

// ------------ method called once each job just after ending the event loop  ------------
void 
CalculateRelIso::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(CalculateRelIso);
