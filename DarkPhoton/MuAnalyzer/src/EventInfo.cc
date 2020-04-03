#include "DarkPhoton/MuAnalyzer/interface/EventInfo.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Common/interface/TriggerNames.h"

EventInfo::EventInfo(){}

bool EventInfo::goodPrimaryVertex(const edm::Event& iEvent, edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVertices_Label){

  using namespace reco;

  edm::Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(primaryVertices_Label, vertices);
	 
  std::vector<const reco::Vertex*> PVertices;
  std::vector<reco::Vertex>::const_iterator firstGoodVertex = vertices->end();

  for (std::vector<reco::Vertex>::const_iterator it=vertices->begin(); it!=firstGoodVertex; ++it) {
    if (!it->isFake() && it->ndof()>4 && it->position().Rho()<2. && std::abs(it->position().Z())<24.) {
      if(firstGoodVertex == vertices->end()){
        firstGoodVertex = it;
	PVertices.push_back(&(*it));
      }
      break;
    }
  }

  if(firstGoodVertex == vertices->end()){
    std::cout<<"NO GOOD VERTEX" << std::endl;
    return false;
  }else{
    bestVtx = *(firstGoodVertex);
    return true;
  }
					      
}

bool EventInfo::passTriggers(const edm::Event& iEvent, edm::EDGetToken m_trigResultsToken, std::vector<std::string> m_muonPathsToPass){

  bool passTriggers;

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(m_trigResultsToken, triggerResults);
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerResults);
  for(size_t i = 0; i < trigNames.size(); ++i) {
     const std::string &name = trigNames.triggerName(i);
     for(auto& pathName : m_muonPathsToPass){
        if((name.find(pathName) != std::string::npos )){
          if(triggerResults->accept(i)){
            passTriggers = true;
          }
        }
     }
  }

  return passTriggers;

}
