#include "SimG4Core/SaveSimTrackAction/interface/SaveSimTrack.h"

#include "SimG4Core/Notification/interface/BeginOfTrack.h"
#include "SimG4Core/Notification/interface/TrackInformation.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "G4Track.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"
#include "G4VParticleChange.hh"
#include <algorithm>

SaveSimTrack::SaveSimTrack(edm::ParameterSet const & p) {

  edm::ParameterSet ps = p.getParameter<edm::ParameterSet>("SaveSimTrack");
  pdgs_ = ps.getUntrackedParameter<std::vector<int>>("PDGCodes");
  MotherId = -1;

  edm::LogInfo("SaveSimTrack") << "SaveSimTrack:: Save Sim Track if PDG code "
			       << "is one from the list of "  << pdgs_.size()
			       << " items";
  for (unsigned int k=0; k<pdgs_.size(); ++k)
    edm::LogInfo("SaveSimTrack") << "[" << k << "] " << pdgs_[k];
}

SaveSimTrack::~SaveSimTrack() {}
 
void SaveSimTrack::update(const BeginOfTrack * trk) {

  G4Track* theTrack = (G4Track*)((*trk)());
  TrackInformation * trkInfo = (TrackInformation *)(theTrack->GetUserInformation());
  if (trkInfo) {
    int pdg = theTrack->GetDefinition()->GetPDGEncoding();
    int mId = theTrack->GetParentID();
    G4ThreeVector Vpos = theTrack->GetVertexPosition();
    const G4VProcess* TrPro = theTrack->GetCreatorProcess();
    if(TrPro!=NULL)
    {
      if ((std::find(pdgs_.begin(),pdgs_.end(),pdg) == pdgs_.end())&&((theTrack->GetCreatorProcess()->GetProcessName())=="biasWrapper(muDBrem)"))
      {
        trkInfo->storeTrack(true);
	if(!theTrack->IsGoodForTracking()){theTrack->SetGoodForTrackingFlag(true);}
	//G4VParticleChange* pc = new G4ParticleChange();
	//pc->Initialize(*theTrack);
	//pc->ProposeWeight(theTrack->GetWeight()*5000);
	//printf("Deflected track ID is %d.\n",theTrack->GetTrackID());
	//printf("Weight after setting is %f.\n",theTrack->GetWeight());
	//theTrack->GetStep()->GetPostStepPoint()->SetWeight(theTrack->GetWeight()*5000);
      }
    }
    if (std::find(pdgs_.begin(),pdgs_.end(),pdg) != pdgs_.end()) {
      trkInfo->storeTrack(true);
      MotherId = mId;
      VertexPos = Vpos;
      LogDebug("SaveSimTrack") << "Save SimTrack the Track " 
			       << theTrack->GetTrackID() << " Type " 
			       << theTrack->GetDefinition()->GetParticleName()
			       << " Momentum " << theTrack->GetMomentum()/MeV 
			       << " MeV/c";
    }
  }
}

void SaveSimTrack::update(const EndOfTrack * trk){

  G4Track* theTrack = (G4Track*)((*trk)());
  TrackInformation * trkInfo = (TrackInformation *)(theTrack->GetUserInformation());
  const G4VProcess* TrPro = theTrack->GetCreatorProcess();
  if (trkInfo&&TrPro!=NULL) 
  {
    int pId = theTrack->GetTrackID();
    int mId = theTrack->GetParentID();
    G4ThreeVector Vpos = theTrack->GetVertexPosition();
    int pdg = theTrack->GetDefinition()->GetPDGEncoding();
    if (std::find(pdgs_.begin(),pdgs_.end(),pdg)==pdgs_.end()&&(theTrack->GetCreatorProcess()->GetProcessName())=="biasWrapper(muDBrem)")
    {
       trkInfo->storeTrack(true);
       std::cout << "Saving sibling track!\n";
    }
    /*if (pId == MotherId)
    {
       trkInfo->storeTrack(true);
       std::cout << "Saving mother track!\n";
    }*/
  }
}

