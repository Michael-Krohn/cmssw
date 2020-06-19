#ifndef SimG4Core_SaveSimTrack_H
#define SimG4Core_SaveSimTrack_H

#include "SimG4Core/Watcher/interface/SimWatcher.h"
#include "SimG4Core/Notification/interface/Observer.h"
#include "SimG4Core/Notification/interface/BeginOfTrack.h"
#include "SimG4Core/Notification/interface/EndOfTrack.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "G4ThreeVector.hh"

#include <vector>

class SaveSimTrack : public SimWatcher,
                     public Observer<const BeginOfTrack *>, 
		     public Observer<const EndOfTrack *> {

public:
  SaveSimTrack(edm::ParameterSet const & p);
  ~SaveSimTrack() override;
  void update(const BeginOfTrack * trk) override;
  void update(const EndOfTrack * trk) override;

private:
  std::vector<int> pdgs_;
  int MotherId;
  G4ThreeVector VertexPos;
};

#endif


