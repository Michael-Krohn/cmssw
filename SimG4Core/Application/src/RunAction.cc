#include "SimG4Core/Application/interface/RunAction.h"
#include "SimG4Core/Application/interface/SimRunInterface.h"

#include "SimG4Core/Notification/interface/BeginOfRun.h"
#include "SimG4Core/Notification/interface/EndOfRun.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimG4Core/Application/interface/XsecBiasingOperator.h"
#include "SimG4Core/Application/interface/DarkBremXsecBiasingOperator.h"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ProductionCuts.hh"
#include "G4Timer.hh" 
#include <iostream>
#include <fstream>

RunAction::RunAction(const edm::ParameterSet& p, SimRunInterface* rm, bool master) 
  : m_runInterface(rm), 
    m_stopFile(p.getParameter<std::string>("StopFile")),
    m_timer(nullptr), m_isMaster(master)
{
   xsecBiasing = nullptr;
}

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run * aRun)
{
  if (std::ifstream(m_stopFile.c_str()))
    {
      edm::LogWarning("SimG4CoreApplication")
        << "RunAction::BeginOfRunAction: termination signal received";
      m_runInterface->abortRun(true);
    }
 
  xsecBiasing = new DarkBremXsecBiasingOperator("DarkBremXsecBiasingOperator");
  for (G4LogicalVolume* volume : *G4LogicalVolumeStore::GetInstance())
  {
     G4String volumeName = volume->GetName();
     //if(volumeName.contains("MBBT")||volumeName.contains("MBAT")||volumeName.contains("CALO")||volumeName.contains("VCAL"))
     if(volumeName.contains("HE")&&(volumeName!="CHEL")&&(volumeName!="RHEX"))
     {
        xsecBiasing->AttachTo(volume);
        std::cout << "Attaching biasing operator " << xsecBiasing->GetName() << " to volume " << volume->GetName() << std::endl;
     }
  }

  xsecBiasing->StartRun(); 
  BeginOfRun r(aRun);
  m_beginOfRunSignal(&r);
  /*
  if (m_isMaster) {
    m_timer = new G4Timer();
    m_timer->Start();
  }
  */
}

void RunAction::EndOfRunAction(const G4Run * aRun)
{
  if (isMaster) {
    edm::LogInfo("SimG4CoreApplication") 
      << "RunAction: total number of events "  << aRun->GetNumberOfEvent();
    if(m_timer) {
      m_timer->Stop();
      edm::LogInfo("SimG4CoreApplication") 
	<< "RunAction: Master thread time  "  << *m_timer;
      // std::cout << "\n" << "Master thread time:  "  << *m_timer << std::endl;
      delete m_timer;
    }
  }
  if (std::ifstream(m_stopFile.c_str()))
    {
      edm::LogWarning("SimG4CoreApplication")
        << "RunAction::EndOfRunAction: termination signal received";
      m_runInterface->abortRun(true);
    }
  xsecBiasing->EndRun();  
  EndOfRun r(aRun);
  m_endOfRunSignal(&r);
}

