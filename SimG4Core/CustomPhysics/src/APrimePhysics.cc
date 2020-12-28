#include "SimG4Core/CustomPhysics/interface/APrimePhysics.h"

// Geant4
#include "G4SystemOfUnits.hh"


    APrimePhysics::APrimePhysics(const G4String& name) :
            G4VPhysicsConstructor(name), aprimeDef_(nullptr) {
    }

    APrimePhysics::~APrimePhysics() {
    }

    void APrimePhysics::ConstructParticle() {

       /**
        * Insert A-prime into the Geant4 particle table.
        * For now we flag it as stable.
        */
    aprimeDef_ = G4APrime::APrime();

    //aprimeDef->SetProcessManager(new G4ProcessManager(aprimeDef));
    }

    void APrimePhysics::ConstructProcess() {
        /*
         G4ProcessManager* pm = aprimeDef->GetProcessManager();
         if (pm != NULL) {
         pm->AddProcess(&scatterProcess, -1, 1, 1);
         pm->AddProcess(&decayProcess, -1, -1, 2);
         } else {
         G4Exception("APrimePhysics::ConstructProcess",
         "InitializationError",
         FatalException,
         "The process manager for APrime is NULL.");
         }
         */
//   G4ParticleDefinition* particle = G4Electron::ElectronDefinition();
   G4ParticleDefinition* muonminus = G4MuonMinus::MuonMinusDefinition();
   G4ParticleDefinition* muonplus = G4MuonPlus::MuonPlusDefinition();
   G4ProcessManager* pmplus = muonplus->GetProcessManager();
   G4ProcessManager* pmminus = muonminus->GetProcessManager();
   pmplus->AddProcess(new G4muDarkBremsstrahlung(), -1, 1, 1);
   pmminus->AddProcess(new G4muDarkBremsstrahlung(), -1, 1, 1);
//   pm->AddProcess(new G4muDarkBremsstrahlung(), -1, 1, 1);
}
