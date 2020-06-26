#ifndef SIMAPPLICATION_APRIMEPHYSICS_H_
#define SIMAPPLICATION_APRIMEPHYSICS_H 1_

// Geant4
#include "G4VPhysicsConstructor.hh"
#include "G4Decay.hh"
#include "G4hMultipleScattering.hh"
#include "G4ProcessManager.hh"
#include "G4Electron.hh"

// Sim Core
#include "G4APrime.h"
#include "G4muDarkBremsstrahlung.h"

class APrimePhysics : public G4VPhysicsConstructor {

   public:

      /**
       * Class constructor.
       * @param name The name of the physics.
       */
      APrimePhysics(const G4String& name = "APrime");

      /**
       * Class destructor.
       */
      virtual ~APrimePhysics();

      /**
       * Construct particles.
       */
      void ConstructParticle();

      /**
       * Construct the process.
       */
      void ConstructProcess();

   private:

      /**
       * Definition of the APrime particle.
       */
      G4ParticleDefinition* aprimeDef_;
      //G4Decay decayProcess;
      //G4hMultipleScattering scatterProcess;
  };

#endif
