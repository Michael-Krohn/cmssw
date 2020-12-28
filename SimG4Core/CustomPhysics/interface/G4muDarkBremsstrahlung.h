/**
 * @file G4muDarkBremsstrahlung.h
 * @brief Class providing the Dark Bremsstrahlung process class.
 * @author Michael Revering, University of Minnesota
 */

#ifndef G4muDarkBremsstrahlung_h
#define G4muDarkBremsstrahlung_h

// Geant
#include "G4VEnergyLossProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4Positron.hh"
#include "G4SystemOfUnits.hh"
#include "G4APrime.h"
#include "G4muDarkBremsstrahlungModel.h"
#include "G4UnitsTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"


class G4Material;

class G4muDarkBremsstrahlung : public G4VEnergyLossProcess
{

   public:

      G4muDarkBremsstrahlung(const G4String& name = "muDBrem");

      virtual ~G4muDarkBremsstrahlung();

      virtual G4bool IsApplicable(const G4ParticleDefinition& p);

      virtual void PrintInfo();

      void SetMethod(std::string method_in);
      
      G4bool IsEnabled();
      void SetEnable(bool active);

   protected:

      virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
                                               const G4ParticleDefinition*);
      G4bool isInitialised;
      G4bool isEnabled;

   private:

      G4muDarkBremsstrahlung & operator=(const G4muDarkBremsstrahlung &right);
      G4muDarkBremsstrahlung(const G4muDarkBremsstrahlung&);

};


#endif
