#include "SimG4Core/CustomPhysics/interface/G4muDarkBremsstrahlung.h"
#include "G4LossTableManager.hh"

using namespace std;

G4muDarkBremsstrahlung::G4muDarkBremsstrahlung(const G4String& name):
   G4VEnergyLossProcess(name),
   isInitialised(false)
{  
   G4int subtype = 63;   
   SetProcessSubType(subtype);
   SetSecondaryParticle(G4APrime::APrime());
   SetIonisation(false);
}

G4muDarkBremsstrahlung::~G4muDarkBremsstrahlung()
{}

G4bool G4muDarkBremsstrahlung::IsApplicable(const G4ParticleDefinition& p)
{
   return (&p == G4MuonPlus::MuonPlus() || &p == G4MuonMinus::MuonMinus());
}

void G4muDarkBremsstrahlung::InitialiseEnergyLossProcess(const G4ParticleDefinition*,
                                                    const G4ParticleDefinition*)
{
   if(!isInitialised)
   {  
      //SetEmModel(new G4muDarkBremsstrahlungModel(),1);
      //if(!EmModel()) {SetEmModel(new G4muDarkBremsstrahlungModel());}
      G4VEmFluctuationModel* fm = 0;
      AddEmModel(0, new G4muDarkBremsstrahlungModel(), fm);

      //EmModel(0)->SetLowEnergyLimit(MinKinEnergy());
      //EmModel(0)->SetHighEnergyLimit(energyLimit);
        
      //AddEmModel(1, EmModel(1), fm);
      isInitialised = true;
      isEnabled = true;
   }

   //EmModel(0)->SetSecondaryThreshold(eth);
   //EmModel(0)->SetLPMFlag(false);
}

void G4muDarkBremsstrahlung::PrintInfo()
{
   if(EmModel(1))
   {
//      G4cout << "    LPM flag: " << "false " << " for E > " << EmModel(1)->HighEnergyLimit()/GeV<< " GeV";
//      G4cout << G4endl;
      
   }
}

void G4muDarkBremsstrahlung::SetMethod(std::string method_in)
{
   ((G4muDarkBremsstrahlungModel*)EmModel(1))->SetMethod(method_in);
   return;
}

G4bool G4muDarkBremsstrahlung::IsEnabled()
{
   return isEnabled;
}

void G4muDarkBremsstrahlung::SetEnable(bool state)
{
   isEnabled = state;
   return;
}
