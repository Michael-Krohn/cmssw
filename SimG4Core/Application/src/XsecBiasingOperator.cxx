/**
 * @file XsecBiasingPlugin.cxx
 * @brief Geant4 Biasing Operator used to bias the occurence of 
 *        events by modifying the cross-section.
 * @author Michael Revering, modified from file by Omar Moreno
 *         University of Minnesota
 */

#include "SimG4Core/Application/interface/XsecBiasingOperator.h"
#include "SimG4Core/CustomPhysics/interface/G4APrime.h"
#include "G4MuonMinus.hh"

XsecBiasingOperator::XsecBiasingOperator(std::string name) :
    G4VBiasingOperator(name) { 
}

XsecBiasingOperator::~XsecBiasingOperator() {
}

void XsecBiasingOperator::StartRun() { 

    if (particleType_.compare("gamma") == 0) {
        processManager_ = G4Gamma::GammaDefinition()->GetProcessManager();
    } else if (particleType_.compare("e-") == 0) { 
        processManager_ = G4Electron::ElectronDefinition()->GetProcessManager();
    } else if (particleType_.compare("mu-") == 0) {
        processManager_ = G4MuonMinus::MuonMinusDefinition()->GetProcessManager();
    }
    else { 
        // Throw an exception
    }

    std::cout << "[ XsecBiasingOperator ]: Biasing particles of type " 
              << particleType_ << std::endl; 

    if (processIsBiased(this->getProcessToBias())) {
        std::cout << "Successfully found the process.\n";
        xsecOperation = new G4BOptnChangeCrossSection("changeXsec-" + this->getProcessToBias());
    } else { 
        std::cout << "Failed to find the process.\n";
        // Throw an exception
    }
}

bool XsecBiasingOperator::processIsBiased(std::string process) { 
    
    // Loop over all processes and check if the given process is being 
    // biased.
    const G4BiasingProcessSharedData* sharedData 
        = G4BiasingProcessInterface::GetSharedData(processManager_);
    if (sharedData) {
        for (size_t iprocess = 0 ; 
                iprocess < (sharedData->GetPhysicsBiasingProcessInterfaces()).size(); ++iprocess) {
            
            const G4BiasingProcessInterface* wrapperProcess
                = (sharedData->GetPhysicsBiasingProcessInterfaces())[iprocess];
            if (wrapperProcess->GetWrappedProcess()->GetProcessName().compareTo(process) == 0) {
                return true; 
            } 
        }
    }
    return false;
}
