/**
 * @file DarkBremXsecBiasingPlugin.h
 *         Modified by Michael Revering from file by Omar Moreno
 */

#ifndef BIASING_DARKBREMXSECBIASINGOPERATOR_H_
#define BIASING_DARKBREMXSECBIASINGOPERATOR_H_

//----------------//
//   C++ StdLib   //
//----------------//
#include <algorithm>

//------------//
//   Geant4   //
//------------//
#include "G4BiasingProcessInterface.hh"
#include "G4BiasingProcessSharedData.hh"
#include "G4BOptnChangeCrossSection.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4VBiasingOperator.hh"

#include "XsecBiasingOperator.h"

//-----------//
//   Root    //
//-----------//

class DarkBremXsecBiasingOperator : public XsecBiasingOperator { 

    public: 

        /** Constructor */
        DarkBremXsecBiasingOperator(std::string name);

        /** Destructor */
        ~DarkBremXsecBiasingOperator();

        /** Method called at the beginning of a run. */
        void StartRun();

        /** Method called at the end of a run. */
	void EndRun();

        /** 
         * @return Method that returns the biasing operation that will be used
         *         to bias the occurence of events.
         */
        G4VBiasingOperation* ProposeOccurenceBiasingOperation(const G4Track* track,
                const G4BiasingProcessInterface* callingProcess);
        

    protected:

        virtual std::string getProcessToBias() { return DARKBREM_PROCESS; }

    private: 

        /** Geant4 photonuclear process name. */
        static const std::string DARKBREM_PROCESS;

        /** Cross-section biasing operation */
        G4BOptnChangeCrossSection* emXsecOperation{nullptr};

        /** Unbiased photonuclear xsec. */
        double dbXsecUnbiased_{0};

        /** Biased photonuclear xsec. */
        double dbXsecBiased_{0}; 

};  // DarkBremXsecBiasingOperator


#endif // SIMPLUGINS_DARKBREMXSECBIASINGOPERATOR_H_ 
