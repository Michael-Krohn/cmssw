#include <memory>
#include <string>
#include <iostream>
#include <sstream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DarkPhoton/MuFilter/interface/ParticleListDrawer.h"

using namespace std;
using namespace reco;
using namespace edm;

void ParticleListDrawer::drawEvent(const edm::Event& iEvent, edm::EDGetToken srcToken_) {
   Handle<reco::CandidateView> particles;
   iEvent.getByToken(srcToken_, particles);
    
   ostringstream out;
   char buf[256];
 
   snprintf(buf,
            sizeof(buf),
            " idx  |    ID -       Name |Stat|  Mo1  Mo2  Da1  Da2 |nMo nDa|    pt       eta     phi   |     px       "
            "  py         pz        m     |");
   out << buf;
   snprintf(buf, sizeof(buf), "        vx       vy        vz     |");
   out << buf;
   out << endl;
 
   int idx = -1;
   int iMo1 = -1;
   int iMo2 = -1;
   int iDa1 = -1;
   int iDa2 = -1;

   vector<const reco::Candidate*> cands;
   vector<const Candidate*>::const_iterator found = cands.begin();
   for (CandidateView::const_iterator p = particles->begin(); p != particles->end(); ++p) {
     cands.push_back(&*p);
   }
 
   for (CandidateView::const_iterator p = particles->begin(); p != particles->end(); p++){
 
     // Particle Name
     int id = p->pdgId();
     string particleName = to_string(id);
 
     // Particle Index
     idx = p - particles->begin();
 
     // Particles Mothers and Daighters
     iMo1 = -1;
     iMo2 = -1;
     iDa1 = -1;
     iDa2 = -1;
     int nMo = p->numberOfMothers();
     int nDa = p->numberOfDaughters();
 
     found = find(cands.begin(), cands.end(), p->mother(0));
     if (found != cands.end())
       iMo1 = found - cands.begin();
 
     found = find(cands.begin(), cands.end(), p->mother(nMo - 1));
     if (found != cands.end())
       iMo2 = found - cands.begin();
 
     found = find(cands.begin(), cands.end(), p->daughter(0));
     if (found != cands.end())
       iDa1 = found - cands.begin();
 
     found = find(cands.begin(), cands.end(), p->daughter(nDa - 1));
     if (found != cands.end())
       iDa2 = found - cands.begin();
 
     char buf[2400];
     snprintf(
         buf,
         sizeof(buf),
         " %4d | %5d - %10s | %2d | %4d %4d %4d %4d | %2d %2d | %7.3f %10.3f %6.3f | %10.3f %10.3f %10.3f %8.3f |",
         idx,
         p->pdgId(),
         particleName.c_str(),
         p->status(),
         iMo1,
         iMo2,
         iDa1,
         iDa2,
         nMo,
         nDa,
         p->pt(),
         p->eta(),
         p->phi(),
         p->px(),
         p->py(),
         p->pz(),
         p->mass());
     out << buf;
 
     snprintf(buf, sizeof(buf), " %10.3f %10.3f %10.3f |", p->vertex().x(), p->vertex().y(), p->vertex().z());
     out << buf;
 
     const reco::GenParticle* gp = dynamic_cast<const reco::GenParticle*>(&*p);
     if (!gp)
       throw cms::Exception("Unsupported", "Status flags can be printed only for reco::GenParticle objects\n");
     if (gp->isPromptFinalState())
       out << "  PromptFinalState";
     if (gp->isDirectPromptTauDecayProductFinalState())
       out << "  DirectPromptTauDecayProductFinalState";
     if (gp->isHardProcess())
       out << "  HardProcess";
     if (gp->fromHardProcessFinalState())
       out << "  HardProcessFinalState";
     if (gp->fromHardProcessBeforeFSR())
       out << "  HardProcessBeforeFSR";
     if (gp->statusFlags().isFirstCopy())
       out << "  FirstCopy";
     if (gp->isLastCopy())
       out << "  LastCopy";
     if (gp->isLastCopyBeforeFSR())
       out << "  LastCopyBeforeFSR";
 
     out << endl;
   }
   LogImportant("ParticleListDrawer") << out.str();
   
 }
