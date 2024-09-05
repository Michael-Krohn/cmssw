#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "DarkPhoton/MuAnalyzer/interface/metFilterTree.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "TTree.h"

metFilterTree::metFilterTree() {}

void metFilterTree::book(TFileDirectory histFolder) {
  failingEvents = histFolder.make<TTree>("failingEventInfo", "Tree containing event number and failed path");
  failingEvents->Branch("run", &m_run);
  failingEvents->Branch("lumi", &m_lumi);
  failingEvents->Branch("event", &m_event);
  failingEvents->Branch("path", &m_path);
  failingEvents->Branch("foundTrig", &m_foundTrig);
}

void metFilterTree::FillTree(int run, int lumi, int event, std::string path, bool foundTrig) {
  m_lumi = lumi;
  m_run = run;
  m_event = event;
  m_path = path;
  m_foundTrig = foundTrig;
  failingEvents->Fill();
}
