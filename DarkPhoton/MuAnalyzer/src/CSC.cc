#include "DarkPhoton/MuAnalyzer/interface/CSC.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DarkPhoton/MuAnalyzer/interface/Histograms.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

CSC::CSC()
{
minDR=10;
recMinDR=100;
CSCdEta=-10;
CSCdPhi=-10;
minTotalImpactParameter=10000;
}

void CSC::ExtrapolateTrackToCSC(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label, std::vector<const reco::Track*>::const_iterator& iTrack, GlobalVector one_momentum, std::vector<reco::TransientTrack> tracksToVertex, GlobalPoint VertexPosition){

  edm::Handle<CSCSegmentCollection> TheCSCSegments;
  iEvent.getByToken(CSCSegment_Label, TheCSCSegments);

  edm::ESHandle<CSCGeometry> TheCSCGeometry;
  iSetup.get<MuonGeometryRecord>().get(TheCSCGeometry);
  if( TheCSCSegments.isValid() ){
    for(CSCSegmentCollection::const_iterator iSegment = TheCSCSegments->begin(); iSegment != TheCSCSegments->end(); iSegment++){

       CSCDetId iDetId = (CSCDetId)(*iSegment).cscDetId();
       if((*iTrack)->eta() < 0 && iDetId.endcap() == 1) continue;
       if((*iTrack)->eta() > 0 && iDetId.endcap() == 2) continue;


       //Only using 1st layer of CSCs
       //if(iDetId.station() != 1) continue;

       DetId TheDetUnitId(iSegment->cscDetId());
       const GeomDetUnit *TheUnit = (*TheCSCGeometry).idToDetUnit(TheDetUnitId);

       double dPhi = fabs(one_momentum.phi() - TheUnit->toGlobal(iSegment->localPosition()).phi());
       //double dPhi = fabs(one_momentum.phi() - TheUnit->position().phi());
       if(dPhi > ROOT::Math::Pi()) dPhi -= 2*ROOT::Math::Pi();

       LocalPoint TheLocalPosition = iSegment->localPosition();
       const BoundPlane& TheSurface = TheUnit->surface();
       GlobalPoint TheGlobalPosition = TheSurface.toGlobal(TheLocalPosition);

       if(minDR > sqrt(( pow((one_momentum.eta() - TheGlobalPosition.eta()),2.0) + pow(dPhi, 2.0)))){
         minDR = sqrt(( pow((one_momentum.eta() - TheGlobalPosition.eta()),2.0) + pow(dPhi, 2.0)));
//       if(minDR > sqrt(( pow((one_momentum.eta() - TheUnit->position().eta()),2.0) + pow(dPhi, 2.0)))){
//         minDR = sqrt(( pow((one_momentum.eta() - TheUnit->position().eta()),2.0) + pow(dPhi, 2.0)));
	 TrackEta_dR = one_momentum.eta();
	 TrackPhi_dR = one_momentum.phi();
         CSCdEta = one_momentum.eta()-TheGlobalPosition.eta();
         CSCdPhi = one_momentum.phi()-TheGlobalPosition.phi();
	 TrackGlobalPoint = GlobalPoint(GlobalPoint::Polar(one_momentum.theta(),one_momentum.phi(),TheGlobalPosition.mag()));
	 TrackP_dR = sqrt(pow(one_momentum.x(), 2) + pow(one_momentum.y(), 2) );
       }

       TrajectoryStateClosestToPoint  traj = tracksToVertex[0].trajectoryStateClosestToPoint(TheGlobalPosition);


       if (minTotalImpactParameter > sqrt(traj.perigeeParameters().transverseImpactParameter()*traj.perigeeParameters().transverseImpactParameter() + traj.perigeeParameters().longitudinalImpactParameter()*traj.perigeeParameters().longitudinalImpactParameter())){
	 minTotalImpactParameter = sqrt(traj.perigeeParameters().transverseImpactParameter()*traj.perigeeParameters().transverseImpactParameter() + traj.perigeeParameters().longitudinalImpactParameter()*traj.perigeeParameters().longitudinalImpactParameter());
	  TrackEta = one_momentum.eta();
	  TrackPhi = one_momentum.phi();
	  TrackP  = sqrt(pow(one_momentum.x(), 2) + pow(one_momentum.y(), 2) + pow(one_momentum.z(), 2));
       }

    }
  }

}

void CSC::ExtrapolateTrackToCSC(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label, edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitToken, const reco::Track* iTrack, reco::TransientTrack tracksToVertex, GlobalPoint VertexPosition){
   ExtrapolateTrackToCSC(iEvent, iSetup, CSCSegment_Label, iTrack, tracksToVertex, VertexPosition);
  
  edm::Handle<RPCRecHitCollection> rpcRecHits;
  iEvent.getByToken(rpcRecHitToken, rpcRecHits);
  edm::ESHandle<RPCGeometry> TheRPCGeometry;
  iSetup.get<MuonGeometryRecord>().get(TheRPCGeometry);

  RPCRecHitCollection::const_iterator recIt;
  for (recIt = rpcRecHits->begin(); recIt != rpcRecHits->end(); recIt++) 
  {
       DetId TheDetUnitId(recIt->rpcId());
       const GeomDetUnit *TheUnit = (*TheRPCGeometry).idToDetUnit(TheDetUnitId);

       LocalPoint TheLocalPosition = recIt->localPosition();
       const BoundPlane& TheSurface = TheUnit->surface();
       GlobalPoint TheGlobalPosition = TheSurface.toGlobal(TheLocalPosition);
       TrajectoryStateClosestToPoint  traj = tracksToVertex.trajectoryStateClosestToPoint(TheGlobalPosition);
       double RHdEta = traj.position().eta()-TheGlobalPosition.eta();
       double RHdPhi = traj.position().phi()-TheGlobalPosition.phi();
       if(recMinDR > sqrt(( pow(RHdEta,2.0) + pow(RHdPhi, 2.0))))
       {
          recMinDR= sqrt(( pow(RHdEta,2.0) + pow(RHdPhi, 2.0)));
       }
  }
}
void CSC::ExtrapolateTrackToCSC(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label, edm::EDGetTokenT<CSCRecHit2DCollection> CSCRecHitToken, const reco::Track* iTrack, reco::TransientTrack tracksToVertex, GlobalPoint VertexPosition){
   ExtrapolateTrackToCSC(iEvent, iSetup, CSCSegment_Label, iTrack, tracksToVertex, VertexPosition);
  
  edm::Handle<CSCRecHit2DCollection> cscRecHits;
  iEvent.getByToken(CSCRecHitToken, cscRecHits);
  edm::ESHandle<CSCGeometry> TheCSCGeometry;
  iSetup.get<MuonGeometryRecord>().get(TheCSCGeometry);

  CSCRecHit2DCollection::const_iterator recIt;
  for (recIt = cscRecHits->begin(); recIt != cscRecHits->end(); recIt++) 
  {
       CSCDetId iDetId = (CSCDetId)(*recIt).cscDetId();
       if(iTrack->eta() < 0 && iDetId.endcap() == 1) continue;
       if(iTrack->eta() > 0 && iDetId.endcap() == 2) continue;
       DetId TheDetUnitId(recIt->cscDetId());
       const GeomDetUnit *TheUnit = (*TheCSCGeometry).idToDetUnit(TheDetUnitId);

       LocalPoint TheLocalPosition = recIt->localPosition();
       const BoundPlane& TheSurface = TheUnit->surface();
       GlobalPoint TheGlobalPosition = TheSurface.toGlobal(TheLocalPosition);
       TrajectoryStateClosestToPoint  traj = tracksToVertex.trajectoryStateClosestToPoint(TheGlobalPosition);
       double RHdEta = traj.position().eta()-TheGlobalPosition.eta();
       double RHdPhi = traj.position().phi()-TheGlobalPosition.phi();
       if(recMinDR > sqrt(( pow(RHdEta,2.0) + pow(RHdPhi, 2.0))))
       {
          recMinDR= sqrt(( pow(RHdEta,2.0) + pow(RHdPhi, 2.0)));
       }
  }

}

void CSC::ExtrapolateTrackToCSC(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label, const reco::Track* iTrack, reco::TransientTrack tracksToVertex, GlobalPoint VertexPosition){

  edm::Handle<CSCSegmentCollection> TheCSCSegments;
  iEvent.getByToken(CSCSegment_Label, TheCSCSegments);
  edm::ESHandle<CSCGeometry> TheCSCGeometry;
  iSetup.get<MuonGeometryRecord>().get(TheCSCGeometry);
  minDR = 10;
  minTotalImpactParameter = 10000;
  if( TheCSCSegments.isValid() ){
    for(CSCSegmentCollection::const_iterator iSegment = TheCSCSegments->begin(); iSegment != TheCSCSegments->end(); iSegment++){

       CSCDetId iDetId = (CSCDetId)(*iSegment).cscDetId();
       if(iTrack->eta() < 0 && iDetId.endcap() == 1) continue;
       if(iTrack->eta() > 0 && iDetId.endcap() == 2) continue;


       //Only using 1st layer of CSCs
       //if(iDetId.station() != 1) continue;

       DetId TheDetUnitId(iSegment->cscDetId());
       const GeomDetUnit *TheUnit = (*TheCSCGeometry).idToDetUnit(TheDetUnitId);

       double dPhi = fabs(iTrack->momentum().phi() - TheUnit->toGlobal(iSegment->localPosition()).phi());
       //double dPhi = fabs(one_momentum.phi() - TheUnit->position().phi());
       if(dPhi > ROOT::Math::Pi()) dPhi -= 2*ROOT::Math::Pi();

       LocalPoint TheLocalPosition = iSegment->localPosition();
       const BoundPlane& TheSurface = TheUnit->surface();
       GlobalPoint TheGlobalPosition = TheSurface.toGlobal(TheLocalPosition);

       if(minDR > sqrt(( pow((iTrack->momentum().eta() - TheGlobalPosition.eta()),2.0) + pow(dPhi, 2.0)))){
          minDR = sqrt(( pow((iTrack->momentum().eta() - TheGlobalPosition.eta()),2.0) + pow(dPhi, 2.0)));
//       if(minDR > sqrt(( pow((one_momentum.eta() - TheUnit->position().eta()),2.0) + pow(dPhi, 2.0)))){
//         minDR = sqrt(( pow((one_momentum.eta() - TheUnit->position().eta()),2.0) + pow(dPhi, 2.0)));
         TrajectoryStateClosestToPoint  traj = tracksToVertex.trajectoryStateClosestToPoint(TheGlobalPosition);
         CSCdEta = traj.position().eta()-TheGlobalPosition.eta();
         CSCdPhi = traj.position().phi()-TheGlobalPosition.phi();

	 TrackEta_dR = iTrack->momentum().eta();
	 TrackPhi_dR = iTrack->momentum().phi();
	 TrackGlobalPoint = GlobalPoint(GlobalPoint::Polar(iTrack->momentum().theta(),iTrack->momentum().phi(),TheGlobalPosition.mag()));
	 TrackP_dR = sqrt(pow(iTrack->momentum().x(), 2) + pow(iTrack->momentum().y(), 2) );
         CSCTraj = TheUnit->toGlobal(iSegment->localDirection());
         CSCChiSq = iSegment->chi2();     
       }

       TrajectoryStateClosestToPoint  traj = tracksToVertex.trajectoryStateClosestToPoint(TheGlobalPosition);
        

       if (minTotalImpactParameter > sqrt(traj.perigeeParameters().transverseImpactParameter()*traj.perigeeParameters().transverseImpactParameter() + traj.perigeeParameters().longitudinalImpactParameter()*traj.perigeeParameters().longitudinalImpactParameter())){
	 minTotalImpactParameter = sqrt(traj.perigeeParameters().transverseImpactParameter()*traj.perigeeParameters().transverseImpactParameter() + traj.perigeeParameters().longitudinalImpactParameter()*traj.perigeeParameters().longitudinalImpactParameter());
	  TrackEta = iTrack->momentum().eta();
	  TrackPhi = iTrack->momentum().phi();
	  TrackP  = sqrt(pow(iTrack->momentum().x(), 2) + pow(iTrack->momentum().y(), 2) + pow(iTrack->momentum().z(), 2));
       }

    }
  }
}

void CSC::ExtrapolateMuonToCSC(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<CSCSegmentCollection > CSCSegment_Label, const reco::Muon* iMuon, GlobalVector two_momentum, std::vector<reco::TransientTrack> tracksToVertex){

  edm::Handle<CSCSegmentCollection> TheCSCSegments;
  iEvent.getByToken(CSCSegment_Label, TheCSCSegments);

  edm::ESHandle<CSCGeometry> TheCSCGeometry;
  iSetup.get<MuonGeometryRecord>().get(TheCSCGeometry);

  if( TheCSCSegments.isValid() ){
    for(CSCSegmentCollection::const_iterator iSegment = TheCSCSegments->begin(); iSegment != TheCSCSegments->end(); iSegment++){
       CSCDetId iDetId = (CSCDetId)(*iSegment).cscDetId();

       //Only using 1st layer of CSCs
//       if(iDetId.station() != 1) continue;
       
       if(iMuon->eta() < 0 && iDetId.endcap() == 1) continue;
       if(iMuon->eta() > 0 && iDetId.endcap() == 2) continue;

       DetId TheDetUnitId(iSegment->cscDetId());
       const GeomDetUnit *TheUnit = (*TheCSCGeometry).idToDetUnit(TheDetUnitId);

       LocalPoint TheLocalPosition = iSegment->localPosition();
//       const BoundPlane& TheSurface = TheUnit->surface();
       GlobalPoint TheGlobalPosition = TheUnit->toGlobal(TheLocalPosition);

       if(fabs(iMuon->eta()) > 1.653){
         TrajectoryStateClosestToPoint  MuonTraj = tracksToVertex[1].trajectoryStateClosestToPoint(TheGlobalPosition);
         if (minTotalImpactParameter_Muon > sqrt(MuonTraj.perigeeParameters().transverseImpactParameter()*MuonTraj.perigeeParameters().transverseImpactParameter() + MuonTraj.perigeeParameters().longitudinalImpactParameter()*MuonTraj.perigeeParameters().longitudinalImpactParameter())){
           minTotalImpactParameter_Muon = sqrt(MuonTraj.perigeeParameters().transverseImpactParameter()*MuonTraj.perigeeParameters().transverseImpactParameter() + MuonTraj.perigeeParameters().longitudinalImpactParameter()*MuonTraj.perigeeParameters().longitudinalImpactParameter());
	   MuonEta = iMuon->eta();
           MuonPhi = iMuon->phi();
           MuonP = sqrt(pow(two_momentum.x(), 2) + pow(two_momentum.y(), 2) + pow(two_momentum.z(), 2));
	   
         }
         double dPhi_Muon = fabs(two_momentum.phi() - TheGlobalPosition.phi());
         if(dPhi_Muon > ROOT::Math::Pi()) dPhi_Muon -= 2*ROOT::Math::Pi();
         if (minDR_Muon > sqrt(( pow((two_momentum.eta() - TheGlobalPosition.eta()),2.0) + pow(dPhi_Muon, 2.0)))){
            minDR_Muon = sqrt(( pow((two_momentum.eta() - TheGlobalPosition.eta()),2.0) + pow(dPhi_Muon, 2.0)));
	    MuonEta_dR = iMuon->eta();
            MuonPhi_dR = iMuon->phi();
	    MuonGlobalPoint = GlobalPoint(GlobalPoint::Polar(two_momentum.theta(),two_momentum.phi(),TheGlobalPosition.mag()));
            MuonP_dR = sqrt(pow(two_momentum.x(), 2) + pow(two_momentum.y(), 2)+pow(two_momentum.z(),2));
	 }
       }
    }
  }

}


