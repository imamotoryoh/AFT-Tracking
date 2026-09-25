#include "AFTTrackEvaluator.hh"

#include <iostream>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

AFTTrackEvaluator::AFTTrackEvaluator(const AFTGeometry& geometry)
  : geometry_(geometry)
{
  nCandidatesXZ_ =
    new TH1D(
      "hNCandidatesXZ",
      "XZ candidates per event;N candidates;Events",
      21, -0.5, 20.5);

  nCandidatesYZ_ =
    new TH1D(
      "hNCandidatesYZ",
      "YZ candidates per event;N candidates;Events",
      21, -0.5, 20.5);

  nCandidatesTotal_ =
    new TH1D(
      "hNCandidatesTotal",
      "All candidates per event;N candidates;Events",
      41, -0.5, 40.5);

  nHitsXZ_ =
    new TH1D(
      "hCandidateNHitsXZ",
      "XZ candidates;N hits;Candidates",
      19, -0.5, 18.5);

  nHitsYZ_ =
    new TH1D(
      "hCandidateNHitsYZ",
      "YZ candidates;N hits;Candidates",
      19, -0.5, 18.5);

  startXZ_ =
    new TH1D(
      "hCandidateStartLayerXZ",
      "XZ candidates;start layer;Candidates",
      18, -0.5, 17.5);

  startYZ_ =
    new TH1D(
      "hCandidateStartLayerYZ",
      "YZ candidates;start layer;Candidates",
      18, -0.5, 17.5);

  endXZ_ =
    new TH1D(
      "hCandidateEndLayerXZ",
      "XZ candidates;end layer;Candidates",
      18, -0.5, 17.5);

  endYZ_ =
    new TH1D(
      "hCandidateEndLayerYZ",
      "YZ candidates;end layer;Candidates",
      18, -0.5, 17.5);

  missingXZ_ =
    new TH1D(
      "hCandidateMissingXZ",
      "XZ candidates;N missing layers;Candidates",
      10, -0.5, 9.5);

  missingYZ_ =
    new TH1D(
      "hCandidateMissingYZ",
      "YZ candidates;N missing layers;Candidates",
      10, -0.5, 9.5);

  residualXZ_ =
    new TH1D(
      "hCandidateResidualXZ",
      "XZ candidates;sum residual^{2} [mm^{2}];Candidates",
      200, 0.0, 200.0);

  residualYZ_ =
    new TH1D(
      "hCandidateResidualYZ",
      "YZ candidates;sum residual^{2} [mm^{2}];Candidates",
      200, 0.0, 200.0);

  hEndHitVsDeltaXZ_ =
    new TH2D(
	     "hEndHitVsDeltaXZ",
	     "XZ track displacement;Last hit layer;X_{last} - X_{first} [mm]",
	     18, -0.5, 17.5,
	     80, -80.0, 80.0
	     );
  
  hEndHitVsDeltaYZ_ =
    new TH2D(
	     "hEndHitVsDeltaYZ",
	     "YZ track displacement;Last hit layer;Y_{last} - Y_{first} [mm]",
	     18, -0.5, 17.5,
	     80, -80.0, 80.0
	     );

  hEndHitPositionXZ_ =
    new TH2D(
	     "hEndHitPositionXZ",
	     "XZ end hit position;Last hit layer;X_{last} [mm]",
	     18, -0.5, 17.5,
	     64, -55.0, 55.0
	     );

  hEndHitPositionYZ_ =
    new TH2D(
	     "hEndHitPositionYZ",
	     "YZ end hit position;Last hit layer;Y_{last} [mm]",
	     18, -0.5, 17.5,
	     64, -55.0, 55.0
	     );
  
  for(auto* hist :
      {nCandidatesXZ_,
       nCandidatesYZ_,
       nCandidatesTotal_,
       nHitsXZ_,
       nHitsYZ_,
       startXZ_,
       startYZ_,
       endXZ_,
       endYZ_,
       missingXZ_,
       missingYZ_,
       residualXZ_,
       residualYZ_}){

    hist->SetDirectory(nullptr);
  }
  for(auto* hist :
	{hEndHitVsDeltaXZ_,
	 hEndHitVsDeltaYZ_,
	 hEndHitPositionXZ_,
	 hEndHitPositionYZ_}){
    hist->SetDirectory(nullptr);
  }

}

AFTTrackEvaluator::~AFTTrackEvaluator()
{
  for(auto* hist :
      {nCandidatesXZ_,
       nCandidatesYZ_,
       nCandidatesTotal_,
       nHitsXZ_,
       nHitsYZ_,
       startXZ_,
       startYZ_,
       endXZ_,
       endYZ_,
       missingXZ_,
       missingYZ_,
       residualXZ_,
       residualYZ_}){

    delete hist;
  }
  for(auto* hist :
	{hEndHitVsDeltaXZ_,
	 hEndHitVsDeltaYZ_,
	 hEndHitPositionXZ_,
	 hEndHitPositionYZ_}){

    delete hist;
  }

}

void
AFTTrackEvaluator::FillCandidate(
    const AFTTrackResult& result)
{
  const bool isXZ =
    result.projection ==
    AFTProjection::XZ;

  TH1D* nHits =
    isXZ ? nHitsXZ_ : nHitsYZ_;

  TH1D* start =
    isXZ ? startXZ_ : startYZ_;

  TH1D* end =
    isXZ ? endXZ_ : endYZ_;

  TH1D* missing =
    isXZ ? missingXZ_ : missingYZ_;

  TH1D* residual =
    isXZ ? residualXZ_ : residualYZ_;

  nHits->Fill(
    static_cast<double>(
      result.track.GetNHits()));

  start->Fill(result.startLayer);
  end->Fill(result.endLayer);
  missing->Fill(result.nMissing);
  residual->Fill(result.geometryResidual);

  const auto& hits =
    result.track.GetHits();
  
  if(hits.empty()){
    return;
  }
  
  const auto& firstFiber =
    geometry_.GetFiber(
		       hits.front().GetFiberID()
		       );
  
  const auto& lastFiber =
    geometry_.GetFiber(
		       hits.back().GetFiberID()
		       );
  
  
  int lastHitLayer = -1;
  
  if(isXZ){
    lastHitLayer =
      2 * lastFiber.layerIndex;
    
    if(lastFiber.layerType ==
       AFTLayerType::XPrime){
      ++lastHitLayer;
    }
    
    const double deltaX =
      lastFiber.center.X() -
      firstFiber.center.X();
    
    hEndHitVsDeltaXZ_->Fill(
			    lastHitLayer,
			    deltaX
			    );
    hEndHitPositionXZ_->Fill(
			     lastHitLayer,			     
			     lastFiber.center.X()			    
			     );
  }
  else{
    lastHitLayer =
      2 * lastFiber.layerIndex;
    
    if(lastFiber.layerType ==
       AFTLayerType::YPrime){
      ++lastHitLayer;
    }

    const double deltaY =
      lastFiber.center.Y() -
      firstFiber.center.Y();
    
    hEndHitVsDeltaYZ_->Fill(
			    lastHitLayer,
			    deltaY
			    );

    hEndHitPositionYZ_->Fill(
			     lastHitLayer,			     
			     lastFiber.center.Y()
			     );

  }
}

void
AFTTrackEvaluator::Process(
    const AFTEvent& event,
    const std::vector<AFTTrackResult>& tracks)
{
  (void)event;

  ++nEvents_;

  int nXZ = 0;
  int nYZ = 0;

  for(const auto& result : tracks){
    if(result.projection ==
       AFTProjection::XZ){
      ++nXZ;
    }
    else{
      ++nYZ;
    }

    FillCandidate(result);
  }

  nCandidatesXZ_->Fill(nXZ);
  nCandidatesYZ_->Fill(nYZ);

  nCandidatesTotal_->Fill(
    static_cast<int>(tracks.size()));

  if(nXZ > 0){
    ++nEventsWithXZ_;
  }

  if(nYZ > 0){
    ++nEventsWithYZ_;
  }

  if(nXZ > 0 && nYZ > 0){
    ++nEventsWithBoth_;
  }

}

void
AFTTrackEvaluator::PrintSummary() const
{
  const double nEvents =
    static_cast<double>(nEvents_);

  std::cout
    << "========== AFT Candidate Evaluation =========="
    << '\n'
    << "Events                  : "
    << nEvents_
    << '\n'
    << "Events with XZ candidate: "
    << nEventsWithXZ_;

  if(nEvents_ > 0){
    std::cout
      << " ("
      << 100.0 *
         static_cast<double>(nEventsWithXZ_) /
         nEvents
      << " %)";
  }

  std::cout
    << '\n'
    << "Events with YZ candidate: "
    << nEventsWithYZ_;

  if(nEvents_ > 0){
    std::cout
      << " ("
      << 100.0 *
         static_cast<double>(nEventsWithYZ_) /
         nEvents
      << " %)";
  }

  std::cout
    << '\n'
    << "Events with both         : "
    << nEventsWithBoth_;

  if(nEvents_ > 0){
    std::cout
      << " ("
      << 100.0 *
         static_cast<double>(nEventsWithBoth_) /
         nEvents
      << " %)";
  }

  std::cout
    << '\n'
    << "=============================================="
    << std::endl;
}

bool
AFTTrackEvaluator::Write(
    const std::string& filename) const
{
  TFile file(
    filename.c_str(),
    "RECREATE");

  if(file.IsZombie()){
    return false;
  }

  for(auto* hist :
      {nCandidatesXZ_,
       nCandidatesYZ_,
       nCandidatesTotal_,
       nHitsXZ_,
       nHitsYZ_,
       startXZ_,
       startYZ_,
       endXZ_,
       endYZ_,
       missingXZ_,
       missingYZ_,
       residualXZ_,
       residualYZ_}){

    hist->Write();
  }
  for(auto* hist :
	{hEndHitVsDeltaXZ_,
	 hEndHitVsDeltaYZ_,
	 hEndHitPositionXZ_,
	 hEndHitPositionYZ_}){
    
    hist->Write();
  }
  file.Close();

  return true;
}
