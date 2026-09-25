#include "ProtonCalibrationTrackSelector.hh"

#include <cmath>
#include <vector>

ProtonCalibrationTrackSelector::
ProtonCalibrationTrackSelector()
  : parameters_()
{
}

ProtonCalibrationTrackSelector::
ProtonCalibrationTrackSelector(
    const Parameters& parameters)
  : parameters_(parameters)
{
}

AFTTrackSelectionResult
ProtonCalibrationTrackSelector::Select(
    const std::vector<AFTTrackResult>& tracks,
    const AFTGeometry& geometry) const
{
  std::vector<const AFTTrackResult*> candidatesXZ;
  std::vector<const AFTTrackResult*> candidatesYZ;

  for(const auto& result : tracks){
    if(result.projection == AFTProjection::XZ){
      candidatesXZ.push_back(&result);
    }
    else if(result.projection == AFTProjection::YZ){
      candidatesYZ.push_back(&result);
    }
  }

  AFTTrackSelectionResult selection;

  // Reject events with too many original candidates.
  if(candidatesXZ.empty() ||
     candidatesYZ.empty()){
    return selection;
  }

  if(candidatesXZ.size() >
       static_cast<std::size_t>(
         parameters_.maxCandidatesPerProjection) ||
     candidatesYZ.size() >
       static_cast<std::size_t>(
         parameters_.maxCandidatesPerProjection)){
    return selection;
  }

  std::vector<const AFTTrackResult*> selectedXZ;
  std::vector<const AFTTrackResult*> selectedYZ;

  for(const auto* result : candidatesXZ){
    if(PassTrackCuts(*result, geometry)){
      selectedXZ.push_back(result);
    }
  }

  for(const auto* result : candidatesYZ){
    if(PassTrackCuts(*result, geometry)){
      selectedYZ.push_back(result);
    }
  }

  // Calibration sample:
  // require one and only one surviving track
  // in each projection.
  if(selectedXZ.size() != 1 ||
     selectedYZ.size() != 1){
    return selection;
  }

  selection.trackXZ = selectedXZ.front();
  selection.trackYZ = selectedYZ.front();

  return selection;
}

bool ProtonCalibrationTrackSelector::PassTrackCuts(
    const AFTTrackResult& result,
    const AFTGeometry& geometry) const
{
  if(result.nMissing >
     parameters_.maxMissingLayers){
    return false;
  }

  if(result.track.GetNHits() < 2){
    return false;
  }

  const double displacement =
    std::abs(
      GetDisplacement(result, geometry)
    );

  if(result.projection == AFTProjection::XZ){
    if(displacement > parameters_.maxDeltaX){
      return false;
    }
  }
  else if(result.projection == AFTProjection::YZ){
    if(displacement > parameters_.maxDeltaY){
      return false;
    }
  }

  return true;
}

double ProtonCalibrationTrackSelector::GetDisplacement(
    const AFTTrackResult& result,
    const AFTGeometry& geometry) const
{
  const auto& hits =
    result.track.GetHits();

  const auto& firstFiber =
    geometry.GetFiber(
      hits.front().GetFiberID()
    );

  const auto& lastFiber =
    geometry.GetFiber(
      hits.back().GetFiberID()
    );

  if(result.projection == AFTProjection::XZ){
    return lastFiber.center.X() -
           firstFiber.center.X();
  }

  return lastFiber.center.Y() -
         firstFiber.center.Y();
}
