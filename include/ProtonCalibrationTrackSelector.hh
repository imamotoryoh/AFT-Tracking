#ifndef PROTON_CALIBRATION_TRACK_SELECTOR_HH
#define PROTON_CALIBRATION_TRACK_SELECTOR_HH

#include "AFTTrackSelector.hh"

class ProtonCalibrationTrackSelector
  : public AFTTrackSelector
{
public:
  struct Parameters
  {
    int maxCandidatesPerProjection = 2;
    int maxMissingLayers = 3;

    double maxDeltaX = 10.0;
    double maxDeltaY = 5.0;
  };

  ProtonCalibrationTrackSelector();

  explicit ProtonCalibrationTrackSelector(
      const Parameters& parameters);

  AFTTrackSelectionResult
  Select(
      const std::vector<AFTTrackResult>& tracks,
      const AFTGeometry& geometry) const override;

private:
  bool PassTrackCuts(
      const AFTTrackResult& result,
      const AFTGeometry& geometry) const;

  double GetDisplacement(
      const AFTTrackResult& result,
      const AFTGeometry& geometry) const;

  Parameters parameters_;
};

#endif
