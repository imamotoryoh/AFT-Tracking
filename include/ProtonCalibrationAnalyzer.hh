#ifndef PROTON_CALIBRATION_ANALYZER_HH
#define PROTON_CALIBRATION_ANALYZER_HH

#include <cstddef>
#include <string>
#include <array>

#include <TH1D.h>
#include <TH2D.h>

#include "AFTGeometry.hh"
#include "AFTTrackSelector.hh"

class ProtonCalibrationAnalyzer
{
public:
  explicit ProtonCalibrationAnalyzer(
				     const AFTGeometry& geometry,
				     const double calib);

  ~ProtonCalibrationAnalyzer();

  void Process(
      const AFTTrackSelectionResult& selection);

  bool Write(
      const std::string& filename) const;

private:
  static constexpr int NGlobalLayers = 36;
  void ProcessTrack(
      const AFTTrackResult& result);

  void FillHit(
      const AFTHit& hit,
      std::size_t indexFromEnd);

  const AFTGeometry& geometry_;

  TH1D* hFirstHitEnergy_;
  TH1D* hSecondHitEnergy_;
  TH1D* hThirdHitEnergy_;

  TH2D* hFirstHitEnergyVsFiberID_;
  TH2D* hSecondHitEnergyVsFiberID_;
  TH2D* hThirdHitEnergyVsFiberID_;

  std::array<TH2D*, NGlobalLayers>
  hFirstHitEnergyVsFiberIDForEachLayer_;
  
  std::array<TH2D*, NGlobalLayers>
  hSecondHitEnergyVsFiberIDForEachLayer_;
  
  std::array<TH2D*, NGlobalLayers>
  hThirdHitEnergyVsFiberIDForEachLayer_;

  double calib_;  
};

#endif
