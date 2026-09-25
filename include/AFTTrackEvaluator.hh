#ifndef AFT_TRACK_EVALUATOR_HH
#define AFT_TRACK_EVALUATOR_HH

#include <cstddef>
#include <string>
#include <vector>

#include "AFTEvent.hh"
#include "AFTTrackResult.hh"
#include "AFTGeometry.hh"

class TH1D;
class TH2D;

class AFTTrackEvaluator
{
public:
  explicit AFTTrackEvaluator(
			     const AFTGeometry& geometry
			     );
  ~AFTTrackEvaluator();

  void Process(
      const AFTEvent& event,
      const std::vector<AFTTrackResult>& tracks);

  void PrintSummary() const;

  bool Write(
      const std::string& filename) const;

  std::size_t GetNEvents() const
  {
    return nEvents_;
  }

  std::size_t GetNEventsWithXZCandidate() const
  {
    return nEventsWithXZ_;
  }

  std::size_t GetNEventsWithYZCandidate() const
  {
    return nEventsWithYZ_;
  }

  std::size_t GetNEventsWithBothCandidates() const
  {
    return nEventsWithBoth_;
  }

private:
  void FillCandidate(
      const AFTTrackResult& result);

  std::size_t nEvents_ = 0;
  std::size_t nEventsWithXZ_ = 0;
  std::size_t nEventsWithYZ_ = 0;
  std::size_t nEventsWithBoth_ = 0;

  TH1D* nCandidatesXZ_ = nullptr;
  TH1D* nCandidatesYZ_ = nullptr;
  TH1D* nCandidatesTotal_ = nullptr;

  TH1D* nHitsXZ_ = nullptr;
  TH1D* nHitsYZ_ = nullptr;

  TH1D* startXZ_ = nullptr;
  TH1D* startYZ_ = nullptr;

  TH1D* endXZ_ = nullptr;
  TH1D* endYZ_ = nullptr;

  TH1D* missingXZ_ = nullptr;
  TH1D* missingYZ_ = nullptr;

  TH1D* residualXZ_ = nullptr;
  TH1D* residualYZ_ = nullptr;

  TH2D* hEndHitVsDeltaXZ_ = nullptr;
  TH2D* hEndHitVsDeltaYZ_ = nullptr;

  TH2D* hEndHitPositionXZ_ = nullptr;
  TH2D* hEndHitPositionYZ_ = nullptr;

  const AFTGeometry& geometry_;
  
};

#endif
