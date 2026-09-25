#ifndef AFT_TRACK_SELECTOR_HH
#define AFT_TRACK_SELECTOR_HH

#include <vector>

#include "AFTGeometry.hh"
#include "AFTTrackResult.hh"

struct AFTTrackSelectionResult
{
  const AFTTrackResult* trackXZ = nullptr;
  const AFTTrackResult* trackYZ = nullptr;

  bool IsSelected() const
  {
    return trackXZ != nullptr &&
           trackYZ != nullptr;
  }
};

class AFTTrackSelector
{
public:
  virtual ~AFTTrackSelector() = default;

  virtual AFTTrackSelectionResult
  Select(
      const std::vector<AFTTrackResult>& tracks,
      const AFTGeometry& geometry) const = 0;
};

#endif
