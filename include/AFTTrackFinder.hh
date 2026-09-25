#ifndef AFT_TRACK_FINDER_HH
#define AFT_TRACK_FINDER_HH

#include <vector>

#include "AFTEvent.hh"
#include "AFTGeometry.hh"
#include "AFTTrackResult.hh"

class AFTTrackFinder
{
public:
  virtual ~AFTTrackFinder() = default;

  virtual std::vector<AFTTrackResult>
  FindTracks(const AFTEvent& event,
             const AFTGeometry& geometry) const = 0;
};

#endif
