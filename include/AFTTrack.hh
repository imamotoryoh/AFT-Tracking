#ifndef AFT_TRACK_HH
#define AFT_TRACK_HH

#include <cstddef>
#include <vector>

#include "AFTHit.hh"

class AFTTrack
{
public:
  AFTTrack() = default;

  void AddHit(const AFTHit& hit);
  void Clear();

  const std::vector<AFTHit>& GetHits() const;

  const AFTHit& GetHit(std::size_t index) const;
  const AFTHit& GetHitFromEnd(std::size_t index) const;

  std::size_t GetNHits() const;

private:
  std::vector<AFTHit> hits_;

};
#endif
