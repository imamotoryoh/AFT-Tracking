#ifndef AFT_EVENT_HH
#define AFT_EVENT_HH

#include <cstddef>
#include <vector>

#include "AFTHit.hh"

class AFTEvent
{
public:
  AFTEvent() = default;

  void AddHit(const AFTHit& hit);
  void Clear();

  const std::vector<AFTHit>& GetHits() const;

  bool HasHit(int fiberID) const;
  const AFTHit* FindHit(int fiberID) const;

  std::size_t GetNHits() const;

private:
  std::vector<AFTHit> hits_;
};

#endif
