#include "AFTEvent.hh"

void AFTEvent::AddHit(const AFTHit& hit)
{
  hits_.push_back(hit);
}

void AFTEvent::Clear()
{
  hits_.clear();
}

const std::vector<AFTHit>& AFTEvent::GetHits() const
{
  return hits_;
}

bool AFTEvent::HasHit(int fiberID) const
{
  for(const auto& hit : hits_){
    if(hit.GetFiberID() == fiberID){
      return true;
    }
  }
  return false;
}

const AFTHit* AFTEvent::FindHit(int fiberID) const
{
  for(const auto& hit : hits_){
    if(hit.GetFiberID() == fiberID){
      return &hit;
    }
  }
  return nullptr;
}

std::size_t AFTEvent::GetNHits() const
{
  return hits_.size();
}
