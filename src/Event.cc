#include "Event.hh"

AFTEvent& Event::GetAFTEvent()
{
  return aftEvent_;
}

const AFTEvent& Event::GetAFTEvent() const
{
  return aftEvent_;
}

void Event::Clear()
{
  aftEvent_.Clear();
}
