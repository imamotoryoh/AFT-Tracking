#ifndef EVENT_HH
#define EVENT_HH

#include "AFTEvent.hh"

class Event
{
public:
  Event() = default;

  AFTEvent& GetAFTEvent();
  const AFTEvent& GetAFTEvent() const;

  void Clear();

private:
  AFTEvent aftEvent_; 
};

#endif
