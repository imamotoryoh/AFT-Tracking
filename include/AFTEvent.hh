#ifndef AFT_EVENT_HH
#define AFT_EVENT_HH

class AFTEvent
{
public:
  int GetEvetID() const;
  
  const std::vector<AFTHit>& GetHits() const;

  void AddHit(const AFTHit& hit);
  void Clear();
  
private:
  int eventID_;
  
  std::vector<AFTHit> hits_;
};
#endif
