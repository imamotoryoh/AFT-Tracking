#ifndef AFT_HIT_HH
#define AFT_HIT_HH

class AFTHit
{
public:
  int GetFiberID() const;
  double GetAdc() const;
  double GetTime() const;
private:
  int fiberID_;

  double adc_;
  double time_;
};

#endif
