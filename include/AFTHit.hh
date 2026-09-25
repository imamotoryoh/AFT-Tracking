#ifndef AFT_HIT_HH
#define AFT_HIT_HH

class AFTHit
{
public:
  AFTHit(int fiberID, double energy)
    : fiberID_(fiberID),
      energy_(energy)
  {}
  
  int GetFiberID() const {return fiberID_ ; }
  double GetEnergy() const {return energy_; }
private:
  int fiberID_;
  double energy_; // with tot cut and time cut, calibrated by 1.8GeV/c Kaon
}; 

#endif
