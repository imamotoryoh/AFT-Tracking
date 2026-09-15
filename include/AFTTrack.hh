#ifndef AFT_TRACK_HH
#define AFT_TRACK_HH

class AFTTrack
{
public:
  size_t GetNHits() const;

  int GetHitIndexFromStart(size_t i) const;
  int GetHitIndexFromEnd(size_t i) const;

  const TVector3& GetStartPosition() const;
  const TVector3& GetEndPosition() const;

  const TVector3& GetDirection() const;
  const TVector3& GetXZDirection() const;
  const TVector3& GetYZDirection() const;

  double GetLength() const;
  double GetChi2() const;

private:
  std::vector<size_t> hitIndices_;

  TVector3 startPosition_;
  TVector3 endPosition_;

  TVector3 direction_;
  TVector3 xzDirection_;
  TVector3 yzDirection_;

  double length_;
  double chi2_;
};
#endif
