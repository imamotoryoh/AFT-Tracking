#ifndef AFT_EVENT_DISPLAY_HH
#define AFT_EVENT_DISPLAY_HH

#include <vector>

#include "AFTEvent.hh"
#include "AFTGeometry.hh"
#include "AFTTrackResult.hh"

class TCanvas;
class TPad;
class TH2Poly;
class TPolyLine;

class AFTEventDisplay
{
public:
  explicit AFTEventDisplay(
      const AFTGeometry& geometry);

  ~AFTEventDisplay();

  void DrawGeometry();

  void DrawEvent(
      const AFTEvent& event);

  void DrawEvent(
      const AFTEvent& event,
      const std::vector<AFTTrackResult>& tracks);

  TCanvas* GetCanvas() const
  {
    return canvas_;
  }

private:
  void Build();
  void Reset();
  void Fill(const AFTEvent& event);
  void DrawTrack(const AFTTrackResult& result);

  const AFTGeometry& geometry_;

  TCanvas* canvas_ = nullptr;

  TPad* padXZ_ = nullptr;
  TPad* padYZ_ = nullptr;

  TH2Poly* geoXZ_ = nullptr;
  TH2Poly* geoYZ_ = nullptr;

  TH2Poly* evtXZ_ = nullptr;
  TH2Poly* evtYZ_ = nullptr;

  std::vector<TPolyLine*> lines_;
};

#endif
