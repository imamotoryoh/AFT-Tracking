#include "AFTEventDisplay.hh"

#include <cmath>

#include <TCanvas.h>
#include <TH2Poly.h>
#include <TMath.h>
#include <TPad.h>
#include <TPolyLine.h>

namespace
{

constexpr int NPolygonPoints = 100;

constexpr double TransverseMin = -60.0;
constexpr double TransverseMax =  60.0;
constexpr double ZMin = 320.0;
constexpr double ZMax = 460.0;

bool IsXZ(AFTLayerType type)
{
  return type == AFTLayerType::X ||
         type == AFTLayerType::XPrime;
}

void AddFiberBin(TH2Poly* hist,
                 double transverse,
                 double z,
                 double radius)
{
  double x[NPolygonPoints + 1];
  double y[NPolygonPoints + 1];

  for(int i = 0; i < NPolygonPoints; ++i){
    const double angle =
      2.0 * TMath::Pi()
      * static_cast<double>(i)
      / static_cast<double>(NPolygonPoints);

    x[i] = transverse + radius * std::cos(angle);
    y[i] = z          + radius * std::sin(angle);
  }

  x[NPolygonPoints] = x[0];
  y[NPolygonPoints] = y[0];

  hist->AddBin(NPolygonPoints + 1, x, y);
}

}

AFTEventDisplay::AFTEventDisplay(
    const AFTGeometry& geometry)
  : geometry_(geometry)
{
  canvas_ =
    new TCanvas("cAFTEventDisplay",
                "AFT Event Display",
                900,
                600);

  canvas_->Divide(2, 1);

  padXZ_ =
    static_cast<TPad*>(canvas_->cd(1));

  padYZ_ =
    static_cast<TPad*>(canvas_->cd(2));

  Build();
}

AFTEventDisplay::~AFTEventDisplay()
{
  for(auto* line : lines_){
    delete line;
  }

  delete evtXZ_;
  delete evtYZ_;

  delete geoXZ_;
  delete geoYZ_;

  delete canvas_;
}

void
AFTEventDisplay::Build()
{
  geoXZ_ = new TH2Poly();
  geoYZ_ = new TH2Poly();

  evtXZ_ = new TH2Poly();
  evtYZ_ = new TH2Poly();

  geoXZ_->SetStats(0);
  geoYZ_->SetStats(0);
  evtXZ_->SetStats(0);
  evtYZ_->SetStats(0);

  for(int fiberID = 0;
      fiberID < static_cast<int>(geometry_.GetNFibers());
      ++fiberID){

    const auto& fiber =
      geometry_.GetFiber(fiberID);

    if(IsXZ(fiber.layerType)){
      AddFiberBin(geoXZ_,
                  fiber.center.X(),
                  fiber.center.Z(),
                  fiber.radius);

      AddFiberBin(evtXZ_,
                  fiber.center.X(),
                  fiber.center.Z(),
                  fiber.radius);
    }
    else{
      AddFiberBin(geoYZ_,
                  fiber.center.Y(),
                  fiber.center.Z(),
                  fiber.radius);

      AddFiberBin(evtYZ_,
                  fiber.center.Y(),
                  fiber.center.Z(),
                  fiber.radius);
    }
  }

  geoXZ_->SetFillStyle(0);
  geoYZ_->SetFillStyle(0);

  evtXZ_->SetMinimum(1.0e-10);
  evtYZ_->SetMinimum(1.0e-10);
}

void
AFTEventDisplay::DrawGeometry()
{
  Reset();

  padXZ_->cd();
  padXZ_->Clear();

  padXZ_->DrawFrame(TransverseMin,
                    ZMin,
                    TransverseMax,
                    ZMax,
                    "AFT XZ;x [mm];z [mm]");

  geoXZ_->Draw("L SAME");

  padYZ_->cd();
  padYZ_->Clear();

  padYZ_->DrawFrame(TransverseMin,
                    ZMin,
                    TransverseMax,
                    ZMax,
                    "AFT YZ;y [mm];z [mm]");

  geoYZ_->Draw("L SAME");

  canvas_->Modified();
  canvas_->Update();
}

void
AFTEventDisplay::Reset()
{
  for(int bin = 1;
      bin <= evtXZ_->GetNumberOfBins();
      ++bin){
    evtXZ_->SetBinContent(bin, 0.0);
  }

  for(int bin = 1;
      bin <= evtYZ_->GetNumberOfBins();
      ++bin){
    evtYZ_->SetBinContent(bin, 0.0);
  }

  for(auto* line : lines_){
    delete line;
  }

  lines_.clear();
}

void
AFTEventDisplay::Fill(
    const AFTEvent& event)
{
  for(const auto& hit : event.GetHits()){
    const auto& fiber =
      geometry_.GetFiber(hit.GetFiberID());

    TH2Poly* hist = nullptr;
    double transverse = 0.0;

    if(IsXZ(fiber.layerType)){
      hist = evtXZ_;
      transverse = fiber.center.X();
    }
    else{
      hist = evtYZ_;
      transverse = fiber.center.Y();
    }

    const int bin =
      hist->FindBin(transverse,
                    fiber.center.Z());

    if(bin > 0){
      hist->SetBinContent(bin,
                          hit.GetEnergy());
    }
  }
}

void
AFTEventDisplay::DrawEvent(
    const AFTEvent& event)
{
  DrawEvent(event, {});
}

void
AFTEventDisplay::DrawEvent(
    const AFTEvent& event,
    const std::vector<AFTTrackResult>& tracks)
{
  Reset();
  Fill(event);

  padXZ_->cd();
  padXZ_->Clear();

  padXZ_->DrawFrame(TransverseMin,
                    ZMin,
                    TransverseMax,
                    ZMax,
                    "AFT XZ;x [mm];z [mm]");

  evtXZ_->Draw("COL SAME");
  geoXZ_->Draw("L SAME");

  padYZ_->cd();
  padYZ_->Clear();

  padYZ_->DrawFrame(TransverseMin,
                    ZMin,
                    TransverseMax,
                    ZMax,
                    "AFT YZ;y [mm];z [mm]");

  evtYZ_->Draw("COL SAME");
  geoYZ_->Draw("L SAME");

  for(const auto& result : tracks){
    DrawTrack(result);
  }

  canvas_->Modified();
  canvas_->Update();
}

void
AFTEventDisplay::DrawTrack(
    const AFTTrackResult& result)
{
  const auto& hits =
    result.track.GetHits();

  if(hits.size() < 2){
    return;
  }

  std::vector<double> transverse(
    hits.size());

  std::vector<double> z(
    hits.size());

  for(std::size_t i = 0;
      i < hits.size();
      ++i){

    const auto& fiber =
      geometry_.GetFiber(
        hits[i].GetFiberID());

    if(result.projection ==
       AFTProjection::XZ){
      transverse[i] =
        fiber.center.X();
    }
    else{
      transverse[i] =
        fiber.center.Y();
    }

    z[i] =
      fiber.center.Z();
  }

  TPad* pad =
    (result.projection == AFTProjection::XZ)
      ? padXZ_
      : padYZ_;

  pad->cd();

  auto* line =
    new TPolyLine(
      static_cast<int>(hits.size()),
      transverse.data(),
      z.data());

  line->SetLineWidth(2);
  line->Draw("L SAME");

  lines_.push_back(line);
}
