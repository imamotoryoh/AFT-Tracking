#include "AFTGeometry.hh"

#include <TApplication.h>
#include <TCanvas.h>
#include <TH2Poly.h>
#include <TMath.h>
#include <TPad.h>

#include <iostream>
#include <vector>

namespace
{
constexpr int NPoly = 100;

bool IsXZ(AFTLayerType type)
{
  return type == AFTLayerType::X ||
         type == AFTLayerType::XPrime;
}

bool IsYZ(AFTLayerType type)
{
  return type == AFTLayerType::Y ||
         type == AFTLayerType::YPrime;
}

void AddFiberBin(TH2Poly* hist,
                 double transverse,
                 double z,
                 double radius)
{
  double pos[NPoly];
  double zz[NPoly];

  // Same convention as the old EventDisp2DwithAFTKK.cc:
  // horizontal axis = x (or y), vertical axis = z.
  for(int i = 0; i < NPoly; ++i){
    const double theta =
      i * 2.0 * TMath::Pi() / static_cast<double>(NPoly);

    pos[i] = transverse + radius * TMath::Cos(theta);
    zz[i]  = z          + radius * TMath::Sin(theta);
  }

  hist->AddBin(NPoly, pos, zz);
}
}

int main(int argc, char** argv)
{
  if(argc != 2){
    std::cerr
      << "Usage: " << argv[0] << " AFTGeometry.dat"
      << std::endl;
    return 1;
  }

  AFTGeometry geometry;

  if(!geometry.Load(argv[1])){
    std::cerr
      << "Cannot open geometry file: " << argv[1]
      << std::endl;
    return 1;
  }

  std::cout
    << "Loaded fibers = " << geometry.GetNFibers()
    << std::endl;

  TApplication app("app", &argc, argv);

  // Same canvas/pad layout as the old EventDisplay:
  // XZ on the left, YZ on the right.
  TCanvas* c1 =
    new TCanvas("c1", "AFT Geometry", 1400, 900);

  c1->Divide(2, 1);

  TPad* p1 = static_cast<TPad*>(c1->cd(1));
  TPad* p2 = static_cast<TPad*>(c1->cd(2));

  TH2Poly* h_geo_x = new TH2Poly();
  TH2Poly* h_geo_y = new TH2Poly();

  // Build XZ and YZ fiber polygons directly from AFTGeometry.
  for(int fiberID = 0;
      fiberID < static_cast<int>(geometry.GetNFibers());
      ++fiberID){

    const auto& fiber = geometry.GetFiber(fiberID);

    if(IsXZ(fiber.layerType)){
      AddFiberBin(h_geo_x,
                  fiber.center.X(),
                  fiber.center.Z(),
                  fiber.radius);
    }

    if(IsYZ(fiber.layerType)){
      AddFiberBin(h_geo_y,
                  fiber.center.Y(),
                  fiber.center.Z(),
                 fiber.radius);
    }
  }

  // Keep the same viewing range and axis convention as the old macro.
  h_geo_x->SetTitle("AFT Geometry XZ;x [mm];z [mm]");
  h_geo_y->SetTitle("AFT Geometry YZ;y [mm];z [mm]");

  h_geo_x->GetXaxis()->SetLimits(-60.0, 60.0);
  h_geo_x->GetYaxis()->SetLimits(320.0, 460.0);

  h_geo_y->GetXaxis()->SetLimits(-60.0, 60.0);
  h_geo_y->GetYaxis()->SetLimits(320.0, 460.0);

  // Geometry only: give every physical fiber the same bin content.
  for(int bin = 1; bin <= h_geo_x->GetNumberOfBins(); ++bin){
    h_geo_x->SetBinContent(bin, 1.0);
  }

  for(int bin = 1; bin <= h_geo_y->GetNumberOfBins(); ++bin){
    h_geo_y->SetBinContent(bin, 1.0);
  }

  p1->cd();
  h_geo_x->SetFillStyle(0);
  h_geo_x->SetLineColor(kBlack);
  h_geo_x->SetLineWidth(1);
  h_geo_x->Draw("L");

  p2->cd();
  h_geo_y->SetFillStyle(0);
  h_geo_y->SetLineColor(kBlack);
  h_geo_y->SetLineWidth(1);
  h_geo_y->Draw("L");

  c1->Modified();
  c1->Update();

  app.Run();
  return 0;
}
