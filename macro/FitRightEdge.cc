#include <iostream>
#include <cmath>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLine.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>


int GetFirstFiberID(int globalLayerID)
{
  const int layerIndex = globalLayerID / 4;
  const int layerType  = globalLayerID % 4;

  int offset = 0;

  switch(layerType){
  case 0: // X
    offset = 0;
    break;

  case 1: // X'
    offset = 32;
    break;

  case 2: // Y
    offset = 64;
    break;

  case 3: // Y'
    offset = 80;
    break;
  }

  return 96 * layerIndex + offset;
}


int GetNFibersInLayer(int globalLayerID)
{
  const int layerType = globalLayerID % 4;

  if(layerType == 0 || layerType == 1)
    return 32;

  return 16;
}


double RightEdgeFunction(double* x, double* par)
{
  const double E = x[0];

  const double C     = par[0];
  const double A     = par[1];
  const double E0    = par[2];
  const double sigma = par[3];

  const double step =
    0.5 * A *
    (1.0 - TMath::Erf(
      (E - E0) /
      (std::sqrt(2.0) * sigma)
    ));

  return C + step;
}


void FitRightEdge(const char* filename,
                  int layer,
                  int fiberIndex,
		  double initial,
                  double fitMin = 8.0,
                  double fitMax = 10.0)
{
  const int nFibers =
    GetNFibersInLayer(layer);

  if(fiberIndex < 0 ||
     fiberIndex >= nFibers){
    std::cerr
      << "Invalid fiber index: "
      << fiberIndex
      << std::endl;
    return;
  }


  // --------------------------------------------------
  // FiberID
  // --------------------------------------------------

  const int firstFiberID =
    GetFirstFiberID(layer);

  const int fiberID =
    firstFiberID + fiberIndex;


  // --------------------------------------------------
  // Open ROOT file
  // --------------------------------------------------

  TFile* file =
    TFile::Open(filename, "READ");

  if(!file || file->IsZombie()){
    std::cerr
      << "Failed to open file: "
      << filename
      << std::endl;
    return;
  }


  // --------------------------------------------------
  // Get TH2D
  // --------------------------------------------------

  TString histName =
    Form("hFirstHitEnergyVsFiberID_Layer%d",
         layer);

  TH2D* hist2D = nullptr;
  file->GetObject(histName, hist2D);

  if(!hist2D){
    std::cerr
      << "Histogram not found: "
      << histName
      << std::endl;

    file->Close();
    return;
  }


  // --------------------------------------------------
  // Find FiberID bin
  // --------------------------------------------------

  const int binX =
    hist2D->GetXaxis()->FindBin(fiberID);

  const double binCenter =
    hist2D->GetXaxis()->GetBinCenter(binX);

  std::cout
    << "FiberID = " << fiberID
    << ", X bin = " << binX
    << ", bin center = " << binCenter
    << std::endl;


  // --------------------------------------------------
  // Projection
  // --------------------------------------------------

  TString projectionName =
    Form("hEnergy_Layer%d_Fiber%d",
         layer,
         fiberID);

  TH1D* hist =
    hist2D->ProjectionY(
      projectionName,
      binX,
      binX
    );

  hist->SetDirectory(nullptr);

  hist->SetTitle(
    Form("Layer %d, FiberID %d;Energy [MeV];Counts",
         layer,
         fiberID)
  );


  // --------------------------------------------------
  // Fit function
  //
  // par[0] : C
  // par[1] : A
  // par[2] : E0
  // par[3] : sigma
  // --------------------------------------------------

  TF1* func =
    new TF1(
      "fRightEdge",
      RightEdgeFunction,
      fitMin,
      fitMax,
      4
    );

  func->SetParNames(
    "C",
    "A",
    "E0",
    "sigma"
  );


  // --------------------------------------------------
  // Initial parameters
  // --------------------------------------------------

  const double E0Initial = initial;

  const int binE0 =
    hist->GetXaxis()->FindBin(E0Initial);

  const double AInitial =
    hist->GetBinContent(binE0);

  func->SetParameters(
    0.0,       // C
    AInitial,  // A
    E0Initial, // E0
    0.15       // sigma
  );


  // --------------------------------------------------
  // Parameter limits
  // --------------------------------------------------

  func->SetParLimits(
    0,
    0.0,
    hist->GetMaximum()
  );

  func->SetParLimits(
    1,
    0.0,
    10.0 * hist->GetMaximum()
  );

  func->SetParLimits(
    2,
    fitMin,
    fitMax
  );

  func->SetParLimits(
    3,
    0.001,
    2.0
  );


  // --------------------------------------------------
  // Fit
  // R : fit range
  // L : likelihood
  // S : TFitResultPtr
  // --------------------------------------------------

  TFitResultPtr fitResult =
    hist->Fit(
      func,
      "RLS"
    );


  // --------------------------------------------------
  // Print result
  // --------------------------------------------------

  std::cout
    << "\n========== Fit result =========="
    << std::endl;

  std::cout
    << "Layer    = " << layer
    << std::endl;

  std::cout
    << "FiberID  = " << fiberID
    << std::endl;

  std::cout
    << "E0       = "
    << func->GetParameter(2)
    << " +/- "
    << func->GetParError(2)
    << " MeV"
    << std::endl;

  std::cout
    << "sigma    = "
    << func->GetParameter(3)
    << " +/- "
    << func->GetParError(3)
    << " MeV"
    << std::endl;

  std::cout
    << "A        = "
    << func->GetParameter(1)
    << std::endl;

  std::cout
    << "C        = "
    << func->GetParameter(0)
    << std::endl;

  std::cout
    << "Fit status = "
    << static_cast<int>(fitResult)
    << std::endl;

  std::cout
    << "================================"
    << std::endl;


  // --------------------------------------------------
  // Draw
  // --------------------------------------------------

  TCanvas* canvas =
    new TCanvas(
      "cFitRightEdge",
      "Right edge fit",
      900,
      700
    );

  hist->Draw("E");
  func->Draw("same");

  const double E0 =
    func->GetParameter(2);

  const double yMin = 0.0;
  const double yMax =
    hist->GetMaximum() * 1.1;

  TLine* edgeLine =
    new TLine(
      E0,
      yMin,
      E0,
      yMax
    );

  edgeLine->SetLineStyle(2);
  edgeLine->SetLineWidth(2);
  edgeLine->Draw("same");

  canvas->Update();

  file->Close();
}
