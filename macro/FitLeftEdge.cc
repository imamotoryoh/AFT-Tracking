#include <iostream>
#include <cmath>

#include <TLine.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

double LeftEdgeFunction(double* x, double* par)
{
  const double E = x[0];

  const double C      = par[0];
  const double mLeft  = par[1];
  const double mRight = par[2];
  const double A      = par[3];
  const double E0     = par[4];
  const double sigma  = par[5];

  // Background:
  // continuous at E0, but its slope changes at E0.
  const double slope =
    (E < E0) ? mLeft : mRight;

  const double bg =
    C + slope * (E - E0);

  // Smeared step
  const double step =
    0.5 * A *
    (1.0 + TMath::Erf(
      (E - E0) /
      (std::sqrt(2.0) * sigma)
    ));

  return bg + step;
}

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

void FitLeftEdge(const char* filename,
                 int layer,
                 int fiberIndex,
		 double initial,
                 double fitMin = 3.0,
                 double fitMax = 7.0)
{
  // --------------------------------------------------
  // Open ROOT file
  // --------------------------------------------------

  TFile* file = TFile::Open(filename, "READ");

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
    Form("hFirstHitEnergyVsFiberID_Layer%d", layer);

  const int firstFiberID = GetFirstFiberID(layer);
  const int fiberID = firstFiberID + fiberIndex;
    
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
  // par[1] : mLeft
  // par[2] : mRight
  // par[3] : A
  // par[4] : E0
  // par[5] : sigma
  // --------------------------------------------------

  TF1* func =
    new TF1(
      "fLeftEdge",
      LeftEdgeFunction,
      fitMin,
      fitMax,
      6
    );

  func->SetParNames(
    "C",
    "mLeft",
    "mRight",
    "A",
    "E0",
    "sigma"
  );


  // --------------------------------------------------
  // Initial parameters
  //
  // These are only starting values.
  // --------------------------------------------------

  const double E0Initial = initial;

  const int binE0 =
    hist->GetXaxis()->FindBin(E0Initial);

  const double CInitial =
    hist->GetBinContent(binE0);

  func->SetParameters(
    CInitial,  // C
    10.0,      // mLeft
    -20.0,     // mRight
    100.0,     // A
    E0Initial, // E0
    0.15       // sigma
  );


  // Keep the important parameters in reasonable regions.
  func->SetParLimits(
    3,
    0.0,
    10000.0
  );

  func->SetParLimits(
    4,
    fitMin,
    fitMax
  );

  func->SetParLimits(
    5,
    0.001,
    2.0
  );


  // --------------------------------------------------
  // Fit
  // --------------------------------------------------

  TFitResultPtr fitResult =
    hist->Fit(
      func,
      "RS"
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
    << func->GetParameter(4)
    << " +/- "
    << func->GetParError(4)
    << " MeV"
    << std::endl;

  std::cout
    << "sigma    = "
    << func->GetParameter(5)
    << " +/- "
    << func->GetParError(5)
    << " MeV"
    << std::endl;

  std::cout
    << "mLeft    = "
    << func->GetParameter(1)
    << std::endl;

  std::cout
    << "mRight   = "
    << func->GetParameter(2)
    << std::endl;

  std::cout
    << "A        = "
    << func->GetParameter(3)
    << std::endl;

  if(func->GetNDF() > 0){
    std::cout
      << "chi2/NDF = "
      << func->GetChisquare()
      / func->GetNDF()
      << std::endl;
  }

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
      "cFitLeftEdge",
      "Left edge fit",
      900,
      700
    );

  hist->Draw("E");
  func->Draw("same");

  const double E0 = func->GetParameter(4);

  const double yMin = 0.0;
  const double yMax = hist->GetMaximum() * 1.1;
  
  TLine* edgeLine =
    new TLine(E0, yMin, E0, yMax);
  
  edgeLine->SetLineStyle(2);
  edgeLine->SetLineWidth(2);
  edgeLine->Draw("same");


  hist->GetXaxis()->SetRangeUser(0,8);
  canvas->Update();


  // ROOT objects needed for display are intentionally
  // left alive here.

  file->Close();
}
