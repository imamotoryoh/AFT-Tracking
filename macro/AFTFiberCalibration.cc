#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TString.h>

struct RunInfo
{
  int runNumber;
  double momentum;
};

// Fill these by hand. Order is irrelevant.
const std::vector<RunInfo> kRuns = {
  // {74618, 0.00},
  // {...,   0.00},
};

const std::string kInputDirectory = ".";
const int kStoppedRunNumber = 74618;

double kLeftFitMin  = 3.0;
double kLeftFitMax  = 7.0;
double kRightFitMin = 8.0;
double kRightFitMax = 10.0;

struct LeftEdgeParameters
{
  double C = 0.0;
  double mLeft = 0.0;
  double mRight = 0.0;
  double A = 0.0;
  double E0 = 5.0;
  double sigma = 0.15;
};

struct RightEdgeParameters
{
  double C = 0.0;
  double A = 0.0;
  double E0 = 9.0;
  double sigma = 0.15;
};

struct EdgeFitResult
{
  bool valid = false;
  int status = -1;
  double position = 0.0;
  double positionError = 0.0;
  LeftEdgeParameters left;
  RightEdgeParameters right;
};

int GetFirstFiberID(int globalLayerID)
{
  const int layerIndex = globalLayerID / 4;
  const int layerType = globalLayerID % 4;

  int offset = 0;
  if(layerType == 1) offset = 32;
  if(layerType == 2) offset = 64;
  if(layerType == 3) offset = 80;

  return 96 * layerIndex + offset;
}

int GetNFibersInLayer(int globalLayerID)
{
  const int layerType = globalLayerID % 4;
  return (layerType == 0 || layerType == 1) ? 32 : 16;
}

std::string GetExpFileName(int runNumber)
{
  return Form("%s/run%d_AFTTrackCalibrationea0c.root",
              kInputDirectory.c_str(), runNumber);
}

std::string GetSimFileName(int runNumber)
{
  return Form("%s/run%d_AFTTrackCalibrationtree.root",
              kInputDirectory.c_str(), runNumber);
}

const RunInfo* FindRunInfo(int runNumber)
{
  for(const auto& run : kRuns){
    if(run.runNumber == runNumber)
      return &run;
  }
  return nullptr;
}

TH1D* GetFiberSpectrum(TFile& file,
                       int globalLayerID,
                       int fiberIndex,
                       const std::string& name,
                       const char* histogramNameFormat)
{
  TString histName =
    Form(histogramNameFormat, globalLayerID);

  TH2D* hist2D = nullptr;
  file.GetObject(histName, hist2D);

  if(!hist2D){
    std::cerr << "Histogram not found: " << histName
              << " in " << file.GetName() << std::endl;
    return nullptr;
  }

  const int fiberID =
    GetFirstFiberID(globalLayerID) + fiberIndex;

  const int binX =
    hist2D->GetXaxis()->FindBin(fiberID);

  TH1D* hist =
    hist2D->ProjectionY(name.c_str(), binX, binX);

  hist->SetDirectory(nullptr);
  return hist;
}

// Left edge:
// BG slope switches at E0; rising step is an error function.
double LeftEdgeFunction(double* x, double* par)
{
  const double E = x[0];
  const double C = par[0];
  const double mLeft = par[1];
  const double mRight = par[2];
  const double A = par[3];
  const double E0 = par[4];
  const double sigma = par[5];

  const double slope = (E < E0) ? mLeft : mRight;
  const double bg = C + slope * (E - E0);
  const double step =
    0.5 * A *
    (1.0 + TMath::Erf((E - E0) /
                      (std::sqrt(2.0) * sigma)));

  return bg + step;
}

// Right edge:
// C + falling error function.
double RightEdgeFunction(double* x, double* par)
{
  const double E = x[0];
  const double C = par[0];
  const double A = par[1];
  const double E0 = par[2];
  const double sigma = par[3];

  return C +
    0.5 * A *
    (1.0 - TMath::Erf((E - E0) /
                      (std::sqrt(2.0) * sigma)));
}

EdgeFitResult FitLeftEdge(TH1D& hist,
                          const LeftEdgeParameters& initial)
{
  EdgeFitResult result;

  TF1 func("fLeftEdge", LeftEdgeFunction,
           kLeftFitMin, kLeftFitMax, 6);

  func.SetParNames("C", "mLeft", "mRight", "A", "E0", "sigma");
  func.SetParameters(initial.C, initial.mLeft, initial.mRight,
                     initial.A, initial.E0, initial.sigma);

  func.SetParLimits(3, 0.0, 1.0e9);
  func.SetParLimits(4, kLeftFitMin, kLeftFitMax);
  func.SetParLimits(5, 0.001, 2.0);

  TFitResultPtr r = hist.Fit(&func, "RLSQ");

  result.status = static_cast<int>(r);
  result.valid = (result.status == 0);
  result.position = func.GetParameter(4);
  result.positionError = func.GetParError(4);

  result.left.C = func.GetParameter(0);
  result.left.mLeft = func.GetParameter(1);
  result.left.mRight = func.GetParameter(2);
  result.left.A = func.GetParameter(3);
  result.left.E0 = func.GetParameter(4);
  result.left.sigma = func.GetParameter(5);

  return result;
}

EdgeFitResult FitRightEdge(TH1D& hist,
                           const RightEdgeParameters& initial)
{
  EdgeFitResult result;

  TF1 func("fRightEdge", RightEdgeFunction,
           kRightFitMin, kRightFitMax, 4);

  func.SetParNames("C", "A", "E0", "sigma");
  func.SetParameters(initial.C, initial.A,
                     initial.E0, initial.sigma);

  func.SetParLimits(0, 0.0, 1.0e9);
  func.SetParLimits(1, 0.0, 1.0e9);
  func.SetParLimits(2, kRightFitMin, kRightFitMax);
  func.SetParLimits(3, 0.001, 2.0);

  TFitResultPtr r = hist.Fit(&func, "RLSQ");

  result.status = static_cast<int>(r);
  result.valid = (result.status == 0);
  result.position = func.GetParameter(2);
  result.positionError = func.GetParError(2);

  result.right.C = func.GetParameter(0);
  result.right.A = func.GetParameter(1);
  result.right.E0 = func.GetParameter(2);
  result.right.sigma = func.GetParameter(3);

  return result;
}

double FindMPV(const TH1D& hist)
{
  const int maxBin = hist.GetMaximumBin();
  return hist.GetXaxis()->GetBinCenter(maxBin);
}

// ============================================================
// Initial parameters for the first fiber
//
// EXP and SIM are intentionally configured independently.
// Tune these four blocks separately.
// ============================================================

// ---------- EXP : left edge ----------
LeftEdgeParameters InitialLeftExp()
{
  LeftEdgeParameters p;

  p.C      = 50.0;
  p.mLeft  = 10.0;
  p.mRight = -10.0;
  p.A      = 100.0;
  p.E0     = 5.0;
  p.sigma  = 0.15;

  return p;
}

// ---------- SIM : left edge ----------
LeftEdgeParameters InitialLeftSim()
{
  LeftEdgeParameters p;

  p.C      = 50.0;
  p.mLeft  = 10.0;
  p.mRight = -10.0;
  p.A      = 100.0;
  p.E0     = 5.0;
  p.sigma  = 0.15;

  return p;
}

// ---------- EXP : right edge ----------
RightEdgeParameters InitialRightExp()
{
  RightEdgeParameters p;

  p.C     = 1.0;
  p.A     = 100.0;
  p.E0    = 9.0;
  p.sigma = 0.15;

  return p;
}

// ---------- SIM : right edge ----------
RightEdgeParameters InitialRightSim()
{
  RightEdgeParameters p;

  p.C     = 1.0;
  p.A     = 100.0;
  p.E0    = 9.0;
  p.sigma = 0.15;

  return p;
}

void AFTFiberCalibration(int globalLayerID,
                         const char* outputFileName =
                           "AFTFiberCalibration.root")
{
  const RunInfo* stoppedRun =
    FindRunInfo(kStoppedRunNumber);

  if(!stoppedRun){
    std::cerr << "Stopped run " << kStoppedRunNumber
              << " is not in kRuns." << std::endl;
    return;
  }

  const double stoppedMomentum = stoppedRun->momentum;

  struct OpenRun
  {
    RunInfo info;
    std::unique_ptr<TFile> exp;
    std::unique_ptr<TFile> sim;
  };

  std::vector<OpenRun> runs;

  for(const auto& run : kRuns){
    if(run.momentum < stoppedMomentum)
      continue;

    auto exp = std::unique_ptr<TFile>(
      TFile::Open(GetExpFileName(run.runNumber).c_str(), "READ"));
    auto sim = std::unique_ptr<TFile>(
      TFile::Open(GetSimFileName(run.runNumber).c_str(), "READ"));

    if(!exp || exp->IsZombie() || !sim || sim->IsZombie()){
      std::cerr << "Failed to open run " << run.runNumber << std::endl;
      continue;
    }

    runs.push_back({run, std::move(exp), std::move(sim)});
  }

  TFile output(outputFileName, "RECREATE");
  if(output.IsZombie())
    return;

  LeftEdgeParameters leftExpInitial = InitialLeftExp();
  LeftEdgeParameters leftSimInitial = InitialLeftSim();
  RightEdgeParameters rightExpInitial = InitialRightExp();
  RightEdgeParameters rightSimInitial = InitialRightSim();

  const int nFibers = GetNFibersInLayer(globalLayerID);
  const int firstFiberID = GetFirstFiberID(globalLayerID);

  for(int fiberIndex = 0; fiberIndex < nFibers; ++fiberIndex){
    const int fiberID = firstFiberID + fiberIndex;

    auto graph = std::make_unique<TGraph>();
    graph->SetName(
      Form("gCalibration_Layer%d_Fiber%d",
           globalLayerID, fiberID));
    graph->SetTitle(
      Form("Layer %d FiberID %d;Simulation feature [MeV];Experimental feature [MeV]",
           globalLayerID, fiberID));

    for(auto& run : runs){
      const std::string expName =
        Form("hExp_L%d_F%d_Run%d",
             globalLayerID, fiberID, run.info.runNumber);
      const std::string simName =
        Form("hSim_L%d_F%d_Run%d",
             globalLayerID, fiberID, run.info.runNumber);

      // The stopped run uses the selected "first hit from track end"
      // spectrum for the edge fits.  Higher-momentum runs use the
      // uncut hit-energy spectrum for the maximum-bin feature.
      const char* histogramNameFormat =
        (run.info.runNumber == kStoppedRunNumber)
        ? "hFirstHitEnergyVsFiberID_Layer%d"
        : "hHitEnergyVsFiberID_Layer%d";

      std::unique_ptr<TH1D> hExp(
        GetFiberSpectrum(*run.exp, globalLayerID,
                         fiberIndex, expName,
                         histogramNameFormat));
      std::unique_ptr<TH1D> hSim(
        GetFiberSpectrum(*run.sim, globalLayerID,
                         fiberIndex, simName,
                         histogramNameFormat));

      if(!hExp || !hSim)
        continue;
      if(hExp->GetEntries() == 0 || hSim->GetEntries() == 0)
        continue;

      if(run.info.runNumber == kStoppedRunNumber){
        const EdgeFitResult leftExp =
          FitLeftEdge(*hExp, leftExpInitial);
        const EdgeFitResult leftSim =
          FitLeftEdge(*hSim, leftSimInitial);

        if(leftExp.valid && leftSim.valid){
          graph->SetPoint(graph->GetN(),
                          leftSim.position,
                          leftExp.position);

          // Seed the next fiber from this successful fit.
          leftExpInitial = leftExp.left;
          leftSimInitial = leftSim.left;
        }

        const EdgeFitResult rightExp =
          FitRightEdge(*hExp, rightExpInitial);
        const EdgeFitResult rightSim =
          FitRightEdge(*hSim, rightSimInitial);

        if(rightExp.valid && rightSim.valid){
          graph->SetPoint(graph->GetN(),
                          rightSim.position,
                          rightExp.position);

          rightExpInitial = rightExp.right;
          rightSimInitial = rightSim.right;
        }
      }
      else{
        const double expMPV = FindMPV(*hExp);
        const double simMPV = FindMPV(*hSim);

        graph->SetPoint(graph->GetN(), simMPV, expMPV);
      }
    }

    output.cd();
    graph->Write();
  }

  output.Close();
}
