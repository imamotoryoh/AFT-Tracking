#include <iostream>
#include <memory>
#include <string>
#include <chrono>

#include <TApplication.h>
#include <TSystem.h>

#include "AFTEvent.hh"
#include "AFTEventDisplay.hh"
#include "AFTEventReader.hh"
#include "AFTGeometry.hh"
#include "AFTTrackEvaluator.hh"
#include "SimpleAFTTrackFinder.hh"
#include "ProtonCalibrationTrackSelector.hh"
#include "ProtonCalibrationAnalyzer.hh"

int
main(int argc, char** argv)
{
  if(argc < 4){
    std::cerr
      << "Usage: "
      << argv[0]
      << " geometry.dat input.root treeName [--display]"
      << std::endl;

    return 1;
  }

  const std::string geometryFile = argv[1];
  const std::string inputFile = argv[2];
  const std::string treeName = argv[3];

  bool displayMode = false;

  for(int i = 4; i < argc; ++i){
    const std::string option = argv[i];

    if(option == "--display"){
      displayMode = true;
    }
    else{
      std::cerr
        << "Unknown option: "
        << option
        << std::endl;

      return 1;
    }
  }


  AFTGeometry geometry;

  if(!geometry.Load(geometryFile)){
    std::cerr
      << "Failed to load geometry file: "
      << geometryFile
      << std::endl;

    return 1;
  }


  AFTEventReader reader(geometry);

  if(!reader.Open(inputFile, treeName)){
    std::cerr
      << "Failed to open input file: "
      << inputFile
      << std::endl;

    return 1;
  }


  SimpleAFTTrackFinder trackFinder;
  AFTTrackEvaluator evaluator(geometry);
  ProtonCalibrationTrackSelector selector;
  double calib = 1;
  if(treeName == "ea0c") calib = 0.4738;
  ProtonCalibrationAnalyzer analyzer(geometry,calib);

  std::unique_ptr<TApplication> app;
  std::unique_ptr<AFTEventDisplay> display;

  if(displayMode){
    int rootArgc = 1;
    char* rootArgv[] = {
      argv[0],
      nullptr
    };

    app =
      std::make_unique<TApplication>(
        "AFTTracking",
        &rootArgc,
        rootArgv
      );

    display =
      std::make_unique<AFTEventDisplay>(
        geometry
      );
  }


  const Long64_t nEntries =
    reader.GetEntries();
  //  const Long64_t nEntries = 1000000;

  std::cout
    << "Number of entries = "
    << nEntries
    << std::endl;

  if(displayMode){
    std::cout
      << "Mode: Event Display"
      << std::endl;
  }
  else{
    std::cout
      << "Mode: Analysis"
      << std::endl;
  }

  double DurationTime = 0.0;
  // Event loop  
  for(Long64_t entry = 0; entry < nEntries; ++entry){
    //    const auto start = std::chrono::steady_clock::now();

    if(entry % 100000 == 0){
      std::cout
	<< "Processing entry "
	<< entry
	<< " / "
	<< nEntries
	<< std::endl;
    }


    AFTEvent event;
    
    if(!reader.GetEvent(entry, event)){
      continue;
    }

    // const auto end = std::chrono::steady_clock::now();
    // DurationTime = std::chrono::duration<double>(end - start).count();
    
    if(event.GetNHits() == 0){
      continue;
    }    
    const auto tracks =
      trackFinder.FindTracks(
        event,
        geometry
      );
   
    evaluator.Process(
      event,
      tracks
    );
    
    const auto selection =
      selector.Select(tracks, geometry);
    
    analyzer.Process(selection);    
    
    if(displayMode){
      display->DrawEvent(
        event,
        tracks
      );

      gSystem->ProcessEvents();
      
      std::cout
        << "Entry "
        << entry
        << " : "
        << tracks.size()
        << " candidates"
        << std::endl;

      std::cout
        << "Press Enter for next event, "
        << "or q to quit: "
        << std::flush;

      std::string input;
      std::getline(std::cin, input);

      if(input == "q" ||
         input == "Q"){
        break;
      }
    }
    
  }


  evaluator.PrintSummary();

  const std::string outputFileEva =
    "root/AFTTrackEvaluation.root";
  const std::string outputFileCal =
    "root/AFTTrackCalibration" + treeName + ".root";

  if(!evaluator.Write(outputFileEva)){
    std::cerr
      << "Failed to write output ROOT file: "
      << outputFileEva
      << std::endl;
    
    return 1;
  }

  if(!analyzer.Write(outputFileCal)){
    std::cerr
      << "Failed to write output ROOT file: "
      << outputFileCal
      << std::endl;
    
    return 1;
  }

  
  std::cout
    << "Output written to "
    << outputFileEva
    << "and "
    << outputFileCal
    << std::endl;

  reader.Close();

  return 0;
}
