#include "ProtonCalibrationAnalyzer.hh"

#include <TFile.h>

ProtonCalibrationAnalyzer::
ProtonCalibrationAnalyzer(
			  const AFTGeometry& geometry,
			  const double calib)
  : geometry_(geometry), calib_(calib)
{
  hFirstHitEnergy_ =
    new TH1D(
      "hFirstHitEnergy",
      "First hit from track end;Energy;Counts",
      400,
      0.0,
      20.0
    );

  hSecondHitEnergy_ =
    new TH1D(
      "hSecondHitEnergy",
      "Second hit from track end;Energy;Counts",
      400,
      0.0,
      20.0
    );

  hThirdHitEnergy_ =
    new TH1D(
      "hThirdHitEnergy",
      "Third hit from track end;Energy;Counts",
      400,
      0.0,
      20.0
    );

  hFirstHitEnergyVsFiberID_ =
    new TH2D(
      "hFirstHitEnergyVsFiberID",
      "First hit from track end;"
      "Fiber ID;Energy",
      864,
      -0.5,
      863.5,
      400,
      0.0,
      20.0
    );


  hSecondHitEnergyVsFiberID_ =
    new TH2D(
      "hSecondHitEnergyVsFiberID",
      "Second hit from track end;"
      "Fiber ID;Energy",
      864,
      -0.5,
      863.5,
      400,
      0.0,
      20.0
    );

  hThirdHitEnergyVsFiberID_ =
    new TH2D(
      "hThirdHitEnergyVsFiberID",
      "Third hit from track end;"
      "Fiber ID;Energy",
      864,
      -0.5,
      863.5,
      400,
      0.0,
      20.0
    );

  for(int layer = 0;
      layer < NGlobalLayers;
      ++layer){
    
    hFirstHitEnergyVsFiberIDForEachLayer_[layer] =
      new TH2D(
	       Form(
		    "hFirstHitEnergyVsFiberID_Layer%02d",
		    layer
		    ),
	       Form(
		    "Last hit - Global Layer %d;"
		    "Fiber ID;Energy",
		    layer
		    ),
	       864,
	       -0.5,
	       863.5,
	       400,
	       0.0,
	       20.0
	       );
    
    hSecondHitEnergyVsFiberIDForEachLayer_[layer] =
      new TH2D(
	       Form(
		    "hSecondHitEnergyVsFiberID_Layer%02d",
		    layer
		    ),
	       Form(
		    "Second hit - Global Layer %d;"
		    "Fiber ID;Energy",
		    layer
		    ),
	       864,
	       -0.5,
	       863.5,
	       400,
	       0.0,
	       20.0
	       );

    hThirdHitEnergyVsFiberIDForEachLayer_[layer] =
      new TH2D(
	       Form(
		    "hThirdHitEnergyVsFiberID_Layer%02d",
		    layer
		    ),
	       Form(
		    "Third hit - Global Layer %d;"
		    "Fiber ID;Energy",
		    layer
		    ),
	       864,
	       -0.5,
	       863.5,
	       400,
	       0.0,
	       20.0
	       );
    
    hFirstHitEnergyVsFiberIDForEachLayer_[layer]->
      SetDirectory(nullptr);
    
    hSecondHitEnergyVsFiberIDForEachLayer_[layer]->
      SetDirectory(nullptr);
    
    hThirdHitEnergyVsFiberIDForEachLayer_[layer]->
      SetDirectory(nullptr);
  }
  
  for(auto* hist :
      {hFirstHitEnergy_,
       hSecondHitEnergy_,
       hThirdHitEnergy_}){
    hist->SetDirectory(nullptr);
  }

  for(auto* hist :
      {hFirstHitEnergyVsFiberID_,
       hSecondHitEnergyVsFiberID_,
       hThirdHitEnergyVsFiberID_}){
    hist->SetDirectory(nullptr);
  }
}

ProtonCalibrationAnalyzer::
~ProtonCalibrationAnalyzer()
{
  delete hFirstHitEnergy_;
  delete hSecondHitEnergy_;
  delete hThirdHitEnergy_;

  delete hFirstHitEnergyVsFiberID_;
  delete hSecondHitEnergyVsFiberID_;
  delete hThirdHitEnergyVsFiberID_;

  for(int layer = 0;
      layer < NGlobalLayers;
      ++layer){
    
   delete hFirstHitEnergyVsFiberIDForEachLayer_[layer];
    
   delete hSecondHitEnergyVsFiberIDForEachLayer_[layer];
    
   delete hThirdHitEnergyVsFiberIDForEachLayer_[layer];
  }

}

void ProtonCalibrationAnalyzer::Process(
    const AFTTrackSelectionResult& selection)
{
  if(!selection.IsSelected()){
    return;
  }

  ProcessTrack(*selection.trackXZ);
  ProcessTrack(*selection.trackYZ);
}

void ProtonCalibrationAnalyzer::ProcessTrack(
    const AFTTrackResult& result)
{
  const auto& track =
    result.track;

  if(track.GetNHits() >= 1){
    FillHit(
      track.GetHitFromEnd(0),
      0
    );
  }

  if(track.GetNHits() >= 2){
    FillHit(
      track.GetHitFromEnd(1),
      1
    );
  }

  if(track.GetNHits() >= 3){
    FillHit(
      track.GetHitFromEnd(2),
      2
    );
  }
}

void ProtonCalibrationAnalyzer::FillHit(
    const AFTHit& hit,
    std::size_t indexFromEnd)
{
  const double energy =
    hit.GetEnergy();
  
  const int fiberID =
    hit.GetFiberID();

  const auto& fiber =
    geometry_.GetFiber(fiberID);
  
  const int globalLayerID =
    fiber.globalLayerID;
  
  if(indexFromEnd == 0){
    hFirstHitEnergy_->Fill(
      energy*calib_
    );

    hFirstHitEnergyVsFiberID_->Fill(
      fiberID,
      energy*calib_
    );

    hFirstHitEnergyVsFiberIDForEachLayer_
      [globalLayerID]->Fill(
			    fiberID,
			    energy*calib_
			    );
  }
  else if(indexFromEnd == 1){
    hSecondHitEnergy_->Fill(
      energy*calib_
    );

    hSecondHitEnergyVsFiberID_->Fill(
      fiberID,
      energy*calib_
    );

    hSecondHitEnergyVsFiberIDForEachLayer_
      [globalLayerID]->Fill(
			    fiberID,
			    energy*calib_
			    );

  }
  else if(indexFromEnd == 2){
    hThirdHitEnergy_->Fill(
      energy*calib_
    );

    hThirdHitEnergyVsFiberID_->Fill(
      fiberID,
      energy*calib_
    );

    hThirdHitEnergyVsFiberIDForEachLayer_
      [globalLayerID]->Fill(
			    fiberID,
			    energy*calib_
			    );

  }
}

bool ProtonCalibrationAnalyzer::Write(
    const std::string& filename) const
{
  TFile output(
    filename.c_str(),
    "RECREATE"
  );

  if(output.IsZombie()){
    return false;
  }

  hFirstHitEnergy_->Write();
  hSecondHitEnergy_->Write();
  hThirdHitEnergy_->Write();

  hFirstHitEnergyVsFiberID_->Write();
  hSecondHitEnergyVsFiberID_->Write();
  hThirdHitEnergyVsFiberID_->Write();

  for(int layer = 0;
      layer < NGlobalLayers;
      ++layer){
    
    hFirstHitEnergyVsFiberIDForEachLayer_[layer]->
      Write();
    
    hSecondHitEnergyVsFiberIDForEachLayer_[layer]->
      Write();
    
    hThirdHitEnergyVsFiberIDForEachLayer_[layer]->
      Write();
  }
  
  output.Close();

  return true;
}
