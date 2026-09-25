#include "SimpleAFTTrackFinder.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

std::vector<AFTTrackResult>
SimpleAFTTrackFinder::FindTracks(
    const AFTEvent& event,
    const AFTGeometry& geometry) const
{
  auto tracksXZ =
    FindProjection(event,
                   geometry,
                   AFTProjection::XZ);

  auto tracksYZ =
    FindProjection(event,
                   geometry,
                   AFTProjection::YZ);

  tracksXZ.insert(tracksXZ.end(),
                  tracksYZ.begin(),
                  tracksYZ.end());

  return tracksXZ;
}

std::vector<SimpleAFTTrackFinder::Layer>
SimpleAFTTrackFinder::Layers(
    AFTProjection projection) const
{
  std::vector<Layer> layers;
  layers.reserve(18);

  for(int layerIndex = 0;
      layerIndex < 9;
      ++layerIndex){

    if(projection == AFTProjection::XZ){
      layers.push_back(
        {AFTLayerType::X, layerIndex});

      layers.push_back(
        {AFTLayerType::XPrime, layerIndex});
    }
    else{
      layers.push_back(
        {AFTLayerType::Y, layerIndex});

      layers.push_back(
        {AFTLayerType::YPrime, layerIndex});
    }
  }

  return layers;
}

std::vector<SimpleAFTTrackFinder::Candidate>
SimpleAFTTrackFinder::MakeSeeds(
    const AFTEvent& event,
    const AFTGeometry& geometry,
    AFTProjection projection) const
{
  const auto layers =
    Layers(projection);

  std::vector<Candidate> seeds;

  const int lastSeedLayer =
    std::min(parameters_.maxSeedLayer,
             static_cast<int>(layers.size()) - 1);

  for(int layerSequenceIndex = 0;
      layerSequenceIndex <= lastSeedLayer;
      ++layerSequenceIndex){

    const auto& layer =
      layers.at(layerSequenceIndex);

    bool foundHit = false;

    const std::size_t nFibers =
      geometry.GetNFibers(layer.type,
                          layer.index);

    for(std::size_t fiberIndex = 0;
        fiberIndex < nFibers;
        ++fiberIndex){

      const auto& fiber =
        geometry.GetFiber(
          layer.type,
          layer.index,
          static_cast<int>(fiberIndex));

      const AFTHit* hit =
        event.FindHit(fiber.fiberID);

      if(hit == nullptr){
        continue;
      }

      foundHit = true;

      Candidate candidate;

      candidate.track.AddHit(*hit);
      candidate.startLayer =
        layerSequenceIndex;
      candidate.currentLayer =
        layerSequenceIndex;

      seeds.push_back(candidate);
    }

    if(foundHit){
      break;
    }
  }

  return seeds;
}

double
SimpleAFTTrackFinder::Position(
    const AFTHit& hit,
    const AFTGeometry& geometry,
    AFTProjection projection) const
{
  const auto& fiber =
    geometry.GetFiber(hit.GetFiberID());

  if(projection == AFTProjection::XZ){
    return fiber.center.X();
  }

  return fiber.center.Y();
}

double
SimpleAFTTrackFinder::Z(
    const AFTHit& hit,
    const AFTGeometry& geometry) const
{
  return geometry.GetFiber(
    hit.GetFiberID()).center.Z();
}

double
SimpleAFTTrackFinder::Predict(
    const Candidate& candidate,
    double searchZ,
    const AFTGeometry& geometry,
    AFTProjection projection) const
{
  const auto& hits =
    candidate.track.GetHits();

  if(hits.empty()){
    throw std::runtime_error(
      "SimpleAFTTrackFinder::Predict: empty track");
  }

  const AFTHit& lastHit =
    hits.back();

  const double lastPosition =
    Position(lastHit,
             geometry,
             projection);

  const double lastZ =
    Z(lastHit, geometry);

  if(hits.size() == 1){
    return lastPosition;
  }

  const AFTHit& previousHit =
    hits[hits.size() - 2];

  const double previousPosition =
    Position(previousHit,
             geometry,
             projection);

  const double previousZ =
    Z(previousHit, geometry);

  const double dz =
    lastZ - previousZ;

  if(std::abs(dz) <
     std::numeric_limits<double>::epsilon()){
    return lastPosition;
  }

  const double slope =
    (lastPosition - previousPosition) / dz;

  return lastPosition +
         slope * (searchZ - lastZ);
}

double
SimpleAFTTrackFinder::Window(
    const Candidate& candidate,
    double searchZ,
    const AFTGeometry& geometry) const
{
  const auto& hits =
    candidate.track.GetHits();

  if(hits.empty()){
    throw std::runtime_error(
      "SimpleAFTTrackFinder::Window: empty track");
  }

  const double lastZ =
    Z(hits.back(), geometry);

  const double dz =
    std::abs(searchZ - lastZ);

  const double angle =
    (hits.size() == 1)
      ? parameters_.seedAngleTolerance
      : parameters_.angleTolerance;

  return parameters_.baseWindow +
         dz * std::tan(angle);
}

std::vector<const AFTHit*>
SimpleAFTTrackFinder::CompatibleHits(
    int searchLayer,
    const AFTEvent& event,
    const AFTGeometry& geometry,
    AFTProjection projection,
    double predictedPosition,
    double window) const
{
  const auto layers =
    Layers(projection);

  const auto& layer =
    layers.at(searchLayer);

  std::vector<const AFTHit*> compatibleHits;

  const std::size_t nFibers =
    geometry.GetNFibers(layer.type,
                        layer.index);

  for(std::size_t fiberIndex = 0;
      fiberIndex < nFibers;
      ++fiberIndex){

    const auto& fiber =
      geometry.GetFiber(
        layer.type,
        layer.index,
        static_cast<int>(fiberIndex));

    const AFTHit* hit =
      event.FindHit(fiber.fiberID);

    if(hit == nullptr){
      continue;
    }

    const double position =
      Position(*hit,
               geometry,
               projection);

    if(std::abs(position -
                predictedPosition) <= window){
      compatibleHits.push_back(hit);
    }
  }

  return compatibleHits;
}

std::vector<SimpleAFTTrackFinder::Candidate>
SimpleAFTTrackFinder::Extend(
    const Candidate& seed,
    const AFTEvent& event,
    const AFTGeometry& geometry,
    AFTProjection projection) const
{
  const auto layers =
    Layers(projection);

  std::vector<Candidate> active{seed};
  std::vector<Candidate> completed;

  while(!active.empty()){
    std::vector<Candidate> nextGeneration;

    for(const auto& candidate : active){
      const int searchLayer =
        candidate.currentLayer + 1;

      if(searchLayer >=
         static_cast<int>(layers.size())){
        completed.push_back(candidate);
        continue;
      }

      const auto& layer =
        layers.at(searchLayer);

      const auto& firstFiber =
        geometry.GetFiber(layer.type,
                          layer.index,
                          0);

      const double searchZ =
        firstFiber.center.Z();

      const double predictedPosition =
        Predict(candidate,
                searchZ,
                geometry,
                projection);

      const double window =
        Window(candidate,
               searchZ,
               geometry);

      const auto compatibleHits =
        CompatibleHits(searchLayer,
                       event,
                       geometry,
                       projection,
                       predictedPosition,
                       window);

      if(!compatibleHits.empty()){
        for(const AFTHit* hit :
            compatibleHits){

          Candidate branch =
            candidate;

          branch.track.AddHit(*hit);

          branch.currentLayer =
            searchLayer;

          branch.consecutiveMissing = 0;

          const double residual =
            Position(*hit,
                     geometry,
                     projection)
            - predictedPosition;

          branch.residual2 +=
            residual * residual;

          nextGeneration.push_back(
            branch);
        }
      }
      else{
        Candidate missed =
          candidate;

        missed.currentLayer =
          searchLayer;

        ++missed.consecutiveMissing;
        ++missed.totalMissing;

        if(missed.consecutiveMissing <=
           parameters_.maxMissingLayers){

          nextGeneration.push_back(
            missed);
        }
        else{
          completed.push_back(
            candidate);
        }
      }
    }

    active.swap(nextGeneration);
  }

  return completed;
}

AFTTrackResult
SimpleAFTTrackFinder::Result(
    const Candidate& candidate,
    AFTProjection projection) const
{
  AFTTrackResult result;

  result.track =
    candidate.track;

  result.projection =
    projection;

  result.startLayer =
    candidate.startLayer;

  result.endLayer =
    candidate.currentLayer;

  result.nMissing =
    candidate.totalMissing;

  result.geometryResidual =
    candidate.residual2;

  return result;
}

std::vector<AFTTrackResult>
SimpleAFTTrackFinder::FindProjection(
    const AFTEvent& event,
    const AFTGeometry& geometry,
    AFTProjection projection) const
{
  std::vector<AFTTrackResult> results;

  const auto seeds =
    MakeSeeds(event,
              geometry,
              projection);

  for(const auto& seed : seeds){
    const auto completed =
      Extend(seed,
             event,
             geometry,
             projection);

    for(const auto& candidate :
        completed){

      results.push_back(
        Result(candidate,
               projection));
    }
  }

  return results;
}
