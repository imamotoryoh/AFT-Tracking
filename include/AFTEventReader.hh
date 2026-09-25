#ifndef AFT_EVENT_READER_HH
#define AFT_EVENT_READER_HH

#include <string>

#include <RtypesCore.h>

#include "AFTEvent.hh"
#include "AFTGeometry.hh"

class TFile;
class TTree;

class AFTEventReader
{
public:
  explicit AFTEventReader(
      const AFTGeometry& geometry);

  ~AFTEventReader();

  bool Open(const std::string& filename,
            const std::string& treeName);

  void Close();

  Long64_t GetEntries() const;

  bool GetEvent(Long64_t entry,
                AFTEvent& event);

private:
  const FiberGeometry*
  FindFiber(int globalLayerID,
            int fiberIndex) const;

  const AFTGeometry& geometry_;

  TFile* file_ = nullptr;
  TTree* tree_ = nullptr;
  std::array<std::array<const FiberGeometry*, 32>, 36> fiberLookup_;
  Double_t aft_de_low_cut_inside_[36][32] = {};
};

#endif
