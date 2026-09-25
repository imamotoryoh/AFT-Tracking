#ifndef AFT_TRACK_RESULT_HH
#define AFT_TRACK_RESULT_HH

#include "AFTTrack.hh"

enum class AFTProjection
{
  XZ,
  YZ
};

struct AFTTrackResult
{
  AFTTrack track;
  AFTProjection projection = AFTProjection::XZ;

  // Index in the 18-layer tracking sequence:
  // XZ: X0, X'0, X1, X'1, ...
  // YZ: Y0, Y'0, Y1, Y'1, ...
  int startLayer = -1;
  int endLayer = -1;

  int nMissing = 0;

  // Sum of squared transverse residuals [mm^2].
  double geometryResidual = 0.0;
};

#endif
