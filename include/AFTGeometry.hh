#ifndef AFT_GEOMETRY_HH
#define AFT_GEOMETRY_HH

enum class AFTLayerType
  {
    X,
    XPrime,
    Y,
    YPrime
  };

struct FiberGeometry
{
  int fiberID = -1;
  int planeID = -1;//easiroc plane id
  int globalLayerID = -1;// 0~36 real layer

  AFTLayerType layerType; //0=X, 1=X', 2=Y, 3=Y'
  int layerIndex = -1; // 0-8 within X,X',Y,Y'
  int fiberIndex = -1; // 0-31 for X/X', 0-15 for Y/Y'
  
  TVector3 center;
  TVector3 direction;

  double length = 0.0;
  double radius = 0.0;
};

class AFTGeometry
{
public:

  static constexpr std::size_t NLayerTypes  = 4;
  static constexpr std::size_t NLayerPerType = 9;
  bool Load(const std::string& filename);

  //Acess by global fiber ID
  const FiberGeometry& GetFiber(int fiberID) const; 
  //Acess by detector structure
  const FiberGeometry& GetFiber(AFTLayerType type,
				int layerIndex,
				int fiberIndex) const; 

  std::size_t GetNFibers() const;
  std::size_t GetNFibers(AFTLayerType type,
			 int layerIndex) const;

  const std::vector<std::size_t>&
  GetFiberIndices(AFTLayerType type,
		  int layerIndex) const;
  
  TVector3 GetCenter(int fiberID) const;
  TVector3 GetDirection(int fiberID) const;

private:
  std::vector<FiberGeometry> fibers_;

  std::array<
    std::array<
    std::vector<std::size_t>,
    NLayerPerType
      >,
    NLayerTypes
    > fiberIndicesForEachlayerType_;

  static constexpr std::size_t
  LayerTypeToIndex(AFTLayerType type)
  {
    return static_cast<std::size_t>(type);
  }

};
#endif
