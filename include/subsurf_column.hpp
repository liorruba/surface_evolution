// Class for subsurf column object
class SubsurfColumn {
public:
  std::vector<Layer> subsurfLayers;
  bool isPermShadow;

  SubsurfColumn();
  double getSurfaceElevation();
  void InitializeColumn();
  void addLayer(Layer layer);
  void removeMaterial(double depthToRemove);
  Layer integrateColumnComposition(double depthToIntergrate);
  void print(bool isNiceInterface = true);

private:
  void removeAllLayers();
  double surfaceElevation;
};
