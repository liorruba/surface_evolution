// Class for subsurf column object
class SubsurfColumn {
public:
  std::vector<Layer> subsurfLayers;

  SubsurfColumn();
  double get_surfaceElevation();
  void InitializeColumn();
  void addLayer(Layer layer);
  void removeMaterial(double depthToRemove);
  Layer integrateColumnComposition(double depthToIntergrate);
  void print();

private:
  void removeAllLayers();
  double surfaceElevation;
};
