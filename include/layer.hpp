// Class for layer object
class Layer {
public:
  double thickness;
  double regolithFraction;
  double iceFraction;
  double sootFraction;

  Layer(double thickness = 0, double regolithFraction = 0, double iceFraction = 0, double sootFraction = 0);

  bool compareComposition(Layer layer);
  void consolidate(Layer layer);
  void shrink(double depthToRemove);
  void print();

private:
  void normalizeComposition();
};
