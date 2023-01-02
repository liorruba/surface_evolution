// Class for layer object
class Layer {
public:
  double thickness;
  double regolithFraction;
  double iceFraction;
  double sootFraction;

  // First constructor for regular layer:
  Layer(double thickness, double regolithFraction, double iceFraction, double sootFraction);
  // Second constructor for a dummy layer used in printing the column. This layer is not normalized:
  Layer(double _thickness, double _regolithFraction);
  bool compareComposition(Layer layer);
  bool isEmpty();
  void consolidate(Layer layer);
  void shrink(double depthToRemove);
  void print(bool isNiceInterface = true);

private:
  void normalizeComposition();
};
