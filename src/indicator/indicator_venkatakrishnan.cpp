#include "indicator_venkatakrishnan.hpp"

IndicatorVenkatakrishnan::IndicatorVenkatakrishnan
(
   Averager& _avgr, 
   ParFiniteElementSpace* _fes,
   ParFiniteElementSpace* _fes_const, 
   const Array<int>& _offsets, 
   int _d
) : IndicatorBJ(_avgr, _fes, _fes_const, _offsets, _d) 
{};

void IndicatorVenkatakrishnan::computeCorrectionFunction(const double& y, double& yMin)
{
   double yCur = (y*y + 2.0 * y) / (y*y + y + 2.0);
   yMin = yCur < yMin ? yCur : yMin;
}