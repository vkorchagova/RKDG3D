#include "indicator_michalak.hpp"
#include "dg_conf.hpp"

IndicatorMichalak::IndicatorMichalak
(
   Averager& _avgr, 
   ParFiniteElementSpace* _fes,
   ParFiniteElementSpace* _fes_const, 
   const Array<int>& _offsets, 
   int _d
) : IndicatorBJ(_avgr, _fes, _fes_const, _offsets, _d) 
{};

void IndicatorMichalak::computeCorrectionFunction(const double& y, double& yMin)
{
   double yStar = DEFAULT_MICHALAK_YSTAR;
   double yRel = y / yStar;
   double yCur = y < 1.0 ? y + (3.0 - 2.0 * yStar) * yRel * yRel + (yStar - 2.0) * yRel * yRel * yRel : 1.0;
   yMin = yCur < yMin ? yCur : yMin;
}