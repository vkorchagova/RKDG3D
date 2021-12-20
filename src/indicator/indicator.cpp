#include "indicator.hpp"


Indicator::Indicator (Averager& _avgr, ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d, BlockVector& _idata) : 
   averager(_avgr),
   fes(_fes), 
   offsets(_offsets), 
   dim(_d),
   values(_idata)
{
   mesh = fes->GetParMesh();
   minValues.SetSize(mesh->GetNE());
};

