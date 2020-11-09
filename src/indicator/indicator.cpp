#include "indicator.hpp"


Indicator::Indicator (ParFiniteElementSpace* _fes, const Array<int>& _offsets, int _d, BlockVector& _idata) : 
   fes(_fes), 
   offsets(_offsets), 
   dim(_d),
   values(_idata)
{
   mesh = fes->GetParMesh();
   parGridX.MakeRef(_fes, NULL);
};

