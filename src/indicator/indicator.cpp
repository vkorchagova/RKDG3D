#include "indicator.hpp"


Indicator::Indicator
(
   Averager& _avgr, 
   ParFiniteElementSpace* _fes,
   ParFiniteElementSpace* _fes_const, 
   const Array<int>& _offsets, 
   int _d
) : 
   averager(_avgr),
   fes(_fes),
   fes_const(_fes_const), 
   offsets(_offsets), 
   dim(_d)
{
   mesh = fes->GetParMesh();
   values = new ParGridFunction(fes_const);
};

Indicator::~Indicator()
{
	delete values;
}

void Indicator::setValue(int iCell, double val)
{
	Vector el_ind(num_equation);

	fe = fes->GetFE(iCell);
	fes_const->GetElementVDofs(iCell, el_vdofs);

	values->GetSubVector(el_vdofs, el_x);

   el_x[0] = val;

   values->SetSubVector(el_vdofs, el_x);
}

const double Indicator::getValue(int iCell)
{
	Vector el_ind(num_equation);

	fe = fes->GetFE(iCell);
	fes_const->GetElementVDofs(iCell, el_vdofs);

	values->GetSubVector(el_vdofs, el_x);

   return el_x[0];
}