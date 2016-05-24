#include "stdafx.h"
#include "FEPressureLoad.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_PARAMETER_LIST(FEPressureLoad, FESurfaceLoad)
	ADD_PARAMETER(m_blinear , FE_PARAM_BOOL  , "linear"  );
	ADD_PARAMETER(m_pressure, FE_PARAM_DOUBLE, "pressure");
	ADD_PARAMETER(m_bsymm   , FE_PARAM_BOOL  , "symmetric_stiffness");
	ADD_PARAMETER(m_PC      , FE_PARAM_SURFACE_MAP, "pressure_map");
END_PARAMETER_LIST()

//-----------------------------------------------------------------------------
//! constructor
FEPressureLoad::FEPressureLoad(FEModel* pfem) : FESurfaceLoad(pfem), m_PC(FE_DOUBLE)
{ 
	m_blinear = false;
	m_pressure = 1.0;
	m_bsymm = true;

	// get the degrees of freedom
	m_dofX = pfem->GetDOFIndex("x");
	m_dofY = pfem->GetDOFIndex("y");
	m_dofZ = pfem->GetDOFIndex("z");
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEPressureLoad::SetSurface(FESurface* ps)
{
	FESurfaceLoad::SetSurface(ps);
	m_PC.Create(ps); 
	m_PC.SetValue(1.0);
}

//-----------------------------------------------------------------------------
//! \deprecated This function is only used by the 1.2 file reader and is to be 
//! considered obsolete.
bool FEPressureLoad::SetAttribute(const char* szatt, const char* szval)
{
	if (strcmp(szatt, "type") == 0)
	{
		if      (strcmp(szval, "linear"   ) == 0) SetLinear(true );
		else if (strcmp(szval, "nonlinear") == 0) SetLinear(false);
		else return false;
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPressureLoad::SetFacetAttribute(int nface, const char* szatt, const char* szval)
{
	if      (strcmp(szatt, "id") == 0) {}
//	else if (strcmp(szatt, "lc") == 0) pc.lc = atoi(szval) - 1;
	else if (strcmp(szatt, "scale") == 0)
	{
		double s = atof(szval);
		m_PC.SetValue(nface, s);
	}
	else return false;

	return true;
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure
void FEPressureLoad::PressureStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn)
{
	// choose the symmetric of unsymmetric formulation
	if (m_bsymm) 
		SymmetricPressureStiffness(el, ke, tn);
	else
		UnsymmetricPressureStiffness(el, ke, tn);
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure
void FEPressureLoad::SymmetricPressureStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn)
{
	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// gauss weights
	double* w = el.GaussWeights();

	// nodal coordinates
	FEMesh& mesh = *m_psurf->GetMesh();
	vec3d rt[FEElement::MAX_NODES];
	for (int j=0; j<neln; ++j) rt[j] = mesh.Node(el.m_node[j]).m_rt;

	// repeat over integration points
	ke.zero();
	for (int n=0; n<nint; ++n)
	{
		double* N = el.H(n);
		double* Gr = el.Gr(n);
		double* Gs = el.Gs(n);

		// traction at integration point
		double tr = 0;
		vec3d dxr(0,0,0), dxs(0,0,0);
		for (int i=0; i<neln; ++i) 
		{
			tr += N[i]*tn[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}
		
		// calculate stiffness component
		for (int i=0; i<neln; ++i)
			for (int j=0; j<neln; ++j)
			{
				vec3d kab = (dxr*(N[j]*Gs[i]-N[i]*Gs[j])
					   -dxs*(N[j]*Gr[i]-N[i]*Gr[j]))*w[n]*0.5*tr;

				ke.add(3*i, 3*j, mat3da(kab));
			}
	}
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure

void FEPressureLoad::UnsymmetricPressureStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn)
{
	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// gauss weights
	double* w = el.GaussWeights();

	// nodal coordinates
	FEMesh& mesh = *m_psurf->GetMesh();
	vec3d rt[FEElement::MAX_NODES];
	for (int j=0; j<neln; ++j) rt[j] = mesh.Node(el.m_node[j]).m_rt;

	// repeat over integration points
	ke.zero();
	for (int n=0; n<nint; ++n)
	{
		double* N = el.H(n);
		double* Gr = el.Gr(n);
		double* Gs = el.Gs(n);

		// traction at integration point
		double tr = 0;
		vec3d dxr(0,0,0), dxs(0,0,0);
		for (int i=0; i<neln; ++i) 
		{
			tr += N[i]*tn[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}
		
		// calculate stiffness component
		for (int i=0; i<neln; ++i)
			for (int j=0; j<neln; ++j)
			{
				vec3d Kab = (dxr*Gs[j] - dxs*Gr[j])*(tr*N[i]*w[n]);
				ke.sub(3*i, 3*j, mat3da(Kab));
			}
	}
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal forces due to hydrostatic pressure

void FEPressureLoad::PressureForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn)
{
	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// nodal coordinates
	FEMesh& mesh = *m_psurf->GetMesh();
	vec3d rt[FEElement::MAX_NODES];
	for (int j=0; j<neln; ++j) rt[j] = mesh.Node(el.m_node[j]).m_rt;

	// repeat over integration points
	zero(fe);
	double* w  = el.GaussWeights();
	for (int n=0; n<nint; ++n)
	{
		double* N  = el.H(n);
		double* Gr = el.Gr(n);
		double* Gs = el.Gs(n);

		// traction at integration points
		double tr = 0;
		vec3d dxr(0,0,0), dxs(0,0,0);
		for (int i=0; i<neln; ++i) 
		{
			tr += N[i]*tn[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}

		// force vector
		vec3d f = (dxr ^ dxs)*tr*w[n];

		for (int i=0; i<neln; ++i)
		{
			fe[3*i  ] += N[i]*f.x;
			fe[3*i+1] += N[i]*f.y;
			fe[3*i+2] += N[i]*f.z;
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal forces due to hydrostatic pressure

void FEPressureLoad::LinearPressureForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn)
{
	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// nodal coordinates
	FEMesh& mesh = *m_psurf->GetMesh();
	vec3d r0[FEElement::MAX_NODES];
	for (int i=0; i<neln; ++i) r0[i] = mesh.Node(el.m_node[i]).m_r0;

	// force vector
	vec3d f;

	// repeat over integration points
	zero(fe);
	double* w  = el.GaussWeights();
	for (int n=0; n<nint; ++n)
	{
		double* N  = el.H(n);
		double* Gr = el.Gr(n);
		double* Gs = el.Gs(n);

		// traction at integration points
		double tr = 0;
		vec3d dxr(0,0,0), dxs(0,0,0);
		for (int i=0; i<neln; ++i) 
		{
			tr += N[i]*tn[i];
			dxr += r0[i]*Gr[i];
			dxs += r0[i]*Gs[i];
		}

		f = (dxr ^ dxs)*tr*w[n];

		for (int i=0; i<neln; ++i)
		{
			fe[3*i  ] += N[i]*f.x;
			fe[3*i+1] += N[i]*f.y;
			fe[3*i+2] += N[i]*f.z;
		}
	}
}

//-----------------------------------------------------------------------------

void FEPressureLoad::Serialize(DumpStream& ar)
{
	FESurfaceLoad::Serialize(ar);

	m_PC.Serialize(ar);
}

//-----------------------------------------------------------------------------
void FEPressureLoad::UnpackLM(FEElement& el, vector<int>& lm)
{
	FEMesh& mesh = *GetSurface().GetMesh();
	int N = el.Nodes();
	lm.resize(N*3);
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = mesh.Node(n);
		vector<int>& id = node.m_ID;

		lm[3*i  ] = id[m_dofX];
		lm[3*i+1] = id[m_dofY];
		lm[3*i+2] = id[m_dofZ];
	}
}

//-----------------------------------------------------------------------------
void FEPressureLoad::StiffnessMatrix(const FETimePoint& tp, FESolver* psolver)
{
	// We only need the stiffness for nonlinear pressure forces
	if (m_blinear) return;

	matrix ke;
	vector<int> lm;

	FESurface& surf = GetSurface();
	int npr = surf.Elements();
	for (int m=0; m<npr; ++m)
	{
		// get the surface element
		FESurfaceElement& el = m_psurf->Element(m);

		// calculate nodal normal tractions
		int neln = el.Nodes();
		vector<double> tn(neln);

		// evaluate the prescribed traction.
		// note the negative sign. This is because this boundary condition uses the 
		// convention that a positive pressure is compressive
		for (int j=0; j<neln; ++j) tn[j] = -m_pressure*m_PC.GetValue<double>(m);

		// get the element stiffness matrix
		int ndof = 3*neln;
		ke.resize(ndof, ndof);

		// calculate pressure stiffness
		PressureStiffness(el, ke, tn);

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble element matrix in global stiffness matrix
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
void FEPressureLoad::Residual(const FETimePoint& tp, FEGlobalVector& R)
{
	vector<double> fe;
	vector<int> lm;

	FESurface& surf = GetSurface();
	int npr = surf.Elements();
	for (int i=0; i<npr; ++i)
	{
		FESurfaceElement& el = m_psurf->Element(i);

		// calculate nodal normal tractions
		int neln = el.Nodes();
		vector<double> tn(neln);

		// evaluate the prescribed traction.
		// note the negative sign. This is because this boundary condition uses the 
		// convention that a positive pressure is compressive
		for (int j=0; j<el.Nodes(); ++j) tn[j] = -m_pressure*m_PC.GetValue<double>(i);
		
		int ndof = 3*neln;
		fe.resize(ndof);

		if (m_blinear) LinearPressureForce(el, fe, tn); else PressureForce(el, fe, tn);

		// get the element's LM vector
		UnpackLM(el, lm);

		// add element force vector to global force vector
		R.Assemble(el.m_node, lm, fe);
	}
}
