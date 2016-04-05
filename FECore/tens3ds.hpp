// NOTE: This file is automatically included from tens3drs.h
// Users should not include this file manually!

inline tens3ds::tens3ds()
{
	zero();
}

inline tens3ds::tens3ds(const double g)
{
	for (int i = 0; i < NNZ; i++)
		d[i] = g;
}

inline tens3ds::tens3ds(double m[10])
{
	for (int i = 0; i < NNZ; i++)
		d[i] = m[i];
}

// operator +
inline tens3ds tens3ds::operator + (const tens3ds& t) const
{
	tens3ds s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] + t.d[i];
	
	return s;
}

// operator -
inline tens3ds tens3ds::operator - (const tens3ds& t) const
{
	tens3ds s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] - t.d[i];

	return s;
}

// operator *
inline tens3ds tens3ds::operator * (double g) const
{
	tens3ds s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = g*d[i];
	
	return s;
}

// operator /
inline tens3ds tens3ds::operator / (double g) const
{
	tens3ds s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i]/g;
	
	return s;
}

// assignment operator +=
inline tens3ds& tens3ds::operator += (const tens3ds& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] += t.d[i];
	
	return (*this);
}

// assignment operator -=
inline tens3ds& tens3ds::operator -= (const tens3ds& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] -= t.d[i];
	
	return (*this);
}

// assignment operator *=
inline tens3ds& tens3ds::operator *= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] *= g;
	
	return (*this);
}

// assignment operator /=
inline tens3ds& tens3ds::operator /= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] /= g;
	
	return (*this);
}

// unary operator -
inline tens3ds tens3ds::operator - () const
{
	tens3ds s;
	s.d[0] = -d[0];
	s.d[1] = -d[1];
	s.d[2] = -d[2];
	s.d[3] = -d[3];
	s.d[4] = -d[4];
	s.d[5] = -d[5];
	s.d[6] = -d[6];
	s.d[7] = -d[7];
	s.d[8] = -d[8];
	s.d[9] = -d[9];
	return s;
}

// intialize to zero
inline void tens3ds::zero()
{
	d[0] = d[1] = d[2] = d[3] = d[4] = d[5] = d[6] = d[7] = d[8] = d[9] = 0.0;
}

// contract the right two legs by the dyad formed by a vector  xi = Tijk*Xi*Xk
inline vec3d tens3ds::contractdyad1(const vec3d& v)
{
    vec3d x;
	x.x = d[0]*v.x*v.x + 2*d[1]*v.x*v.y + 2*d[2]*v.x*v.z + d[3]*v.y*v.y + 2*d[4]*v.y*v.z + d[5]*v.z*v.z;
	x.y = d[1]*v.x*v.x + 2*d[3]*v.x*v.y + 2*d[4]*v.x*v.z + d[6]*v.y*v.y + 2*d[7]*v.y*v.z + d[8]*v.z*v.z;
	x.z = d[2]*v.x*v.x + 2*d[4]*v.x*v.y + 2*d[5]*v.x*v.z + d[7]*v.y*v.y + 2*d[8]*v.y*v.z + d[9]*v.z*v.z;

	return x;
}

// triple contraction by a similar 3o tensor m = Tijk*Hijk
inline double tens3ds::tripledot(const tens3ds& H)
{
	const double* h = H.d;
	return d[0]*h[0] + 3*d[1]*h[1] + 3*d[2]*h[2] + 3*d[3]*h[3] + 6*d[4]*h[4] + 3*d[5]*h[5] + d[6]*h[6] + 3*d[7]*h[7] + 3*d[8]*h[8] + d[9]*h[9];
}

// calculates the symmetric tensor A_ijk = (l_i*m_j*r_k + perm(i,j,k))/6
inline tens3ds dyad3s(const vec3d& l, const vec3d& m, const vec3d& r)
{
	tens3ds a;
	a.d[0] = (l.x*m.x*r.x); 
	a.d[1] = (l.x*m.x*r.y + l.x*m.y*r.x + l.y*m.x*r.x)/3.0; 
	a.d[2] = (l.x*m.x*r.z + l.x*m.z*r.x + l.z*m.x*r.x)/3.0;
	a.d[3] = (l.x*m.y*r.y + l.y*m.x*r.y + l.y*m.y*r.x)/3.0; 
	a.d[4] = (l.x*m.y*r.z + l.y*m.x*r.z + l.z*m.y*r.x + l.x*m.z*r.y + l.z*m.x*r.y + l.y*m.z*r.x)/6.0; 
	a.d[5] = (l.x*m.z*r.z + l.z*m.x*r.z + l.z*m.z*r.x)/3.0;
	a.d[6] = (l.y*m.y*r.y); 
	a.d[7] = (l.y*m.y*r.z + l.y*m.z*r.y + l.z*m.y*r.y)/3.0;
	a.d[8] = (l.y*m.z*r.z + l.z*m.y*r.z + l.z*m.z*r.y)/3.0;
	a.d[9] = (l.z*m.z*r.z);
	return a;
}
