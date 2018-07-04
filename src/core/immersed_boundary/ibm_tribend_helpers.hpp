
#ifndef IBM_TRIBEND_HELPERS_H
#define IBM_TRIBEND_HELPERS_H

#ifdef IMMERSED_BOUNDARY

// ******* Geometry stuff ********

// Structs

typedef struct {
  double el[3];
} Vector3D;

typedef struct {
  double el[3][3];
} Matrix3D;

inline void get_mi_vector_const(double res[3], const double a[3], const double b[3])
{
  int i;
  
  for(i=0;i<3;i++) {
    res[i] = a[i] - b[i];
#ifdef PARTIAL_PERIODIC
    if (PERIODIC(i))
#endif
      res[i] -= std::round(res[i]*box_l_i[i])*box_l[i];
  }
}

// Subtract: return a-b
inline void Subtr(Vector3D &res, const Vector3D &a, const Vector3D &b)
{
  get_mi_vector_const(res.el, a.el, b.el);
}

inline void Subtr(Vector3D &res, const Particle *const a, const Particle *const b)
{
  get_mi_vector_const(res.el, a->r.p, b->r.p);
}

// LengthSqr
inline double LengthSqr(const Vector3D &a)
{
  return a.el[0]*a.el[0] + a.el[1]*a.
  
  el[1] + a.el[2]*a.el[2];
}

// Length
inline double Length(const Vector3D &a)
{
  return sqrt(LengthSqr(a));
}

// ScalarProduct
inline double ScalarProduct(const Vector3D &a, const Vector3D &b)
{
  return a.el[0]*b.el[0] + a.el[1]*b.el[1] + a.el[2]*b.el[2];
}

// DyadicProduct

inline void DyadicProduct(Matrix3D &res, const Vector3D &a, const Vector3D &b)
{
  res.el[0][0] = a.el[0] * b.el[0];
  res.el[0][1] = a.el[0] * b.el[1];
  res.el[0][2] = a.el[0] * b.el[2];
  res.el[1][0] = a.el[1] * b.el[0];
  res.el[1][1] = a.el[1] * b.el[1];
  res.el[1][2] = a.el[1] * b.el[2];
  res.el[2][0] = a.el[2] * b.el[0];
  res.el[2][1] = a.el[2] * b.el[1];
  res.el[2][2] = a.el[2] * b.el[2];
}

// Add
inline void AddTo(Matrix3D &a, const Matrix3D &b)
{
  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      a.el[i][j] += b.el[i][j];
}

inline void AddScalarTo(Matrix3D &a, const double s)
{
  // Multiply scalar with unit matrix and add
  a.el[0][0] += s;
  a.el[1][1] += s;
  a.el[2][2] += s;
}

inline void AddTo(Vector3D &a, const Vector3D &b)
{
  a.el[0] += b.el[0];
  a.el[1] += b.el[1];
  a.el[2] += b.el[2];
}

// Multiply
inline void Multiply(Vector3D &res, const Vector3D &a, const double f)
{
  res.el[0] = a.el[0] * f;
  res.el[1] = a.el[1] * f;
  res.el[2] = a.el[2] * f;
}

// LeftVectorMatrix
inline void LeftVectorMatrix(Vector3D &res, const Vector3D &v, const Matrix3D &m)
{
  res.el[0] = v.el[0]*m.el[0][0] + v.el[1]*m.el[1][0] + v.el[2] * m.el[2][0];
  res.el[1] = v.el[0]*m.el[0][1] + v.el[1]*m.el[1][1] + v.el[2] * m.el[2][1];
  res.el[2] = v.el[0]*m.el[0][2] + v.el[1]*m.el[1][2] + v.el[2] * m.el[2][2];
}

// ******** More advanced routines *********

/***************
 CalcCosTheta
 ***********/
/// \brief Calculates the cosine of the angle between the vectors xi-xm and xj-xm.

inline double CalcCosTheta(const Particle *const xi, const Particle *const xj, const Particle *const xm)
{
  Vector3D dxim;
  Subtr(dxim, xi, xm);
  Vector3D dxjm;
  Subtr(dxjm, xj, xm);
  return ScalarProduct(dxim, dxjm) / (Length(dxim) * Length(dxjm));
}

/***************
 CalcCot
 ***********/

/// \brief Given the cosine of the angle, this function returns the cotangens. Note: This assumes that the angle is in [0;Pi].
inline double CalcCot(double cosTheta)
{
  return cosTheta / sqrt(1.0 - cosTheta*cosTheta);
}

/***************
 CalcCotDerivativeGompperAnalyt
 ***********/

/// \brief Calculates the gradient of the cotangens of the angle between xi-xm and xj-xm with respect to the node with ID \p derivativeNodeID.

void CalcCotDerivativeGompperAnalyt(Vector3D &cosThetaDeriv, const Particle *const xi, const Particle *const xj, const Particle *const xm, const int derivativeNodeID)
{
  // Goal: Calculate the gradient with respect to the node with ID "derivativeNodeID" of the cotangens of the angle between xi-xj and xi-xm.
  // Note: Compared with Maple, calculation is correct.
  
  const double cosTheta = CalcCosTheta(xi, xj, xm);
  
  // Get sine. Note: The angle is in [0;Pi], so this is possible.
  const double sinTheta = sqrt(1.0 - cosTheta*cosTheta);
  
  const double prefactor = 1.0 / (sinTheta*sinTheta*sinTheta);
  
  Vector3D dxim;
  Subtr(dxim, xi, xm);
  Vector3D dxjm;
  Subtr(dxjm, xj, xm);
  const double dximLen = Length(dxim);
  const double dxjmLen = Length(dxjm);
  
  // Kronecker-Deltas. "l" stands for "derivativeNodeID".
  const int ilDelta = (derivativeNodeID == xi->p.identity) ? 1 : 0;
  const int mlDelta = (derivativeNodeID == xm->p.identity) ? 1 : 0;
  const int jlDelta = (derivativeNodeID == xj->p.identity) ? 1 : 0;
  
  // Gradient of cosine theta.
  for (int i=0; i < 3; i++)
  {
    cosThetaDeriv.el[i]=	1.0/(dximLen*dxjmLen) * ((ilDelta - mlDelta) * dxjm.el[i] +
                                                   (jlDelta - mlDelta) * dxim.el[i] -
                                                   dxjmLen/dximLen * cosTheta * (ilDelta - mlDelta) * dxim.el[i] -
                                                   dximLen/dxjmLen * cosTheta * (jlDelta - mlDelta) * dxjm.el[i] );
  }
  cosThetaDeriv.el[0] *= prefactor;
  cosThetaDeriv.el[1] *= prefactor;
  cosThetaDeriv.el[2] *= prefactor;
}






#endif

#endif
