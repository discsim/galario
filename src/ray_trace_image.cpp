
#include "assert.h"

#include<cmath>
#include<limits>
#include<memory>
#include<iostream>
#include<fstream>
#include<vector>


template<class T>
inline T my_min(T x, T y) {
  return x < y ? x : y ;
}
template<class T>
inline T my_max(T x, T y) {
  return x > y ? x : y ;
}

template<class T>
class Vector3 {
public:
  Vector3() : _x{0,0,0}  { } ;
  Vector3(T x, T y, T z)
    : _x{x,y,z}
  { } ;
  
  T& operator[](int i) {
    return _x[i] ;
  }
  const T& operator[](int i) const {
    return _x[i] ;
  }
  
private:
  T _x[3] ;
} ;

using Vector3d = Vector3<double> ;
using Vector3f = Vector3<float> ;
using Vector3i = Vector3<int>;

struct Box {
  Vector3d min, max ;
} ;


struct Ray {

  Ray(const Vector3d& orig, const Vector3d&  direction)
    : origin(orig),
      dir(direction)
  {
    double norm =
      std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    for (int i=0; i < 3; i++) {
      dir[i] /= norm ;
      inv_dir[i] = 1 / dir[i] ;      
      side[i] = inv_dir[i] < 0 ; 
    }      
  }

  Vector3d pos(double l) const {
    return Vector3d(origin[0] + l*dir[0],
		    origin[1] + l*dir[1],
		    origin[2] + l*dir[2]) ;		    
  }
  
    
  
  Vector3d origin ;
  Vector3d dir ;
  Vector3d inv_dir ;
  Vector3i side ;
} ;

struct RayBoxIntersection {

  RayBoxIntersection(const Ray& ray, const Box& box) {
    Vector3d x0((box.min[0] - ray.origin[0])*ray.inv_dir[0],
		(box.min[1] - ray.origin[1])*ray.inv_dir[1],
		(box.min[2] - ray.origin[2])*ray.inv_dir[2]) ;
    
    Vector3d x1((box.max[0] - ray.origin[0])*ray.inv_dir[0],
		(box.max[1] - ray.origin[1])*ray.inv_dir[1],
		(box.max[2] - ray.origin[2])*ray.inv_dir[2]) ;

    // Nan-safe max/min
    l0 = x0[0] ; l1 = x1[0] ;
    for (int i=1; i < 3; i++) {
      l0 = my_max(l0, my_max(x0[i],-std::numeric_limits<double>::infinity()));
      l1 = my_min(l1, my_min(x1[i], std::numeric_limits<double>::infinity()));
    }

    assert(l1 == l1 && l0 == l0) ;
    
    dl = std::max(l1 - l0, 0.) ;

    for (int i=0; i < 3; i++) {
      entry[i] = (l0 == x0[i]) ;
      exit [i] = (l1 == x1[i]) ;

      if (ray.dir[i] < 0) 
	exit[i] *= -1 ;
      else
	entry[i] *= -1 ;
    }
  }

  double l0, l1, dl ;
  Vector3i entry, exit ;  
} ;


class Grid {
public:
  Grid(const Vector3i& Nx_,
       const Vector3d& xmin,
       const Vector3d& dx_={1,1,1})
    : Nx(Nx_),
      dx(dx_)
  {
    for (int i(0); i < 3; i++) {
      bounds.min[i] = xmin[i] ;
      bounds.max[i] = xmin[i] + Nx[i] * dx[i] ;
    }
  }

  Vector3i cell_index(const Vector3d& pos) const {
    return Vector3i( (pos[0] - bounds.min[0]) / dx[0],
		     (pos[1] - bounds.min[1]) / dx[1],
		     (pos[2] - bounds.min[2]) / dx[2]) ;
    
  }

  Vector3i cell_index_safe(const Vector3d& pos) const {
    
    Vector3i idx = cell_index(pos) ;
    
    for (int k=0; k < 3; k++) {
      if (dx[k] > 0) {
	if (pos[k] > (bounds.min[k] +  idx[k]*dx[k]))
	  idx[k]++ ;
	if (pos[k] < (bounds.min[k] +  (idx[k]+1)*dx[k]))
	  idx[k]--;
      }
      else {
	if (pos[k] > (bounds.min[k] +  idx[k]*dx[k]))
	  idx[k]-- ;
	if (pos[k] < (bounds.min[k] +  (idx[k]+1)*dx[k]))
	  idx[k]++ ;
      }
    }
    return idx ;
  }

  int cell_index(const Vector3i& idx) const {
    return idx[0]*Nx[1]*Nx[2] + idx[1]*Nx[1] + idx[2] ;
  }
  
  Box cell_box(const Vector3i idx) const {

    Box box ;
    for (int i(0); i < 3; i++) {
      box.min[i] = bounds.min[i] + dx[i] *  idx[i] ;
      box.max[i] = bounds.min[i] + dx[i] * (idx[i] + 1) ;
    }
    return box ;
  }

  const Box& domain() const {
    return bounds ;
  }

  int size(int i) const {
    return Nx[i] ;
  }

  bool in_domain(const Vector3i idx) const {
    return
      (idx[0] >= 0 && idx[0] < Nx[0]) &&
      (idx[1] >= 0 && idx[1] < Nx[1]) &&
      (idx[2] >= 0 && idx[2] < Nx[2]) ;
  }

  void flip(int i) {
    std::swap(bounds.min[i], bounds.max[i]) ;
    dx[i] *= -1 ;
  }

  void interp_corner(const Vector3i idx, const Vector3d pos,
		     const float* __restrict__ data,
		     float* __restrict__ result, int Nqty) const {

    // Fractional position in the cell
    Vector3f f((pos[0] - bounds.min[0])/dx[0] - idx[0],
	       (pos[1] - bounds.min[1])/dx[1] - idx[1],
	       (pos[2] - bounds.min[2])/dx[2] - idx[2]) ;

    Vector3i Nxc(Nx[0] + 1, Nx[1] + 1, Nx[2] + 1) ;
    
    int idx3D = ((idx[0]*Nxc[1] + idx[1])*Nxc[2] + idx[2]) * Nqty ;
    int i ;
    
    i = idx3D ;
    for (int n=0; n < Nqty; n++)
      result[n]  = (1-f[0])*(1-f[1])*(1-f[2]) * data[i + n] ;  
    i = idx3D + Nqty ;
    for (int n=0; n < Nqty; n++)
      result[n] += (1-f[0])*(1-f[1])*(  f[2]) * data[i + n] ;

    i = idx3D + Nqty*Nxc[2];
    for (int n=0; n < Nqty; n++)
      result[n] += (1-f[0])*(  f[1])*(1-f[2]) * data[i + n] ;
    i = idx3D + Nqty*(Nxc[2] + 1); 
    for (int n=0; n < Nqty; n++)
      result[n] += (1-f[0])*(  f[1])*(  f[2]) * data[i + n] ;

    
    i = idx3D + Nqty*(Nxc[2]*Nxc[1]) ;
    for (int n=0; n < Nqty; n++)
      result[n] += (  f[0])*(1-f[1])*(1-f[2]) * data[i + n] ;
    i = idx3D + Nqty*(Nxc[2]*Nxc[1] + 1);
    for (int n=0; n < Nqty; n++)
      result[n] += (  f[0])*(1-f[1])*(  f[2]) * data[i + n] ;

    i = idx3D + Nqty*Nxc[2]*(Nxc[1]+1);
    for (int n=0; n < Nqty; n++)
      result[n] += (  f[0])*(  f[1])*(1-f[2]) * data[i + n] ;
    i = idx3D + Nqty*(Nxc[2]*(Nxc[1]+1) + 1);
    for (int n=0; n < Nqty; n++)
      result[n] += (  f[0])*(  f[1])*(  f[2]) * data[i + n] ;
  }

private:
  Vector3i Nx ;
  Box bounds ;
  Vector3d dx ;
} ;


void _integrate_ray(Ray ray, const Grid& grid,
		    int Nfreq, double* __restrict__ I_nu, double* __restrict__ tau_nu,
		    const float* __restrict__ emiss, const float* __restrict__ kappa) {


  bool second_order = true ;

  // Zero out the result
  for (int n=0; n < Nfreq; n++)
    I_nu[n] = tau_nu[n] = 0 ;
    
  // Work out where the ray meets the grid
  RayBoxIntersection rbi(ray, grid.domain()) ;

  // Ray misses the domain
  if (rbi.dl <= 0) return ;

  // Find the first cell intersected by the ray
  Vector3d pos = ray.pos(rbi.l0) ;  
  Vector3i idx = grid.cell_index_safe(pos) ;

  // Check that we are going in the right direction
  if (!grid.in_domain(idx)) 
    for (int k(0); k < 3; k++) {
      if (idx[k] <  0            && ray.inv_dir[k] < 0) return ;
      if (idx[k] >= grid.size(k) && ray.inv_dir[k] > 0) return ;
    }

      
  // Walk along the ray until we hit the domain
  //  This protects against round-off placing us in the wrong cell
  while (!grid.in_domain(idx)) {
    rbi = RayBoxIntersection(ray, grid.cell_box(idx)) ;

    Vector3d ray_pos = ray.pos(rbi.l0) ;

    // Exit if we will miss the domain
    for (int k(0); k < 3; k++)
      if (ray.inv_dir[k] > 0) {
	if (ray_pos[k] > grid.domain().max[k]) return ;
      }
      else {
	if (ray_pos[k] < grid.domain().max[k]) return ;
      }

    for (int k(0); k < 3; k++)
      idx[k] += rbi.exit[k] ;
  }


  std::vector<float> e0, e1, k0, k1 ;
  if (second_order) {
    e0.resize(Nfreq) ;
    k0.resize(Nfreq) ;
    e1.resize(Nfreq) ;
    k1.resize(Nfreq) ;
  }
  
  if (second_order) {
    grid.interp_corner(idx, ray.pos(rbi.l0), emiss, e0.data(), Nfreq) ;
    grid.interp_corner(idx, ray.pos(rbi.l0), kappa, k0.data(), Nfreq) ;
  }
      
      
  // Now loop through the domain, building up the intensity
  do {
    int idx3d = grid.cell_index(idx) * Nfreq ;

    rbi = RayBoxIntersection(ray, grid.cell_box(idx)) ;

    if (second_order) {
      grid.interp_corner(idx, ray.pos(rbi.l1), emiss, e1.data(), Nfreq) ;
      grid.interp_corner(idx, ray.pos(rbi.l1), kappa, k1.data(), Nfreq) ;

      for (int n=0; n < Nfreq; n++) {
	double dtau = 0.5*(k0[n] + k1[n]) * rbi.dl ;
	double e_dl = 0.5*(e0[n] + e1[n]) * rbi.dl ;

	double e_tau = std::exp(-tau_nu[n]) ;
	    
	if (dtau == 0) {
	  I_nu[n] += e_tau * e_dl ;
	}
	else {
	  if (dtau > 1e-6) {
	    double edt = exp(-dtau) ;
	    double f0 = 1 - (1 - edt)/dtau ;
	    double f1 = (1 - edt) / dtau - edt ;

	    I_nu[n] += e_tau * std::min(e_dl, (e0[n]/k0[n])*f0 + (e1[n]/k1[n])*f1) ;
	  }
	  else {
	    I_nu[n] += e_tau * 
	      std::min(e_dl, 0.25 * (e0[n]*(1 + k1[n]/k0[n]) + e1[n]*(1 + k0[n]/k1[n]))*rbi.dl) ;
	  }
	      
	  tau_nu[n] += dtau ;
	}
	    
	e0[n] = e1[n] ;
	k0[n] = k1[n] ;
      }
	  
    }
    else {	
      for (int n=0; n < Nfreq; n++) {

	double e = emiss[idx3d + n] ;
	double k = kappa[idx3d + n] ;      

	double dtau = k * rbi.dl ;

	if (k > 0)
	  I_nu[n] -= (e/k) * std::exp(-tau_nu[n]) * std::expm1(-dtau) ;
	else
	  I_nu[n] += e * rbi.dl ;

	tau_nu[n] += dtau ;
      }
    }

	
    // Find the next cell
    for (int k(0); k < 3; k++)
      idx[k] += rbi.exit[k] ;
  } while (grid.in_domain(idx)) ;
}


/* Compute an image by ray tracing through a volume.
 *
 * The image follows FFTW conventions, where the origin of the domain is at floor(image_Nx/2).
 *
 * args:
 *  image_Nx    : Size of image (image_Nx*image_Nx)
 *  image_dx    : Image cell size (assumed square)
 *  alpha       : First Euler rotation angle (about Z)
 *  beta        : Second Euler roation angle (about Y)
 *  gamma       : Third Euler rotation angle (about Z again)
 *  Nfreq       : Number of frequency bins
 *  I_nu        : Computed intensity (image_Nx*image_Nx*Nfreq)
 *  tau_nu      : Computed optical depth (image_Nx*image_Nx*Nfreq)
 *  grid_Nx     : 3D index specifiying the logical size of the emission/opacity grid
 *  grid_origin : 3D vector specifying the left corner of the grid
 *  grid_dx     : 3D vector of grid sizes
 *  emiss       : Emissivity data for each cell and frequency, specified at cell corners
 *  kappa       : Opacity data for each cell and frequency, specified at cell corners
 */
void RayTraceImage(int image_Nx, double image_dx, double alpha, double beta, double gamma,
		   int Nfreq, double* __restrict__ I_nu, double* __restrict__ tau_nu,
		   Vector3i grid_Nx, Vector3d grid_origin, Vector3d grid_dx,
		   const float* __restrict__ emiss, const float* __restrict__ kappa) {

  using std::cos ;
  using std::sin ;

  /* Compute the ray direction from angles
   *   TODO: Check definition against galario 
   */
  double
    c1 = cos(alpha), s1 = sin(alpha),
    c2 = cos(beta),  s2 = sin(beta),
    c3 = cos(gamma), s3 = sin(gamma) ;
  
  Vector3d ray_dir(s2*c1, s2*s1, c2) ;
  
  // Setup the grid index space.
  Grid grid(grid_Nx, grid_origin, grid_dx) ;

  // Flip the grid (left-right) to match the rays' direction
  // Division handles negative zero
  for (int k=0; k < 3; k++)
    if (1. / ray_dir[k] < 0)
      grid.flip(k)  ;

	  
  for (int i(0); i < image_Nx; i++)
    for (int j(0); j < image_Nx; j++) {
      int idx2D = (i*image_Nx + j) * Nfreq ;
      
      // Compute a position along the ray via rotation
      double x0 = (i - image_Nx / 2) * image_dx ;
      double y0 = (j - image_Nx / 2) * image_dx ;

      Vector3d ray_pos(+ x0*(c1*c2*c3-s1*s3) - y0*(c3*s1+c1*c2*s3),
		       + x0*(c1*s3+c2*c3*s1) + y0*(c1*c3-c2*s1*s3),
		       - x0*(c3*s2)          + y0*(s2*s3) ) ;

      Ray ray(ray_pos, ray_dir) ;

      // Do the integration along the ray
      _integrate_ray(ray, grid,
		     Nfreq, I_nu + idx2D, tau_nu + idx2D,
		     emiss, kappa) ;
    }
}



		     
int main() {

  int Nx = 128;
  int Nfreq = 3 ;

  int Npix = 64 ;
  
  double dx = 1. / Nx ;

  Grid g({Nx, Nx, Nx}, {0.,0.,0.}, {dx, dx, dx}) ;

  int Ntot = (Nx+1)*(Nx+1)*(Nx+1) ;
  std::vector<float> emiss(Ntot*Nfreq, 0);
  std::vector<float> kappa(Ntot*Nfreq, 0);

  std::vector<double> I_nu(Npix*Npix*Nfreq, 0);
  std::vector<double> tau_nu(Npix*Npix*Nfreq, 0);

  for (int j=0; j <Ntot; j++) {
    int i = j*Nfreq ;
    
    emiss[i] = 1 ;
    kappa[i] = 0 ;
    
    emiss[i+1] = 0 ;
    kappa[i+1] = 1 ;
    
    emiss[i+2] = 1 ;
    kappa[i+2] = 1 ;	  
  }

  {
    RayTraceImage(Npix, 2./Npix, 0., 0., 0.,
		  Nfreq, I_nu.data(), tau_nu.data(),
		  Vector3i(Nx,Nx,Nx), Vector3d(-0.5,-0.5,-0.5), Vector3d(dx,dx,dx),
		  emiss.data(), kappa.data()) ;

    int i = Npix/2 ;
    int j = Npix/2 ;
    int idx = (i*Npix + j) * Nfreq ;
    
    std::cout << "I_nu: " ;
    for (int i=0; i < Nfreq; i++)
      std::cout << I_nu[idx+i] << " " ;
    std::cout << "\n" ;

    std::cout << "tau_nu: " ;
    for (int i=0; i < Nfreq; i++)
      std::cout << tau_nu[idx+i] << " " ;
    std::cout << "\n" ;
    std::cout << std::endl;
  }

  {
    RayTraceImage(Npix, 2./Npix, M_PI/2, M_PI/4, 0., 
		  Nfreq, I_nu.data(), tau_nu.data(),
		  Vector3i(Nx,Nx,Nx), Vector3d(-0.5,-0.5,-0.5), Vector3d(dx,dx,dx),
		  emiss.data(), kappa.data()) ;

    int i = Npix/2 ;
    int j = Npix/2 ;
    int idx = (i*Npix + j) * Nfreq ;
    
    std::cout << "I_nu: " ;
    for (int i=0; i < Nfreq; i++)
      std::cout << I_nu[idx+i] << " " ;
    std::cout << "\n" ;

    std::cout << "tau_nu: " ;
    for (int i=0; i < Nfreq; i++)
      std::cout << tau_nu[idx+i] << " " ;
    std::cout << "\n\n" ;

  }

  {
    RayTraceImage(Npix, 2./Npix, M_PI, 0., 0.,
		  Nfreq, I_nu.data(), tau_nu.data(),
		  Vector3i(Nx,Nx,Nx), Vector3d(-0.5,-0.5,-0.5), Vector3d(dx,dx,dx),
		  emiss.data(), kappa.data()) ;
    
    int i = Npix/2 ;
    int j = Npix/2 ;
    int idx = (i*Npix + j) * Nfreq ;
    
    std::cout << "I_nu: " ;
    for (int i=0; i < Nfreq; i++)
      std::cout << I_nu[idx+i] << " " ;
    std::cout << "\n" ;

    std::cout << "tau_nu: " ;
    for (int i=0; i < Nfreq; i++)
      std::cout << tau_nu[idx+i] << " " ;
    std::cout << "\n\n" ;
  }

  {
    RayTraceImage(Npix, 2./Npix, M_PI/2., M_PI/4, 0.,
		  Nfreq, I_nu.data(), tau_nu.data(),
		  Vector3i(Nx,Nx,Nx), Vector3d(-0.5,-0.5,-0.5), Vector3d(dx,dx,dx),
		  emiss.data(), kappa.data()) ;
    
    int i = (-0.5 + 1) * Npix/2 ;
    int j = (-0.0 + 1) * Npix/2 ;
    int idx = (i*Npix + j) * Nfreq ;
    
    std::cout << "I_nu: " ;
    for (int i=0; i < Nfreq; i++)
      std::cout << I_nu[idx+i] << " " ;
    std::cout << "\n" ;

    std::cout << "tau_nu: " ;
    for (int i=0; i < Nfreq; i++)
      std::cout << tau_nu[idx+i] << " " ;
    std::cout << "\n\n" ;

    std::ofstream f("tmp") ;
    for (int i=0; i < Npix; i++) 
      for (int j=0; j < Npix; j++) {
	idx = (i*Npix + j) * Nfreq ;
	f << (i - Npix/2) * 2./Npix << " " << (j - Npix/2) * 2./Npix << " "  ;
	for (int n=0; n < Nfreq; n++)
	  f << I_nu[idx+n] << " " << tau_nu[idx+n] << " " ;
	f << "\n" ;	
      }
    
  }
    
}
