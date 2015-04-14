/* Low-level dynamics integrators
 *
 * Written by Konrad Hinsen
 */

#include "MMTK/universe.h"
#include "MMTK/forcefield.h"
#include "MMTK/trajectory.h"

#include <unistd.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#define DEBUG 0



/* Global variables */

double kB, temperature_factor;

/* Allocate and initialize Output variable descriptors */

static PyTrajectoryVariable *
get_data_descriptors(int n, PyUniverseSpecObject *universe_spec,
		     PyArrayObject *configuration, PyArrayObject *velocities,
		     PyArrayObject *gradients, PyArrayObject *masses,
		     int *ndf,
		     double *time, double *p_energy, double *k_energy,
		     double *n_energy, double *a_energy,
		     double *temperature, double *xi,
		     double *pressure, double *volume, double *alpha,
		     double *box_size)
{
  PyTrajectoryVariable *vars = (PyTrajectoryVariable *)
                               malloc((n+1)*sizeof(PyTrajectoryVariable));
  int i = 0;
  if (vars == NULL)
    return (PyTrajectoryVariable *)PyErr_NoMemory();
  if (time != NULL && i < n) {
    vars[i].name = "time";
    vars[i].text = "Time: %lf\n";
    vars[i].unit = time_unit_name;
    vars[i].type = PyTrajectory_Scalar;
    vars[i].class = PyTrajectory_Time;
    vars[i].value.dp = time;
    i++;
  }
  if (p_energy != NULL && i < n) {
    vars[i].name = "potential_energy";
    vars[i].text = "Potential energy: %lf, ";
    vars[i].unit = energy_unit_name;
    vars[i].type = PyTrajectory_Scalar;
    vars[i].class = PyTrajectory_Energy;
    vars[i].value.dp = p_energy;
    i++;
  }
  if (k_energy != NULL && i < n) {
    vars[i].name = "kinetic_energy";
    vars[i].text = "Kinetic energy: %lf\n";
    vars[i].unit = energy_unit_name;
    vars[i].type = PyTrajectory_Scalar;
    vars[i].class = PyTrajectory_Energy;
    vars[i].value.dp = k_energy;
    i++;
  }
  if (n_energy != NULL && i < n) {
    vars[i].name = "nose_energy";
    vars[i].text = "Nose energy: %lf\n";
    vars[i].unit = energy_unit_name;
    vars[i].type = PyTrajectory_Scalar;
    vars[i].class = PyTrajectory_Energy;
    vars[i].value.dp = n_energy;
    i++;
  }
  if (a_energy != NULL && i < n) {
    vars[i].name = "andersen_energy";
    vars[i].text = "Andersen energy: %lf\n";
    vars[i].unit = energy_unit_name;
    vars[i].type = PyTrajectory_Scalar;
    vars[i].class = PyTrajectory_Energy;
    vars[i].value.dp = a_energy;
    i++;
  }
  if (temperature != NULL && i < n) {
    vars[i].name = "temperature";
    vars[i].text = "Temperature: %lf\n";
    vars[i].unit = temperature_unit_name;
    vars[i].type = PyTrajectory_Scalar;
    vars[i].class = PyTrajectory_Thermodynamic;
    vars[i].value.dp = temperature;
    i++;
  }
  if (xi != NULL && i < n) {
    vars[i].name = "thermostat_coordinate";
    vars[i].text = "Thermostat coordinate: %lf\n";
    vars[i].unit = frequency_unit_name;
    vars[i].type = PyTrajectory_Scalar;
    vars[i].class = PyTrajectory_Auxiliary;
    vars[i].value.dp = xi;
    i++;
  }
  if (pressure != NULL && i < n) {
    vars[i].name = "pressure";
    vars[i].text = "Pressure: %lf\n";
    vars[i].unit = pressure_unit_name;
    vars[i].type = PyTrajectory_Scalar;
    vars[i].class = PyTrajectory_Thermodynamic;
    vars[i].value.dp = pressure;
    i++;
  }
  if (alpha != NULL && i < n) {
    vars[i].name = "barostat_coordinate";
    vars[i].text = "Barostat coordinate: %lf\n";
    vars[i].unit = frequency_unit_name;
    vars[i].type = PyTrajectory_Scalar;
    vars[i].class = PyTrajectory_Auxiliary;
    vars[i].value.dp = alpha;
    i++;
  }
  if (volume != NULL && i < n) {
    vars[i].name = "volume";
    vars[i].text = "Volume: %lf\n";
    vars[i].unit = volume_unit_name;
    vars[i].type = PyTrajectory_Scalar;
    vars[i].class = PyTrajectory_Thermodynamic;
    vars[i].value.dp = volume;
    i++;
  }
  if (configuration != NULL && i < n) {
    vars[i].name = "configuration";
    vars[i].text = "Configuration:\n";
    vars[i].unit = length_unit_name;
    vars[i].type = PyTrajectory_ParticleVector;
    vars[i].class = PyTrajectory_Configuration;
    vars[i].value.array = configuration;
    i++;
  }
  if (box_size != NULL && i < n) {
    vars[i].name = "box_size";
    vars[i].text = "Box size:";
    vars[i].unit = length_unit_name;
    vars[i].type = PyTrajectory_BoxSize;
    vars[i].class = PyTrajectory_Configuration;
    vars[i].value.dp = box_size;
    vars[i].length = universe_spec->geometry_data_length;
    i++;
  }
  if (velocities != NULL && i < n) {
    vars[i].name = "velocities";
    vars[i].text = "Velocities:\n";
    vars[i].unit = velocity_unit_name;
    vars[i].type = PyTrajectory_ParticleVector;
    vars[i].class = PyTrajectory_Velocities;
    vars[i].value.array = velocities;
    i++;
  }
  if (gradients != NULL && i < n) {
    vars[i].name = "gradients";
    vars[i].text = "Energy gradients:\n";
    vars[i].unit = energy_gradient_unit_name;
    vars[i].type = PyTrajectory_ParticleVector;
    vars[i].class = PyTrajectory_Gradients;
    vars[i].value.array = gradients;
    i++;
  }
  if (masses != NULL && i < n) {
    vars[i].name = "masses";
    vars[i].text = "Masses:\n";
    vars[i].unit = mass_unit_name;
    vars[i].type = PyTrajectory_ParticleScalar;
    vars[i].class = PyTrajectory_Internal;
    vars[i].value.array = masses;
    i++;
  }
  if (ndf != NULL && i < n) {
    vars[i].name = "degrees_of_freedom";
    vars[i].text = "Degrees of freedom: %d\n";
    vars[i].unit = "";
    vars[i].type = PyTrajectory_IntScalar;
    vars[i].class = PyTrajectory_Internal;
    vars[i].value.ip = ndf;
    i++;
  }
  vars[i].name = NULL;
  return vars;
}

/* Constraints: SHAKE */

static void
shake(long *const_pairs, int from, int to, vector3 *x,
      double *m, vector3 *const_vect, double *const_dist,
      distance_fn *d_fn, double *d_data)
{
  const int max_iter = 500;
  const double tolerance = 1.e-8;
  int i, j;

  for (i = 0; i < max_iter; i++) {
    double max_dev = 0.;
    for (j = from; j < to; j++) {
      vector3 d;
      double s, dev, l;
      int i1, i2;
      i1 = const_pairs[2*j];
      i2 = const_pairs[2*j+1];
      (*d_fn)(d, x[i1], x[i2], d_data);
      s = 0.5*(vector_length_sq(d)-const_dist[j]);
      dev = fabs(s)/const_dist[j];
      if (dev > max_dev) max_dev = dev;
      if (dev > tolerance) {
	l = -s*m[i1]*m[i2]/((m[i1]+m[i2])*dot(d, const_vect[j]));
	x[i1][0] -= l*const_vect[j][0]/m[i1];
	x[i1][1] -= l*const_vect[j][1]/m[i1];
	x[i1][2] -= l*const_vect[j][2]/m[i1];
	x[i2][0] += l*const_vect[j][0]/m[i2];
	x[i2][1] += l*const_vect[j][1]/m[i2];
	x[i2][2] += l*const_vect[j][2]/m[i2];
      }
    }
    if (max_dev < tolerance) {
      break;
    }
  }
#if 0
  if (i == max_iter)
    printf("No convergence in SHAKE\n");
  else
    printf("%d SHAKE iterations\n", i);
#endif
}


/* Constraints: Lagrange multipliers */

typedef struct {
  double multiplier[4];  /* mu1, mu2, g, gdot */
  double diag;
  double scratch;
} projection_data;

enum multipliers {MU1=0, MU2=1, G=2, GDOT=3};

static void
project(int n_const, long *const_pairs,
	double *const_dist, vector3* const_vect,
	projection_data *p, int mindex, double *mass,
	vector3 *in, vector3 *out, int natoms)
{
  const int max_iter = 1000;
  const double tolerance = 1.e-8;
  int converged, niter;
  int i;

  for (i = 0; i < natoms; i++)
    out[i][0] = out[i][1] = out[i][2] = 0.;
  for (i = 0; i < n_const; i++) {
    int a1 = const_pairs[2*i];
    int a2 = const_pairs[2*i+1];
    double m = p[i].multiplier[mindex];
    out[a2][0] += m*const_vect[i][0]/mass[a2];
    out[a2][1] += m*const_vect[i][1]/mass[a2];
    out[a2][2] += m*const_vect[i][2]/mass[a2];
    out[a1][0] -= m*const_vect[i][0]/mass[a1];
    out[a1][1] -= m*const_vect[i][1]/mass[a1];
    out[a1][2] -= m*const_vect[i][2]/mass[a1];
  }
  niter = 0;
  while (1) {
    converged = 0;
    for (i = 0; i < n_const; i++) {
      int a1 = const_pairs[2*i];
      int a2 = const_pairs[2*i+1];
      double dm;
      if (mindex == G)
	dm = (- const_dist[i]
	      - const_vect[i][0]*(out[a2][0]-out[a1][0])
	      - const_vect[i][1]*(out[a2][1]-out[a1][1])
	      - const_vect[i][2]*(out[a2][2]-out[a1][2]))/p[i].diag;
      else
	dm =   (const_vect[i][0]*(in[a2][0]-in[a1][0]-out[a2][0]+out[a1][0])
	      + const_vect[i][1]*(in[a2][1]-in[a1][1]-out[a2][1]+out[a1][1])
	      + const_vect[i][2]*(in[a2][2]-in[a1][2]-out[a2][2]+out[a1][2]))
	      / p[i].diag;
      if (fabs(dm) < tolerance*fabs(p[i].multiplier[mindex]))
	converged++;
      p[i].multiplier[mindex] += dm;
      out[a2][0] += dm*const_vect[i][0]/mass[a2];
      out[a2][1] += dm*const_vect[i][1]/mass[a2];
      out[a2][2] += dm*const_vect[i][2]/mass[a2];
      out[a1][0] -= dm*const_vect[i][0]/mass[a1];
      out[a1][1] -= dm*const_vect[i][1]/mass[a1];
      out[a1][2] -= dm*const_vect[i][2]/mass[a1];
    }
    niter++;
    if (converged == n_const) {
#if DEBUG
      printf("%d projection iterations for %d\n", niter, mindex);
#endif
      break;
    }
    if (niter > max_iter) {
#if DEBUG
      printf("No convergence in projection for %d\n", mindex);
#endif
      break;
    }
  }
}

static void
project2(int n_const, long *const_pairs,
	 double *const_dist, vector3* const_vect,
	 projection_data *p, int mindex, double *mass,
	 double *in, vector3 *out, int natoms)
{
  const int max_iter = 1000;
  const double tolerance = 1.e-8;
  int converged, niter;
  int i;

  for (i = 0; i < natoms; i++)
    out[i][0] = out[i][1] = out[i][2] = 0.;
  for (i = 0; i < n_const; i++) {
    int a1 = const_pairs[2*i];
    int a2 = const_pairs[2*i+1];
    double m = p[i].multiplier[mindex];
    out[a2][0] += m*const_vect[i][0]/mass[a2];
    out[a2][1] += m*const_vect[i][1]/mass[a2];
    out[a2][2] += m*const_vect[i][2]/mass[a2];
    out[a1][0] -= m*const_vect[i][0]/mass[a1];
    out[a1][1] -= m*const_vect[i][1]/mass[a1];
    out[a1][2] -= m*const_vect[i][2]/mass[a1];
  }
  niter = 0;
  while (1) {
    converged = 0;
    for (i = 0; i < n_const; i++) {
      int a1 = const_pairs[2*i];
      int a2 = const_pairs[2*i+1];
      double dm;
      dm =   ((const_vect[i][0]*(out[a2][0]-out[a1][0])
	       + const_vect[i][1]*(out[a2][1]-out[a1][1])
	       + const_vect[i][2]*(out[a2][2]-out[a1][2])) + in[i])
 	     / p[i].diag;
      if (fabs(dm) < tolerance*fabs(p[i].multiplier[mindex]))
	converged++;
      p[i].multiplier[mindex] -= dm;
      out[a2][0] -= dm*const_vect[i][0]/mass[a2];
      out[a2][1] -= dm*const_vect[i][1]/mass[a2];
      out[a2][2] -= dm*const_vect[i][2]/mass[a2];
      out[a1][0] += dm*const_vect[i][0]/mass[a1];
      out[a1][1] += dm*const_vect[i][1]/mass[a1];
      out[a1][2] += dm*const_vect[i][2]/mass[a1];
    }
    niter++;
    if (converged == n_const) {
#if DEBUG
      printf("%d projection iterations for %d\n", niter, mindex);
#endif
      break;
    }
    if (niter > max_iter) {
#if DEBUG
      printf("No convergence in projection for %d\n", mindex);
#endif
      break;
    }
  }
}

static void
mult_by_h_plus_one(vector3 *in, vector3 *out, int natoms, double *mass,
		   long *const_pairs, projection_data *p, int n_const)
{
  int i;

  for (i = 0; i < natoms; i++)
    vector_copy(out[i], in[i]);
  for (i = 0; i < n_const; i++) {
    int a1 = const_pairs[2*i];
    int a2 = const_pairs[2*i+1];
    double g = p[i].multiplier[0];
    out[a1][0] += g*(in[a1][0]-in[a2][0])/mass[a1];
    out[a1][1] += g*(in[a1][1]-in[a2][1])/mass[a1];
    out[a1][2] += g*(in[a1][2]-in[a2][2])/mass[a1];
    out[a2][0] += g*(in[a2][0]-in[a1][0])/mass[a2];
    out[a2][1] += g*(in[a2][1]-in[a1][1])/mass[a2];
    out[a2][2] += g*(in[a2][2]-in[a1][2])/mass[a2];
  }
}

/* Runge-Kutta-Merson integrator */

static PyObject *
integrateRKM(PyObject *dummy, PyObject *args)
{
  PyObject *universe;
  PyUniverseSpecObject *universe_spec;
  PyArrayObject *configuration;
  PyArrayObject *velocities;
  PyArrayObject *masses;
  PyArrayObject *fixed;
  PyArrayObject *gradients;
  PyArrayObject *constraints, *constraint_distances_squared, *c_blocks;
  PyArrayObject *t_parameters, *t_coordinates;
  PyArrayObject *b_parameters, *b_coordinates;
  PyListObject *spec_list;
  PyFFEvaluatorObject *evaluator;
  PyTrajectoryOutputSpec *output;
  PyTrajectoryVariable *data_descriptors = NULL;
  double delta_t, dth;
  int first_step, last_step;
  char *description;
  vector3 *x, *v, *f;
  double *m, *const_dist;
  long *fix, *const_pairs, *const_blocks;
  vector3 *const_vect;
  vector3 *scratch;
  vector3 *xp, *xph, *vh, *v1, *v2, *xold;
  double *temp;
  projection_data *pdata;
  double time;
  energy_data p_energy;
  double k_energy, n_energy, a_energy;
  double temperature, volume, const_virial1, const_virial2, pressure;
  double t_temp, t_tau, t_mass, t_energy, *t_xi, *t_lns;
  double b_press, b_tau, b_mass, *b_alpha;
  double dxi, dalpha, factor1, factor2;
  int atoms, n_const, n_const_blocks, df;
  int pressure_available, thermostat, barostat;
  int i, j;

  /* Parse and check arguments */
  if (!PyArg_ParseTuple(args, "OO!O!O!O!O!O!O!O!O!O!O!O!diiO!s", &universe,
			&PyArray_Type, &configuration,
			&PyArray_Type, &velocities,
			&PyArray_Type, &masses,
			&PyArray_Type, &fixed,
			&PyFFEvaluator_Type, &evaluator,
			&PyArray_Type, &constraints,
			&PyArray_Type, &constraint_distances_squared,
			&PyArray_Type, &c_blocks,
			&PyArray_Type, &t_parameters,
			&PyArray_Type, &t_coordinates,
			&PyArray_Type, &b_parameters,
			&PyArray_Type, &b_coordinates,
			&delta_t, &first_step, &last_step,
			&PyList_Type, &spec_list,
			&description))
    return NULL;
  universe_spec = (PyUniverseSpecObject *)
                   PyObject_GetAttrString(universe, "_spec");
  if (universe_spec == NULL)
    return NULL;

  /* Create gradient array */
#if defined(NUMPY)
  gradients = (PyArrayObject *)PyArray_Copy(configuration);
#else
  gradients = (PyArrayObject *)PyArray_FromDims(configuration->nd,
						configuration->dimensions,
						PyArray_DOUBLE);
#endif
  if (gradients == NULL)
    return NULL;

  /* Set some convenient variables */
  atoms = configuration->dimensions[0];
  n_const = constraints->dimensions[0];
  n_const_blocks = c_blocks->dimensions[0]-1;
  thermostat = (t_parameters->dimensions[0] > 0);
  barostat = (b_parameters->dimensions[0] > 0);
  dth = 0.5*delta_t;
  x = (vector3 *)configuration->data;
  v = (vector3 *)velocities->data;
  f = (vector3 *)gradients->data;
  m = (double *)masses->data;
  fix = (long *)fixed->data;
  const_pairs = (long *)constraints->data;
  const_dist = (double *)constraint_distances_squared->data;
  const_blocks = (long *)c_blocks->data;

  /* Calculate number of degrees of freedom */
  df = 3*atoms;
  for (j = 0; j < atoms; j++) {
    if (fix[j])
      df -= 3;
  }
  for (j = 0; j < n_const; j++) {
    if (fix[const_pairs[2*j]] || fix[const_pairs[2*j+1]]) {
      PyErr_SetString(PyExc_ValueError,
		      "distance constraint on a fixed atom");
      return NULL;
    }
    df--;
  }

  /* Thermostat and barostat information.
     Note: the coordinates are guaranteed to be zero when there is
     no thermostat and/or barostat. */
  if (thermostat) {
    t_temp = *(double *)t_parameters->data;
    t_tau = *(((double *)t_parameters->data)+1);
    t_energy = df*t_temp*kB;
    t_mass = t_energy*t_tau*t_tau;
  }
  else {
    t_mass = 0.;
    t_energy = 0.;
    t_temp = 0.; /* unused, initialize just to make gcc happy */
    t_mass = 0;  /* unused, initialize just to make gcc happy */
  }
  t_xi = (double *)t_coordinates->data;
  t_lns = t_xi + 1;
  if (barostat) {
    b_press = *(double *)b_parameters->data;
    b_tau = *(((double *)b_parameters->data)+1);
  }
  else {
    b_mass = 0.;
    b_press = 0.;
    b_tau = 0.; /* unused, initialize just to make gcc happy */
  }
  b_alpha = (double *)b_coordinates->data;

  /* Allocate arrays for temporary data */
  pdata = NULL;
  if (n_const > 0) {
    scratch = (vector3 *)malloc((5*atoms+n_const)*sizeof(vector3));
    if (scratch == NULL) {
      PyErr_NoMemory();
      goto error2;
    }
    v1 = scratch;
    v2 = v1 + atoms;
    xp = v2 + atoms;
    xph = xp + atoms;
    vh = xph + atoms;
    const_vect = vh + atoms;
    xold = v1;
    temp = (double *)v2;
    pdata = (projection_data *)malloc(n_const*sizeof(projection_data));
    if (pdata == NULL) {
      PyErr_NoMemory();
      goto error2;
    }
  }
  else {
    scratch = NULL;
    v1 = v2 = NULL;
    xp = xph = x;
    vh = v;
    xold = NULL;
    temp = NULL;
    const_vect = NULL;  /* unused, initialize just to make gcc happy */
  }

  /* Enforce constraints and initialize constraint data */
  Py_BEGIN_ALLOW_THREADS;
  PyUniverseSpec_StateLock(universe_spec, -1);
  universe_spec->correction_function(x, atoms, universe_spec->geometry_data);
  if (n_const > 0) {
    for (j = 0; j < n_const; j++) {
      universe_spec->distance_function(const_vect[j],
				       x[const_pairs[2*j]],
				       x[const_pairs[2*j+1]],
				       universe_spec->geometry_data);
      pdata[j].multiplier[MU1] = 0.;
      pdata[j].multiplier[MU2] = 0.;
      pdata[j].multiplier[G] = 0.;
      pdata[j].multiplier[GDOT] = 0.;
      pdata[j].diag = const_dist[j] *
	(1./m[const_pairs[2*j]] + 1./m[const_pairs[2*j+1]]);
    }
    for (j = 0; j < n_const_blocks; j++)
      shake(const_pairs, const_blocks[j], const_blocks[j+1],
	    x, m, const_vect, const_dist, universe_spec->distance_function,
	    universe_spec->geometry_data);
    for (j = 0; j < n_const; j++)
      universe_spec->distance_function(const_vect[j],
				       x[const_pairs[2*j]],
				       x[const_pairs[2*j+1]],
				       universe_spec->geometry_data);
    if (barostat) {
      project(n_const, const_pairs, const_dist, const_vect,
	      pdata, G, m, NULL, xp, atoms);
      mult_by_h_plus_one(v, vh, atoms, m, const_pairs, pdata, n_const);
      project(n_const, const_pairs, const_dist, const_vect,
	      pdata, MU2, m, vh, v2, atoms);
      for (j = 0; j < n_const; j++) {
	int a1 = const_pairs[2*j];
	int a2 = const_pairs[2*j+1];
	vector3 diff1, diff2;
	diff1[0] = (*b_alpha)*(const_vect[j][0] + xp[a2][0]-xp[a1][0])
	           + v[a2][0] - v[a1][0];
	diff1[1] = (*b_alpha)*(const_vect[j][1] + xp[a2][1]-xp[a1][1])
	           + v[a2][1] - v[a1][1];
	diff1[2] = (*b_alpha)*(const_vect[j][2] + xp[a2][2]-xp[a1][2])
	           + v[a2][2] - v[a1][2];
	diff2[0] = v[a2][0] - v[a1][0];
	diff2[1] = v[a2][1] - v[a1][1];
	diff2[2] = v[a2][2] - v[a1][2];
	temp[j] = dot(diff1, diff2);
	diff1[0] = f[a1][0]/m[a1] - f[a2][0]/m[a1];
	diff1[1] = f[a1][1]/m[a1] - f[a2][1]/m[a1];
	diff1[2] = f[a1][2]/m[a1] - f[a2][2]/m[a1];
	temp[j] += dot(const_vect[j], diff1);
	temp[j] *= dth/(dth*(*t_xi)-1.);
      }
      project2(n_const, const_pairs, const_dist, const_vect,
	       pdata, MU1, m, temp, v1, atoms);
    }
    const_virial1 = const_virial2 = 0.;
    for (j = 0; j < n_const; j++) {
      const_virial1 += pdata[j].multiplier[MU1]*const_dist[j];
      const_virial2 += pdata[j].multiplier[MU2]*const_dist[j];
    }
  }
  else
    const_virial1 = const_virial2 = 0.;
  PyUniverseSpec_StateLock(universe_spec, -2);
  Py_END_ALLOW_THREADS;

  /* Initial force calculation */
  p_energy.gradients = (PyObject *)gradients;
  p_energy.gradient_fn = NULL;
  p_energy.force_constants = NULL;
  p_energy.fc_fn = NULL;
#ifdef WITH_THREAD
  evaluator->tstate_save = PyEval_SaveThread();
#endif
  PyUniverseSpec_StateLock(universe_spec, 1);
  (*evaluator->eval_func)(evaluator, &p_energy, configuration, 0);
  PyUniverseSpec_StateLock(universe_spec, 2);
#ifdef WITH_THREAD
  PyEval_RestoreThread(evaluator->tstate_save);
#endif
  if (p_energy.error)
    goto error2;

  /* Check if pressure can be calculated */
  Py_BEGIN_ALLOW_THREADS;
  PyUniverseSpec_StateLock(universe_spec, 1);
  volume = universe_spec->volume_function(1., universe_spec->geometry_data);
  PyUniverseSpec_StateLock(universe_spec, 2);
  Py_END_ALLOW_THREADS;
  pressure_available = volume > 0. && p_energy.virial_available;
  if (thermostat && barostat)
    b_mass = df*t_temp*kB*b_tau*b_tau/(volume*volume);

  /* Initialize output */
  data_descriptors =
     get_data_descriptors(9 + pressure_available
			  + 2*thermostat
			  + 3*barostat
			  + (universe_spec->geometry_data_length > 0),
			  universe_spec,
			  configuration, velocities,
			  gradients, masses, &df,
			  &time, &p_energy.energy, &k_energy,
			  thermostat ? &n_energy : NULL,
			  barostat ? &a_energy : NULL,
			  &temperature,
			  thermostat ? t_xi : NULL,
			  pressure_available ? &pressure:NULL,
			  barostat ? &volume : NULL,
			  barostat ? b_alpha : NULL,
			  (universe_spec->geometry_data_length > 0) ?
			              universe_spec->geometry_data : NULL);
  if (data_descriptors == NULL)
    goto error2;
  output = PyTrajectory_OutputSpecification(universe, spec_list,
					    description,
					    data_descriptors);
  if (output == NULL)
    goto error2;

  evaluator->tstate_save = PyEval_SaveThread();
  
  /* Get write access for the integration, switching to
     read access only during energy evaluation */
  PyUniverseSpec_StateLock(universe_spec, -1);

  /** Main integration loop **/
  time  = first_step*delta_t;
  for (i = first_step; i < last_step; i++) {

    /* Calculation of thermodynamic properties */
    k_energy = 0.;
    for (j = 0; j < atoms; j++)
      if (!fix[j])
	k_energy += m[j]*dot(v[j], v[j]);
    k_energy *= 0.5;
    n_energy = 0.5*(*t_xi)*(*t_xi)*t_mass*exp(2*(*t_lns))
                 + t_energy*(*t_lns);
    a_energy = 4.5*b_mass*volume*volume*(*b_alpha)*(*b_alpha) + volume*b_press;
    temperature = 2.*k_energy*temperature_factor/df;
    pressure = (2.*k_energy+p_energy.virial
		+ ((*t_xi)-1./dth)*const_virial1
		+ (*b_alpha)*const_virial2)/(3.*volume);
    if (barostat && !thermostat)
      b_mass = 0.5*k_energy*b_tau*b_tau/(volume*volume);
    /* Trajectory and log output */
    if (PyTrajectory_Output(output, i, data_descriptors,
			    &evaluator->tstate_save) == -1) {
      PyUniverseSpec_StateLock(universe_spec, -2);
      PyEval_RestoreThread(evaluator->tstate_save);
      goto error;
    }

    /* First part of integration step */
    if (barostat) {
      if (n_const > 0) {
	project(n_const, const_pairs, const_dist, const_vect,
		pdata, G, m, NULL, xp, atoms);
	for (j = 0; j < atoms; j++)
	  v1[j][0] = v1[j][1] = v1[j][2] = 0.;
	for (j = 0; j < n_const; j++) {
	  int a1 = const_pairs[2*j];
	  int a2 = const_pairs[2*j+1];
	  vector3 diff1, diff2;
	  diff1[0] = (*b_alpha)*(const_vect[j][0] + xp[a2][0]-xp[a1][0])
	              + v[a2][0] - v[a1][0];
	  diff1[1] = (*b_alpha)*(const_vect[j][1] + xp[a2][1]-xp[a1][1])
	              + v[a2][1] - v[a1][1];
	  diff1[2] = (*b_alpha)*(const_vect[j][2] + xp[a2][2]-xp[a1][2])
	              + v[a2][2] - v[a1][2];
	  diff2[0] = xp[a2][0] - xp[a1][0];
	  diff2[1] = xp[a2][1] - xp[a1][1];
	  diff2[2] = xp[a2][2] - xp[a1][2];
	  temp[j] = dot(diff1, diff2);
	  v1[a2][0] += diff1[0]*pdata[j].multiplier[G]/m[a2];
	  v1[a2][1] += diff1[1]*pdata[j].multiplier[G]/m[a2];
	  v1[a2][2] += diff1[2]*pdata[j].multiplier[G]/m[a2];
	  v1[a1][0] -= diff1[0]*pdata[j].multiplier[G]/m[a1];
	  v1[a1][1] -= diff1[1]*pdata[j].multiplier[G]/m[a1];
	  v1[a1][2] -= diff1[2]*pdata[j].multiplier[G]/m[a1];
	}
	for (j = 0; j < n_const; j++) {
	  int a1 = const_pairs[2*j];
	  int a2 = const_pairs[2*j+1];
	  vector3 diff;
	  diff[0] = v1[a2][0] - v1[a1][0];
	  diff[1] = v1[a2][1] - v1[a1][1];
	  diff[2] = v1[a2][2] - v1[a1][2];
	  temp[j] += dot(const_vect[j], diff);
	}
	project2(n_const, const_pairs, const_dist, const_vect,
		 pdata, GDOT, m, temp, v1, atoms);
	for (j = 0; j < atoms; j++) {
	  xp[j][0] += x[j][0];
	  xp[j][1] += x[j][1];
	  xp[j][2] += x[j][2];
	}
	mult_by_h_plus_one(xp, xph, atoms, m, const_pairs, pdata, n_const);
      }
      dalpha = (pressure-b_press)/(3.*volume*b_mass)
	             -(*b_alpha)*((*t_xi)+3.*(*b_alpha));
      factor1 = (*b_alpha) + dth*dalpha;
      factor2 = dth*(*b_alpha)*(*b_alpha);
    }
    else {
      factor1 = factor2 = 0.;
      dalpha = 0.; /* unused, initialize just to make gcc happy */
    }
    for (j = 0; j < atoms; j++)
      if (!fix[j]) {
	vector3 dv;
	dv[0] = -dth*(f[j][0]/m[j]+(*t_xi)*v[j][0]);
	dv[1] = -dth*(f[j][1]/m[j]+(*t_xi)*v[j][1]);
	dv[2] = -dth*(f[j][2]/m[j]+(*t_xi)*v[j][2]);
	x[j][0] += delta_t*(v[j][0] + dv[0]
			    + factor1*xp[j][0] + factor2*xph[j][0]);
	x[j][1] += delta_t*(v[j][1] + dv[1]
			    + factor1*xp[j][1] + factor2*xph[j][1]);
	x[j][2] += delta_t*(v[j][2] + dv[2]
			    + factor1*xp[j][2] + factor2*xph[j][2]);
	if (barostat && n_const > 0) {
	  x[j][0] += delta_t*dth*v1[j][0];
	  x[j][1] += delta_t*dth*v1[j][1];
	  x[j][2] += delta_t*dth*v1[j][2];
	}
	v[j][0] += dv[0]-dth*(*b_alpha)*vh[j][0];
	v[j][1] += dv[1]-dth*(*b_alpha)*vh[j][1];
	v[j][2] += dv[2]-dth*(*b_alpha)*vh[j][2];
      }
    if (barostat) {
      volume = universe_spec->volume_function(1.+delta_t*(factor1+factor2),
					      universe_spec->geometry_data);
      (*b_alpha) += dth*dalpha;
    }
    if (thermostat) {
      dxi = (2.*k_energy - t_energy +
	     9.*b_mass*volume*volume*(*b_alpha)*(*b_alpha))/t_mass;
      (*t_xi) += dth*dxi;
      (*t_lns) += delta_t*(*t_xi);
    }

    /* Constraints: SHAKE */
    if (n_const > 0) {
      memcpy(xold, x, atoms*sizeof(vector3));
      for (j = 0; j < n_const_blocks; j++)
	shake(const_pairs, const_blocks[j], const_blocks[j+1],
	      x, m, const_vect, const_dist, universe_spec->distance_function,
	      universe_spec->geometry_data);
      for (j = 0; j < atoms; j++)
	if (!fix[j]) {
	  v[j][0] += (x[j][0]-xold[j][0])/delta_t;
	  v[j][1] += (x[j][1]-xold[j][1])/delta_t;
	  v[j][2] += (x[j][2]-xold[j][2])/delta_t;
	}
      for (j = 0; j < n_const; j++)
	universe_spec->distance_function(const_vect[j],
					 x[const_pairs[2*j]],
					 x[const_pairs[2*j+1]],
					 universe_spec->geometry_data);
    }

    /* Coordinate correction (for periodic universes etc.) */
    universe_spec->correction_function(x, atoms, universe_spec->geometry_data);

    /* Mid-step energy evaluation */
    PyUniverseSpec_StateLock(universe_spec, -2);
    PyUniverseSpec_StateLock(universe_spec, 1);
    (*evaluator->eval_func)(evaluator, &p_energy, configuration, 1);
    PyUniverseSpec_StateLock(universe_spec, 2);
    if (p_energy.error) {
      PyEval_RestoreThread(evaluator->tstate_save);
      goto error;
    }
    PyUniverseSpec_StateLock(universe_spec, -1);

    /* Second part of integration step */
    for (j = 0; j < atoms; j++)
      if (!fix[j]) {
	double factor = dth/m[j];
	v[j][0] -= factor*f[j][0];
	v[j][1] -= factor*f[j][1];
	v[j][2] -= factor*f[j][2];
      }

    /* Constraints: constraint forces */
    if (n_const > 0) {
      project(n_const, const_pairs, const_dist, const_vect,
	      pdata, MU1, m, v, v1, atoms);
      if (barostat) {
	mult_by_h_plus_one(v, vh, atoms, m, const_pairs, pdata, n_const);
	project(n_const, const_pairs, const_dist, const_vect,
		pdata, MU2, m, vh, v2, atoms);
      }
      const_virial1 = const_virial2 = 0.;
      for (j = 0; j < n_const; j++) {
	const_virial1 += pdata[j].multiplier[MU1]*const_dist[j];
	const_virial2 += pdata[j].multiplier[MU2]*const_dist[j];
      }
    }

    /* Iterative determination of xi and alpha */
    if (thermostat || barostat) {
      double tol = 1.e-15;
      int niter = 1000;
      double alpha, xi;
      double ke1, ke2, ke3;

      ke1 = ke2 = ke3 = 0.;
      if (barostat && n_const > 0)
	for (j = 0; j < atoms; j++) {
	  vector3 va, vb;
	  va[0] = v[j][0]-v1[j][0];
	  va[1] = v[j][1]-v1[j][1];
	  va[2] = v[j][2]-v1[j][2];
	  vb[0] = vh[j][0]-v2[j][0];
	  vb[1] = vh[j][1]-v2[j][1];
	  vb[2] = vh[j][2]-v2[j][2];
	  ke1 += m[j]*dot(va, va);
	  ke2 += m[j]*dot(vb, vb);
	  ke3 += m[j]*dot(va, vb);
	}
      else {
	for (j = 0; j < atoms; j++)
	  ke1 += m[j]*dot(v[j], v[j]);
	ke2 = ke3 = ke1;
      }

      if (barostat)
	(*b_alpha) += dth*(p_energy.virial/(3.*volume)-b_press)
	                            / (3.*volume*b_mass);
      if (thermostat)
	(*t_xi) -= dth*t_energy/t_mass;
      alpha = (*b_alpha);
      xi = (*t_xi);
      while (1) {
	double ke = (1.-dth*xi)*(1.-dth*xi)*ke1 + dth*dth*alpha*alpha*ke2
	              - 2*(1.-dth*xi)*dth*alpha*ke3;
	double kea = delta_t*(dth*alpha*ke2-(1.-dth*xi)*ke3);
	double kex = delta_t*((dth*xi-1.)*ke1+dth*alpha*ke3);
	double cv = (xi-1./dth)*const_virial1 + alpha*const_virial2;
	double vf = 9.*volume*volume*b_mass;
	double h1 = alpha - (*b_alpha)
                          - dth*((ke+cv)/vf - (xi+3.*alpha)*alpha);
	double h1a = -1+dth*((kea+const_virial2)/vf-xi-6.*alpha);
	double h1x = dth*((kex+const_virial1)/vf-alpha);
	double h2 = xi - (*t_xi) - dth*(ke+vf*alpha*alpha)/t_mass;
	double h2a = dth*(kea+2.*vf*alpha)/t_mass;
	double h2x = -1+dth*kex/t_mass;
	double da = 0.;
	double dx = 0.;
	double converged;
	if (thermostat && barostat) {
	  da = (h2x*h1-h1x*h2)/(h2x*h1a-h1x*h2a);
	  dx = (h2-h2a*da)/h2x;
	}
	else if (thermostat)
	  dx = h2/h2x;
	else
	  da = h1/h1a;
	xi += dx;
	alpha += da;
	converged = 1;
	if (thermostat)
	  converged = converged && fabs(dx) < tol*fabs(xi);
	if (barostat)
	  converged = converged && fabs(da) < tol*fabs(alpha);
	if (converged) {
#if DEBUG
	  printf("%d iterations, %lf/%lf\n", 1000-niter, h1, h2);
#endif
	  break;
	}
	if (niter-- == 0) {
#if DEBUG
	  printf("No convergence\n");
#endif
	  break;
	}
      }
      if (thermostat) (*t_xi) = xi;
      if (barostat) (*b_alpha) = alpha;
    }
    /* Final velocity calculation */
    for (j = 0; j < atoms; j++)
      if (!fix[j]) {
	if (n_const > 0) {
	  v[j][0] -= v1[j][0];
	  v[j][1] -= v1[j][1];
	  v[j][2] -= v1[j][2];
	  if (barostat) {
	    vh[j][0] -= v2[j][0];
	    vh[j][1] -= v2[j][1];
	    vh[j][2] -= v2[j][2];
	  }
	}
	v[j][0] = (1.-dth*(*t_xi))*v[j][0] - dth*(*b_alpha)*vh[j][0];
	v[j][1] = (1.-dth*(*t_xi))*v[j][1] - dth*(*b_alpha)*vh[j][1];
	v[j][2] = (1.-dth*(*t_xi))*v[j][2] - dth*(*b_alpha)*vh[j][2];
      }

    /* The End - next time step! */
    time += delta_t;
  }
  /** End of main integration loop **/

  /* Final thermodynamic property evaluation */
  k_energy = 0.;
  for (j = 0; j < atoms; j++)
    if (!fix[j])
      k_energy += m[j]*dot(v[j], v[j]);
  k_energy *= 0.5;
  n_energy = 0.5*(*t_xi)*(*t_xi)*t_mass*exp(2*(*t_lns))
                + t_energy*(*t_lns);
  a_energy = 4.5*b_mass*volume*volume*(*b_alpha)*(*b_alpha) + volume*b_press;
  temperature = 2.*k_energy*temperature_factor/df;
  pressure = (2.*k_energy+p_energy.virial
	      + ((*t_xi)-1./dth)*const_virial1
	      + (*b_alpha)*const_virial2)/(3.*volume);

  /* Final trajectory and log output */
  if (PyTrajectory_Output(output, i, data_descriptors,
			  &evaluator->tstate_save) == -1) {
    PyUniverseSpec_StateLock(universe_spec, -2);
    PyEval_RestoreThread(evaluator->tstate_save);
    goto error;
  }

  /* Cleanup */
  PyUniverseSpec_StateLock(universe_spec, -2);
  PyEval_RestoreThread(evaluator->tstate_save);
  PyTrajectory_OutputFinish(output, i, 0, 1, data_descriptors);
  free(scratch);
  free(pdata);
  free(data_descriptors);
  Py_DECREF(gradients);
  Py_INCREF(Py_None);
  return Py_None;

  /* Error return */
error:
  PyTrajectory_OutputFinish(output, i, 1, 1, data_descriptors);
error2:
  free(scratch);
  free(pdata);
  free(data_descriptors);
  Py_DECREF(gradients);
  return NULL;
}




/* Runge-Kutta-Merson integrator - modified little by little*/

static PyObject *
integrateRKM_EVOL(PyObject *dummy, PyObject *args)
{
  PyObject *universe;
  PyUniverseSpecObject *universe_spec;
  PyArrayObject *configuration;
  PyArrayObject *velocities;
  PyArrayObject *masses;
  PyArrayObject *fixed;
  PyArrayObject *gradients;
  PyArrayObject *constraints, *constraint_distances_squared, *c_blocks;
  PyArrayObject *t_parameters, *t_coordinates;
  PyArrayObject *b_parameters, *b_coordinates;
  PyListObject *spec_list;
  PyFFEvaluatorObject *evaluator;
  PyTrajectoryOutputSpec *output;
  PyTrajectoryVariable *data_descriptors = NULL;
  double delta_t, dth;
  int first_step, last_step;
  char *description;
  vector3 *x, *v, *f;
  double *m, *const_dist;
  long *fix, *const_pairs, *const_blocks;
  vector3 *const_vect;
  vector3 *scratch;
  vector3 *xp, *xph, *vh, *v1, *v2, *xold;
  double *temp;
  projection_data *pdata;
  double time;
  energy_data p_energy;
  double k_energy, n_energy, a_energy;
  double temperature, volume, const_virial1, const_virial2, pressure;
  double t_temp, t_tau, t_mass, t_energy, *t_xi, *t_lns;
  double b_press, b_tau, b_mass, *b_alpha;
  double dxi, dalpha, factor1, factor2;
  int atoms, n_const, n_const_blocks, df;
  int pressure_available, thermostat, barostat;
  int i, j;

  /* Parse and check arguments */
  if (!PyArg_ParseTuple(args, "OO!O!O!O!O!O!O!O!O!O!O!O!diiO!s", &universe,
			&PyArray_Type, &configuration,
			&PyArray_Type, &velocities,
			&PyArray_Type, &masses,
			&PyArray_Type, &fixed,
			&PyFFEvaluator_Type, &evaluator,
			&PyArray_Type, &constraints,
			&PyArray_Type, &constraint_distances_squared,
			&PyArray_Type, &c_blocks,
			&PyArray_Type, &t_parameters,
			&PyArray_Type, &t_coordinates,
			&PyArray_Type, &b_parameters,
			&PyArray_Type, &b_coordinates,
			&delta_t, &first_step, &last_step,
			&PyList_Type, &spec_list,
			&description))
    return NULL;
  universe_spec = (PyUniverseSpecObject *)
                   PyObject_GetAttrString(universe, "_spec");
  if (universe_spec == NULL)
    return NULL;

  /* Create gradient array */
#if defined(NUMPY)
  gradients = (PyArrayObject *)PyArray_Copy(configuration);
#else
  gradients = (PyArrayObject *)PyArray_FromDims(configuration->nd,
						configuration->dimensions,
						PyArray_DOUBLE);
#endif
  if (gradients == NULL)
    return NULL;

  /* Set some convenient variables */
  atoms = configuration->dimensions[0];
  n_const = constraints->dimensions[0];
  n_const_blocks = c_blocks->dimensions[0]-1;
  thermostat = (t_parameters->dimensions[0] > 0);
  barostat = (b_parameters->dimensions[0] > 0);
  dth = 0.5*delta_t;
  x = (vector3 *)configuration->data;
  v = (vector3 *)velocities->data;
  f = (vector3 *)gradients->data;
  m = (double *)masses->data;
  fix = (long *)fixed->data;
  const_pairs = (long *)constraints->data;
  const_dist = (double *)constraint_distances_squared->data;
  const_blocks = (long *)c_blocks->data;

  printf(" -- integrateRKM_EVOL BEGIN\n");

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  /* Check MMTK  &&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  char prefix[11];
  char preid[3];
  char sufix[5];
  prefix[0] = 'p'; prefix[1] = 'd'; prefix[2] = 'b';
  prefix[3] = 's'; prefix[4] = '/';
  prefix[5] = 'm'; prefix[6] = 'm'; prefix[7] = 't';
  prefix[8] = 'k'; prefix[9] = '_'; prefix[10] = '\0';
  preid[0] = preid[1] = preid[2] = '\0';
  sufix[0] = '.'; sufix[1] = 'p'; sufix[2] = 'd';
  sufix[3] = 'b'; sufix[4] = '\0';

  int EVOL_i = 0;
  char elem;
  printf("integrate_EVOL: Structure got from Python\n");
  for(EVOL_i=0; EVOL_i<atoms; EVOL_i++){
    if(((m[EVOL_i] - 1) < 0.2) && ((m[EVOL_i] - 1) > -0.2) )
      elem = 'H';
    else if(((m[EVOL_i] - 12) < 0.2) && ((m[EVOL_i] - 12) > -0.2) )
      elem = 'C';
    else if(((m[EVOL_i] - 14) < 0.2) && ((m[EVOL_i] - 14) > -0.2) )
      elem = 'N';
    else if(((m[EVOL_i] - 16) < 0.2) && ((m[EVOL_i] - 16) > -0.2) )
      elem = 'O';
    else if(((m[EVOL_i] - 32) < 0.2) && ((m[EVOL_i] - 32) > -0.2) )
      elem = 'S';
    else if(((m[EVOL_i] - 30) < 0.2) && ((m[EVOL_i] - 30) > -0.2) )
      elem = 'P';
    else elem = 'X';

    printf("HETATM %4d %c%-3d UNK     1    %8.3f%8.3f%8.3f  1.00  0.00\n",
      EVOL_i, elem, EVOL_i, x[EVOL_i][0] * 10, x[EVOL_i][1] * 10, x[EVOL_i][2] * 10);
  } 




  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  /* Convert description to array of floats */
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  printf("integrateRKM_EVOL description from descorder: %s\n", description);
  int order[atoms];
  int ci = 0;
  int oi = 0;
  int bi = 0;
  char buff[7];
  for(bi=0; bi<7; bi++){
    buff[bi] = '\0';
  }
  bi = 0;
  while(description[ci] != '\0'){
    if(description[ci] == '|'){
      order[oi++] = atof(buff);
      for(bi=0; bi<7; bi++){
        buff[bi] = '\0';
      }
      bi = 0;
      ci++;
    }
    else{
      buff[bi] = description[ci];
      bi++;
      ci++;
    }
  }
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

  /* Debug a lot */
  printf(" -- integrateRKM_EVOL order from descorder: ");
  for(oi=0; oi<atoms; oi++)
    printf("%d|", order[oi]);
  printf("\n");
  
  printf(" -- integrateRKM_EVOL delta_t = %.5f; first_step = %d; last_step = %d; ",\
              delta_t,        first_step,      last_step);
  printf(" -- integrateRKM_EVOL atoms = %d; n_const = %d; n_const_blocks = %d; thermostat  = %d; barostat = %d; dth = %.5f ", \
              atoms,      n_const,      n_const_blocks,      thermostat,       barostat,      dth);
  printf("\n -- integrateRKM_EVOL configuration: ");
  for(EVOL_i=0; EVOL_i<atoms; EVOL_i++)
    printf("[%.2f %.2f %.2f] ", x[EVOL_i][0], x[EVOL_i][1], x[EVOL_i][2]);
  printf("\n -- integrateRKM_EVOL velocities: ");
  for(EVOL_i=0; EVOL_i<atoms; EVOL_i++)
    printf("[%.2f %.2f %.2f] ", v[EVOL_i][0], v[EVOL_i][1], v[EVOL_i][2]);
  printf("\n -- integrateRKM_EVOL gradients: ");
  for(EVOL_i=0; EVOL_i<atoms; EVOL_i++)
    printf("[%.2f %.2f %.2f] ", f[EVOL_i][0], f[EVOL_i][1], f[EVOL_i][2]);
  printf("\n -- integrateRKM_EVOL masses: ");
  for(EVOL_i=0; EVOL_i<atoms; EVOL_i++){printf("%.1f ", m[EVOL_i]);}
  printf("\n -- integrateRKM_EVOL fix: ");
  for(EVOL_i=0; EVOL_i<atoms; EVOL_i++){printf("%ld ", fix[EVOL_i]);}
  printf("\n -- integrateRKM_EVOL const pairs: ");
  for(EVOL_i=0; EVOL_i<n_const; EVOL_i++){printf("%ld ", const_pairs[EVOL_i]);}
  printf("\n -- integrateRKM_EVOL const dists: ");
  for(EVOL_i=0; EVOL_i<n_const; EVOL_i++){printf("%.2f ", const_dist[EVOL_i]);}
  printf("\n -- integrateRKM_EVOL const blocks: ");
  for(EVOL_i=0; EVOL_i<n_const_blocks; EVOL_i++){printf("%ld ", const_blocks[EVOL_i]);}


  /* Calculate number of degrees of freedom */
  df = 3*atoms;
  for (j = 0; j < atoms; j++) {
    if (fix[j])
      df -= 3;
  }
  for (j = 0; j < n_const; j++) {
    if (fix[const_pairs[2*j]] || fix[const_pairs[2*j+1]]) {
      PyErr_SetString(PyExc_ValueError,
		      "distance constraint on a fixed atom");
      return NULL;
    }
    df--;
  }

  /* Thermostat and barostat information.
     Note: the coordinates are guaranteed to be zero when there is
     no thermostat and/or barostat. */
  if (thermostat) {
    t_temp = *(double *)t_parameters->data;
    t_tau = *(((double *)t_parameters->data)+1);
    t_energy = df*t_temp*kB;
    t_mass = t_energy*t_tau*t_tau;
  }
  else {
    t_mass = 0.;
    t_energy = 0.;
    t_temp = 0.; /* unused, initialize just to make gcc happy */
    t_mass = 0;  /* unused, initialize just to make gcc happy */
  }
  t_xi = (double *)t_coordinates->data;
  t_lns = t_xi + 1;
  if (barostat) {
    b_press = *(double *)b_parameters->data;
    b_tau = *(((double *)b_parameters->data)+1);
  }
  else {
    b_mass = 0.;
    b_press = 0.;
    b_tau = 0.; /* unused, initialize just to make gcc happy */
  }
  b_alpha = (double *)b_coordinates->data;

  /* Allocate arrays for temporary data */
  pdata = NULL;
  if (n_const > 0) {
    scratch = (vector3 *)malloc((5*atoms+n_const)*sizeof(vector3));
    if (scratch == NULL) {
      PyErr_NoMemory();
      goto error2;
    }
    v1 = scratch;
    v2 = v1 + atoms;
    xp = v2 + atoms;
    xph = xp + atoms;
    vh = xph + atoms;
    const_vect = vh + atoms;
    xold = v1;
    temp = (double *)v2;
    pdata = (projection_data *)malloc(n_const*sizeof(projection_data));
    if (pdata == NULL) {
      PyErr_NoMemory();
      goto error2;
    }
  }
  else {
    scratch = NULL;
    v1 = v2 = NULL;
    xp = xph = x;
    vh = v;
    xold = NULL;
    temp = NULL;
    const_vect = NULL;  /* unused, initialize just to make gcc happy */
  }

  /* Enforce constraints and initialize constraint data */
  Py_BEGIN_ALLOW_THREADS;
  PyUniverseSpec_StateLock(universe_spec, -1);
  universe_spec->correction_function(x, atoms, universe_spec->geometry_data);
  if (n_const > 0) {
    for (j = 0; j < n_const; j++) {
      universe_spec->distance_function(const_vect[j],
				       x[const_pairs[2*j]],
				       x[const_pairs[2*j+1]],
				       universe_spec->geometry_data);
      pdata[j].multiplier[MU1] = 0.;
      pdata[j].multiplier[MU2] = 0.;
      pdata[j].multiplier[G] = 0.;
      pdata[j].multiplier[GDOT] = 0.;
      pdata[j].diag = const_dist[j] *
	(1./m[const_pairs[2*j]] + 1./m[const_pairs[2*j+1]]);
    }
    for (j = 0; j < n_const_blocks; j++)
      shake(const_pairs, const_blocks[j], const_blocks[j+1],
	    x, m, const_vect, const_dist, universe_spec->distance_function,
	    universe_spec->geometry_data);
    for (j = 0; j < n_const; j++)
      universe_spec->distance_function(const_vect[j],
				       x[const_pairs[2*j]],
				       x[const_pairs[2*j+1]],
				       universe_spec->geometry_data);
    if (barostat) {
      project(n_const, const_pairs, const_dist, const_vect,
	      pdata, G, m, NULL, xp, atoms);
      mult_by_h_plus_one(v, vh, atoms, m, const_pairs, pdata, n_const);
      project(n_const, const_pairs, const_dist, const_vect,
	      pdata, MU2, m, vh, v2, atoms);
      for (j = 0; j < n_const; j++) {
	int a1 = const_pairs[2*j];
	int a2 = const_pairs[2*j+1];
	vector3 diff1, diff2;
	diff1[0] = (*b_alpha)*(const_vect[j][0] + xp[a2][0]-xp[a1][0])
	           + v[a2][0] - v[a1][0];
	diff1[1] = (*b_alpha)*(const_vect[j][1] + xp[a2][1]-xp[a1][1])
	           + v[a2][1] - v[a1][1];
	diff1[2] = (*b_alpha)*(const_vect[j][2] + xp[a2][2]-xp[a1][2])
	           + v[a2][2] - v[a1][2];
	diff2[0] = v[a2][0] - v[a1][0];
	diff2[1] = v[a2][1] - v[a1][1];
	diff2[2] = v[a2][2] - v[a1][2];
	temp[j] = dot(diff1, diff2);
	diff1[0] = f[a1][0]/m[a1] - f[a2][0]/m[a1];
	diff1[1] = f[a1][1]/m[a1] - f[a2][1]/m[a1];
	diff1[2] = f[a1][2]/m[a1] - f[a2][2]/m[a1];
	temp[j] += dot(const_vect[j], diff1);
	temp[j] *= dth/(dth*(*t_xi)-1.);
      }
      project2(n_const, const_pairs, const_dist, const_vect,
	       pdata, MU1, m, temp, v1, atoms);
    }
    const_virial1 = const_virial2 = 0.;
    for (j = 0; j < n_const; j++) {
      const_virial1 += pdata[j].multiplier[MU1]*const_dist[j];
      const_virial2 += pdata[j].multiplier[MU2]*const_dist[j];
    }
  }
  else
    const_virial1 = const_virial2 = 0.;
  PyUniverseSpec_StateLock(universe_spec, -2);
  Py_END_ALLOW_THREADS;

  /* Initial force calculation */
  p_energy.gradients = (PyObject *)gradients;
  p_energy.gradient_fn = NULL;
  p_energy.force_constants = NULL;
  p_energy.fc_fn = NULL;
#ifdef WITH_THREAD
  evaluator->tstate_save = PyEval_SaveThread();
#endif
  PyUniverseSpec_StateLock(universe_spec, 1);
  (*evaluator->eval_func)(evaluator, &p_energy, configuration, 0);
  PyUniverseSpec_StateLock(universe_spec, 2);
#ifdef WITH_THREAD
  PyEval_RestoreThread(evaluator->tstate_save);
#endif
  if (p_energy.error)
    goto error2;

  /* Check if pressure can be calculated */
  Py_BEGIN_ALLOW_THREADS;
  PyUniverseSpec_StateLock(universe_spec, 1);
  volume = universe_spec->volume_function(1., universe_spec->geometry_data);
  PyUniverseSpec_StateLock(universe_spec, 2);
  Py_END_ALLOW_THREADS;
  pressure_available = volume > 0. && p_energy.virial_available;
  if (thermostat && barostat)
    b_mass = df*t_temp*kB*b_tau*b_tau/(volume*volume);

  /* Initialize output */
  data_descriptors =
     get_data_descriptors(9 + pressure_available
			  + 2*thermostat
			  + 3*barostat
			  + (universe_spec->geometry_data_length > 0),
			  universe_spec,
			  configuration, velocities,
			  gradients, masses, &df,
			  &time, &p_energy.energy, &k_energy,
			  thermostat ? &n_energy : NULL,
			  barostat ? &a_energy : NULL,
			  &temperature,
			  thermostat ? t_xi : NULL,
			  pressure_available ? &pressure:NULL,
			  barostat ? &volume : NULL,
			  barostat ? b_alpha : NULL,
			  (universe_spec->geometry_data_length > 0) ?
			              universe_spec->geometry_data : NULL);
  if (data_descriptors == NULL)
    goto error2;
  output = PyTrajectory_OutputSpecification(universe, spec_list,
					    description,
					    data_descriptors);
  if (output == NULL)
    goto error2;

  evaluator->tstate_save = PyEval_SaveThread();
  
  /* Get write access for the integration, switching to
     read access only during energy evaluation */
  PyUniverseSpec_StateLock(universe_spec, -1);

  /** Main integration loop **/
  time  = first_step*delta_t;
  for (i = first_step; i < last_step; i++) {

    printf("integrateRKM BEGIN - step  %d\n", i+1);

    /* Calculation of thermodynamic properties */
    k_energy = 0.;
    for (j = 0; j < atoms; j++)
      if (!fix[j])
	k_energy += m[j]*dot(v[j], v[j]);
    k_energy *= 0.5;
    n_energy = 0.5*(*t_xi)*(*t_xi)*t_mass*exp(2*(*t_lns))
                 + t_energy*(*t_lns);
    a_energy = 4.5*b_mass*volume*volume*(*b_alpha)*(*b_alpha) + volume*b_press;
    temperature = 2.*k_energy*temperature_factor/df;
    pressure = (2.*k_energy+p_energy.virial
		+ ((*t_xi)-1./dth)*const_virial1
		+ (*b_alpha)*const_virial2)/(3.*volume);
    if (barostat && !thermostat)
      b_mass = 0.5*k_energy*b_tau*b_tau/(volume*volume);
    /* Trajectory and log output */
    if (PyTrajectory_Output(output, i, data_descriptors,
			    &evaluator->tstate_save) == -1) {
      PyUniverseSpec_StateLock(universe_spec, -2);
      PyEval_RestoreThread(evaluator->tstate_save);
      goto error;
    }

    /* First part of integration step */
    if (barostat) {
      if (n_const > 0) {
	project(n_const, const_pairs, const_dist, const_vect,
		pdata, G, m, NULL, xp, atoms);
	for (j = 0; j < atoms; j++)
	  v1[j][0] = v1[j][1] = v1[j][2] = 0.;
	for (j = 0; j < n_const; j++) {
	  int a1 = const_pairs[2*j];
	  int a2 = const_pairs[2*j+1];
	  vector3 diff1, diff2;
	  diff1[0] = (*b_alpha)*(const_vect[j][0] + xp[a2][0]-xp[a1][0])
	              + v[a2][0] - v[a1][0];
	  diff1[1] = (*b_alpha)*(const_vect[j][1] + xp[a2][1]-xp[a1][1])
	              + v[a2][1] - v[a1][1];
	  diff1[2] = (*b_alpha)*(const_vect[j][2] + xp[a2][2]-xp[a1][2])
	              + v[a2][2] - v[a1][2];
	  diff2[0] = xp[a2][0] - xp[a1][0];
	  diff2[1] = xp[a2][1] - xp[a1][1];
	  diff2[2] = xp[a2][2] - xp[a1][2];
	  temp[j] = dot(diff1, diff2);
	  v1[a2][0] += diff1[0]*pdata[j].multiplier[G]/m[a2];
	  v1[a2][1] += diff1[1]*pdata[j].multiplier[G]/m[a2];
	  v1[a2][2] += diff1[2]*pdata[j].multiplier[G]/m[a2];
	  v1[a1][0] -= diff1[0]*pdata[j].multiplier[G]/m[a1];
	  v1[a1][1] -= diff1[1]*pdata[j].multiplier[G]/m[a1];
	  v1[a1][2] -= diff1[2]*pdata[j].multiplier[G]/m[a1];
	}
	for (j = 0; j < n_const; j++) {
	  int a1 = const_pairs[2*j];
	  int a2 = const_pairs[2*j+1];
	  vector3 diff;
	  diff[0] = v1[a2][0] - v1[a1][0];
	  diff[1] = v1[a2][1] - v1[a1][1];
	  diff[2] = v1[a2][2] - v1[a1][2];
	  temp[j] += dot(const_vect[j], diff);
	}
	project2(n_const, const_pairs, const_dist, const_vect,
		 pdata, GDOT, m, temp, v1, atoms);
	for (j = 0; j < atoms; j++) {
	  xp[j][0] += x[j][0];
	  xp[j][1] += x[j][1];
	  xp[j][2] += x[j][2];
	}
	mult_by_h_plus_one(xp, xph, atoms, m, const_pairs, pdata, n_const);
      }
      dalpha = (pressure-b_press)/(3.*volume*b_mass)
	             -(*b_alpha)*((*t_xi)+3.*(*b_alpha));
      factor1 = (*b_alpha) + dth*dalpha;
      factor2 = dth*(*b_alpha)*(*b_alpha);
    }
    else {
      factor1 = factor2 = 0.;
      dalpha = 0.; /* unused, initialize just to make gcc happy */
    }



    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
    /* Client stuff 1 - EVOL &&&&&&&&&&&&&&&&&&&&&&&&&&&*/
    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
#ifndef TARGET_TYPE
#define TARGET_TYPE     double
#endif
      float myxold[atoms][3];
      float xfromshm[atoms][3];
  
    unsigned long int SHMSZ = (
      2*sizeof(TARGET_TYPE) + // round and flag
      3*atoms*sizeof(TARGET_TYPE) + // Positions
      3*atoms*sizeof(TARGET_TYPE) + // Velocities
      3*atoms*sizeof(TARGET_TYPE) + // indexMap
      3*atoms*sizeof(TARGET_TYPE) + // Gradients
      3*sizeof(TARGET_TYPE)
    );
    printf("Client round %d - Shared memory size = %ld\n", i, SHMSZ); fflush(stdout);
  
    float CLIENT_FLAG = 1.0;
    float SERVER_FLAG = 0.0;
    //long int CNT_START = 1945;
    //long int CNT_STOP = 1947;
  
    int shmid;
    unsigned int cli=0;
    key_t key;
    TARGET_TYPE *shm;
    float myround = 0;
  
    ////// Get the segment
    key = 5678;
    if ((shmid = shmget(key, SHMSZ, IPC_CREAT | 0666)) < 0) {
      perror("shmget");
      exit(1);
    }
    printf("Got %d values\n", (SHMSZ)/sizeof(TARGET_TYPE));
  
    ////// Map the segment
    if ((shm = (TARGET_TYPE *)shmat(shmid, NULL, 0)) == (TARGET_TYPE *) -1) {
      perror("shmat");
      exit(1);
    }
  
  
    while(1){ ////// Wait until server is done
      usleep(1000);
      if(shm[1] == SERVER_FLAG){
        // Write coords and vels to memory
        printf("CLIENT: CALL SERVER step %d", i);
        int a;
        cli = 2;
        for (a=0; a<atoms; a++){
          shm[cli++] = x[a][0];
          shm[cli++] = x[a][1];
          shm[cli++] = x[a][2];
        }
        for (a=0; a<atoms; a++){
          shm[cli++] = v[a][0];
          shm[cli++] = v[a][1];
          shm[cli++] = v[a][2];
        }
        for (a=0; a<atoms; a++){
          cli += 2;
          shm[cli++] = order[a];
        }
        for (a=0; a<atoms; a++){
          shm[cli++] = f[a][0];
          shm[cli++] = f[a][1];
          shm[cli++] = f[a][2];
        }
  
        printf("shm index cli went up to = %d\n", cli); fflush(stdout);
        int j;
  
        printf("CLIENT SENT coords: [");
        j=0;
        for(cli=2; cli<2+3*atoms; cli++){
          printf("%.2f ", shm[cli]);
          if(j%3==0){printf("][");}
          j++;
        }
        printf("\n"); fflush(stdout);
  
        printf("CLIENT SENT vels: [");
        j=0;
        for(cli=2+3*atoms; cli<2+3*atoms+3*atoms; cli++){
          printf("%.2f ", shm[cli]);
          if(j%3==0){printf("][");}
          j++;
        }
        printf("\n"); fflush(stdout);
  
        printf("CLIENT SENT grads: [");
        j=0;
        for(cli=2+3*atoms+3*atoms+3*atoms; cli<2+3*atoms+3*atoms+3*atoms+3*atoms; cli++){
          printf("%.2f ", shm[cli]);
          if(j%3==0){printf("][");}
          j++;
        }
        printf("\n"); fflush(stdout);
  
  
        ////// Change flag
        shm[1] = CLIENT_FLAG;
        printf("CLIENT: CLIENT_FLAG ON\n");
        break;
      
      }
    }
    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
    /* End of client stuff 1 - EVOL &&&&&&&&&&&&&&&&&&&&*/
    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  
    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
    /* Client stuff 2 - EVOL &&&&&&&&&&&&&&&&&&&&&&&&&&&*/
    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
    int waiting_rounds = 0;
    int PATIENCE = 10000;
    while(1){ ////// Wait until server is done
      usleep(1000);
      ////// If server is done read memory
      if(shm[1] == SERVER_FLAG){
        printf("CLIENT: SERVER IS DONE\n");
  #ifdef DEBUGALOT
        for(cli=2; cli<(SHMSZ)/sizeof(TARGET_TYPE); cli++){
          printf("CLIENT RECEIEVED i=%d: %.2f\n", cli, shm[cli]);
        }
  #endif
  
        /* Read positions and velocities from memory*/
        cli = 2;
        int a;
  
        for (a=0; a<atoms; a++){
          myxold[a][0] = x[a][0];
          myxold[a][1] = x[a][1];
          myxold[a][2] = x[a][2];
        }
  
        for (a=0; a<atoms; a++){
          xfromshm[a][0] = shm[cli++];
          xfromshm[a][1] = shm[cli++];
          xfromshm[a][2] = shm[cli++];
        }
  
  
  /*
        for (a=0; a<atoms; a++){
          v[a][0] = shm[cli++];
          v[a][1] = shm[cli++];
          v[a][2] = shm[cli++];
        }
        for (a=0; a<atoms; a++){ // skip indexMap
          cli += 3;
        }
        for (a=0; a<atoms; a++){
          f[a][0] = shm[cli++];
          f[a][1] = shm[cli++];
          f[a][2] = shm[cli++];
        }
  */
        printf("CLIENT RECEIEVED coords:");
        for(a=0; a<atoms; a++){
          printf("[%.2f %.2f %.2f] ", x[a][0], x[a][1], x[a][2]);
        }
        printf("\n");
  
        printf("CLIENT RECEIEVED vels:");
        for(a=0; a<atoms; a++){
          printf("[%.5f %.5f %.5f] ", v[a][0], v[a][1], v[a][2]);
        }
        printf("\n");
  
        printf("CLIENT RECEIEVED grads:");
        for(a=0; a<atoms; a++){
          printf("[%.5f %.5f %.5f] ", f[a][0], f[a][1], f[a][2]);
        }
        printf("\n");
  
        //shmdt(shm);
        break;
      }
      if(waiting_rounds >= PATIENCE){
        printf("Client waited more then %d seconds. Exiting...\n", PATIENCE/1000);
        exit(1);
      }
      waiting_rounds++;
    }
  
    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
    /* End of client stuff 2 - EVOL &&&&&&&&&&&&&&&&&&&&*/
    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/


    /*&&&&&&&&&&&&&&&&&&&&&&&*/
    cli = 2;
    /*&&&&&&&&&&&&&&&&&&&&&&&*/

    for (j = 0; j < atoms; j++){
      if (!fix[j]) {
	vector3 dv;
	dv[0] = -dth*(f[j][0]/m[j]+(*t_xi)*v[j][0]);
	dv[1] = -dth*(f[j][1]/m[j]+(*t_xi)*v[j][1]);
	dv[2] = -dth*(f[j][2]/m[j]+(*t_xi)*v[j][2]);


/*
	x[j][0] += delta_t*(v[j][0] + dv[0]
			    + factor1*xp[j][0] + factor2*xph[j][0]);
	x[j][1] += delta_t*(v[j][1] + dv[1]
			    + factor1*xp[j][1] + factor2*xph[j][1]);
	x[j][2] += delta_t*(v[j][2] + dv[2]
			    + factor1*xp[j][2] + factor2*xph[j][2]);
*/
        /*&&&&&&&&&&&&&&&&&&&&&&&*/
        x[j][0] = xfromshm[j][0];
        x[j][1] = xfromshm[j][1];
        x[j][2] = xfromshm[j][2];
        /*&&&&&&&&&&&&&&&&&&&&&&&*/


	if (barostat && n_const > 0) {
	  x[j][0] += delta_t*dth*v1[j][0];
	  x[j][1] += delta_t*dth*v1[j][1];
	  x[j][2] += delta_t*dth*v1[j][2];
	}
	v[j][0] += dv[0]-dth*(*b_alpha)*vh[j][0];
	v[j][1] += dv[1]-dth*(*b_alpha)*vh[j][1];
	v[j][2] += dv[2]-dth*(*b_alpha)*vh[j][2];
      }
    }


        int a;
        int cti;
        int ct;
        float tres;
        for(cti=0; cti<atoms; cti++){
          ct = cti;
          for (a=0; a<atoms; a++){
            tres = fabs(myxold[ct][0]-x[a][0]) + fabs(myxold[ct][1]-x[a][1]) + fabs(myxold[ct][2]-x[a][2]);
            if(tres < 0.18){
              printf("integrateRKM_EVOL: myxold %d - x: %3d [%9.4f%9.4f%9.4f] [%9.4f%9.4f%9.4f] = [%9.4f%9.4f%9.4f]\n",
                ct,a, myxold[ct][0],         myxold[ct][1],         myxold[ct][2], x[a][0], x[a][1], x[a][2],
                      myxold[ct][0]-x[a][0], myxold[ct][1]-x[a][1], myxold[ct][2]-x[a][2]);
            }
          }
          printf("-----------------------------------------------------------------------------\n");
        }

    /*&&&&&&&& Detach the memory &&&&&&&&&&&&&&&*/
    shmdt(shm);
    /*&&&&&&&&&&&&&&&&&&&&&&&*/

    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
    /* Write structure to file - EVOL &&&&&&&&&&&&&&&&&&*/
    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
    char filename1[18];
    char filename2[18];
    FILE *fpo;
  
    if(i<10){
      preid[0] = '0'; preid[1] = '0'; preid[2] = '\0';
    }
    else if((i>=10) && (i<100)){
      preid[0] = '0'; preid[1] = '\0'; preid[2] = '\0';
    }
  
    sprintf(filename1, "%s%s%d%s", prefix, preid, i, sufix);
    fpo = fopen(filename1, "w");
    if(fpo == NULL){
      printf("File not opened...\n");
      exit(2);
    }
  
    fprintf(fpo, "REMARK integrate_EVOL: Structure after client stuff\n");
    for(EVOL_i=0; EVOL_i<atoms; EVOL_i++){
      if(((m[EVOL_i] - 1) < 0.2) && ((m[EVOL_i] - 1) > -0.2) )
        elem = 'H';
      else if(((m[EVOL_i] - 12) < 0.2) && ((m[EVOL_i] - 12) > -0.2) )
        elem = 'C';
      else if(((m[EVOL_i] - 14) < 0.2) && ((m[EVOL_i] - 14) > -0.2) )
        elem = 'N';
      else if(((m[EVOL_i] - 16) < 0.2) && ((m[EVOL_i] - 16) > -0.2) )
        elem = 'O';
      else if(((m[EVOL_i] - 32) < 0.2) && ((m[EVOL_i] - 32) > -0.2) )
        elem = 'S';
      else if(((m[EVOL_i] - 30) < 0.2) && ((m[EVOL_i] - 30) > -0.2) )
        elem = 'P';
      else elem = 'X';
  
      fprintf(fpo, "HETATM %4d %c%-3d UNK     1    %8.3f%8.3f%8.3f  1.00  0.00\n",
        EVOL_i, elem, EVOL_i, x[EVOL_i][0] * 10, x[EVOL_i][1] * 10, x[EVOL_i][2] * 10);
      printf("HETATM %4d %c%-3d UNK     1    %8.3f%8.3f%8.3f  1.00  0.00\n",
        EVOL_i, elem, EVOL_i, x[EVOL_i][0] * 10, x[EVOL_i][1] * 10, x[EVOL_i][2] * 10);
    }
  
    fclose(fpo);
  
    sprintf(filename2, "%s%s%d%s", "pdbs/ante_", preid, i, sufix);
    fpo = fopen(filename2, "w");
    if(fpo == NULL){
      printf("File not opened...\n");
      exit(2);
    }
  
    fprintf(fpo, "REMARK integrate_EVOL: Structure after client stuff\n");
    for(EVOL_i=0; EVOL_i<atoms; EVOL_i++){
      if(((m[EVOL_i] - 1) < 0.2) && ((m[EVOL_i] - 1) > -0.2) )
        elem = 'H';
      else if(((m[EVOL_i] - 12) < 0.2) && ((m[EVOL_i] - 12) > -0.2) )
        elem = 'C';
      else if(((m[EVOL_i] - 14) < 0.2) && ((m[EVOL_i] - 14) > -0.2) )
        elem = 'N';
      else if(((m[EVOL_i] - 16) < 0.2) && ((m[EVOL_i] - 16) > -0.2) )
        elem = 'O';
      else if(((m[EVOL_i] - 32) < 0.2) && ((m[EVOL_i] - 32) > -0.2) )
        elem = 'S';
      else if(((m[EVOL_i] - 30) < 0.2) && ((m[EVOL_i] - 30) > -0.2) )
        elem = 'P';
      else elem = 'X';
  
      fprintf(fpo, "HETATM %4d %c%-3d UNK     1    %8.3f%8.3f%8.3f  1.00  0.00\n",
        EVOL_i, elem, EVOL_i, myxold[EVOL_i][0] * 10, myxold[EVOL_i][1] * 10, myxold[EVOL_i][2] * 10);
    }
  
    fclose(fpo);
    printf("integrateRKM - Wrote files %s and %s\n", filename1, filename2);
    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
    /* End of writing structure - EVOL &&&&&&&&&&&&&&&&&*/
    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/


    if (barostat) {
      volume = universe_spec->volume_function(1.+delta_t*(factor1+factor2),
					      universe_spec->geometry_data);
      (*b_alpha) += dth*dalpha;
    }
    if (thermostat) {
      dxi = (2.*k_energy - t_energy +
	     9.*b_mass*volume*volume*(*b_alpha)*(*b_alpha))/t_mass;
      (*t_xi) += dth*dxi;
      (*t_lns) += delta_t*(*t_xi);
    }

    /* Constraints: SHAKE */
    if (n_const > 0) {
      memcpy(xold, x, atoms*sizeof(vector3));
      for (j = 0; j < n_const_blocks; j++)
	shake(const_pairs, const_blocks[j], const_blocks[j+1],
	      x, m, const_vect, const_dist, universe_spec->distance_function,
	      universe_spec->geometry_data);
      for (j = 0; j < atoms; j++)
	if (!fix[j]) {
	  v[j][0] += (x[j][0]-xold[j][0])/delta_t;
	  v[j][1] += (x[j][1]-xold[j][1])/delta_t;
	  v[j][2] += (x[j][2]-xold[j][2])/delta_t;
	}
      for (j = 0; j < n_const; j++)
	universe_spec->distance_function(const_vect[j],
					 x[const_pairs[2*j]],
					 x[const_pairs[2*j+1]],
					 universe_spec->geometry_data);
    }

    /* Coordinate correction (for periodic universes etc.) */
    universe_spec->correction_function(x, atoms, universe_spec->geometry_data);

    /* Mid-step energy evaluation */
    PyUniverseSpec_StateLock(universe_spec, -2);
    PyUniverseSpec_StateLock(universe_spec, 1);
    (*evaluator->eval_func)(evaluator, &p_energy, configuration, 1);
    PyUniverseSpec_StateLock(universe_spec, 2);
    if (p_energy.error) {
      PyEval_RestoreThread(evaluator->tstate_save);
      goto error;
    }
    PyUniverseSpec_StateLock(universe_spec, -1);

    /* Second part of integration step */
    for (j = 0; j < atoms; j++)
      if (!fix[j]) {
	double factor = dth/m[j];
	v[j][0] -= factor*f[j][0];
	v[j][1] -= factor*f[j][1];
	v[j][2] -= factor*f[j][2];
      }

    /* Constraints: constraint forces */
    if (n_const > 0) {
      project(n_const, const_pairs, const_dist, const_vect,
	      pdata, MU1, m, v, v1, atoms);
      if (barostat) {
	mult_by_h_plus_one(v, vh, atoms, m, const_pairs, pdata, n_const);
	project(n_const, const_pairs, const_dist, const_vect,
		pdata, MU2, m, vh, v2, atoms);
      }
      const_virial1 = const_virial2 = 0.;
      for (j = 0; j < n_const; j++) {
	const_virial1 += pdata[j].multiplier[MU1]*const_dist[j];
	const_virial2 += pdata[j].multiplier[MU2]*const_dist[j];
      }
    }

    /* Iterative determination of xi and alpha */
    if (thermostat || barostat) {
      double tol = 1.e-15;
      int niter = 1000;
      double alpha, xi;
      double ke1, ke2, ke3;

      ke1 = ke2 = ke3 = 0.;
      if (barostat && n_const > 0)
	for (j = 0; j < atoms; j++) {
	  vector3 va, vb;
	  va[0] = v[j][0]-v1[j][0];
	  va[1] = v[j][1]-v1[j][1];
	  va[2] = v[j][2]-v1[j][2];
	  vb[0] = vh[j][0]-v2[j][0];
	  vb[1] = vh[j][1]-v2[j][1];
	  vb[2] = vh[j][2]-v2[j][2];
	  ke1 += m[j]*dot(va, va);
	  ke2 += m[j]*dot(vb, vb);
	  ke3 += m[j]*dot(va, vb);
	}
      else {
	for (j = 0; j < atoms; j++)
	  ke1 += m[j]*dot(v[j], v[j]);
	ke2 = ke3 = ke1;
      }

      if (barostat)
	(*b_alpha) += dth*(p_energy.virial/(3.*volume)-b_press)
	                            / (3.*volume*b_mass);
      if (thermostat)
	(*t_xi) -= dth*t_energy/t_mass;
      alpha = (*b_alpha);
      xi = (*t_xi);
      while (1) {
	double ke = (1.-dth*xi)*(1.-dth*xi)*ke1 + dth*dth*alpha*alpha*ke2
	              - 2*(1.-dth*xi)*dth*alpha*ke3;
	double kea = delta_t*(dth*alpha*ke2-(1.-dth*xi)*ke3);
	double kex = delta_t*((dth*xi-1.)*ke1+dth*alpha*ke3);
	double cv = (xi-1./dth)*const_virial1 + alpha*const_virial2;
	double vf = 9.*volume*volume*b_mass;
	double h1 = alpha - (*b_alpha)
                          - dth*((ke+cv)/vf - (xi+3.*alpha)*alpha);
	double h1a = -1+dth*((kea+const_virial2)/vf-xi-6.*alpha);
	double h1x = dth*((kex+const_virial1)/vf-alpha);
	double h2 = xi - (*t_xi) - dth*(ke+vf*alpha*alpha)/t_mass;
	double h2a = dth*(kea+2.*vf*alpha)/t_mass;
	double h2x = -1+dth*kex/t_mass;
	double da = 0.;
	double dx = 0.;
	double converged;
	if (thermostat && barostat) {
	  da = (h2x*h1-h1x*h2)/(h2x*h1a-h1x*h2a);
	  dx = (h2-h2a*da)/h2x;
	}
	else if (thermostat)
	  dx = h2/h2x;
	else
	  da = h1/h1a;
	xi += dx;
	alpha += da;
	converged = 1;
	if (thermostat)
	  converged = converged && fabs(dx) < tol*fabs(xi);
	if (barostat)
	  converged = converged && fabs(da) < tol*fabs(alpha);
	if (converged) {
#if DEBUG
	  printf("%d iterations, %lf/%lf\n", 1000-niter, h1, h2);
#endif
	  break;
	}
	if (niter-- == 0) {
#if DEBUG
	  printf("No convergence\n");
#endif
	  break;
	}
      }
      if (thermostat) (*t_xi) = xi;
      if (barostat) (*b_alpha) = alpha;
    }
    /* Final velocity calculation */
    for (j = 0; j < atoms; j++)
      if (!fix[j]) {
	if (n_const > 0) {
	  v[j][0] -= v1[j][0];
	  v[j][1] -= v1[j][1];
	  v[j][2] -= v1[j][2];
	  if (barostat) {
	    vh[j][0] -= v2[j][0];
	    vh[j][1] -= v2[j][1];
	    vh[j][2] -= v2[j][2];
	  }
	}
	v[j][0] = (1.-dth*(*t_xi))*v[j][0] - dth*(*b_alpha)*vh[j][0];
	v[j][1] = (1.-dth*(*t_xi))*v[j][1] - dth*(*b_alpha)*vh[j][1];
	v[j][2] = (1.-dth*(*t_xi))*v[j][2] - dth*(*b_alpha)*vh[j][2];
      }

    printf("integrateRKM END - step  %d\n", i+1);

    /* The End - next time step! */
    time += delta_t;
  }
  /** End of main integration loop **/

  /* Final thermodynamic property evaluation */
  k_energy = 0.;
  for (j = 0; j < atoms; j++)
    if (!fix[j])
      k_energy += m[j]*dot(v[j], v[j]);
  k_energy *= 0.5;
  n_energy = 0.5*(*t_xi)*(*t_xi)*t_mass*exp(2*(*t_lns))
                + t_energy*(*t_lns);
  a_energy = 4.5*b_mass*volume*volume*(*b_alpha)*(*b_alpha) + volume*b_press;
  temperature = 2.*k_energy*temperature_factor/df;
  pressure = (2.*k_energy+p_energy.virial
	      + ((*t_xi)-1./dth)*const_virial1
	      + (*b_alpha)*const_virial2)/(3.*volume);

  /* Final trajectory and log output */
  if (PyTrajectory_Output(output, i, data_descriptors,
			  &evaluator->tstate_save) == -1) {
    PyUniverseSpec_StateLock(universe_spec, -2);
    PyEval_RestoreThread(evaluator->tstate_save);
    goto error;
  }

  /* Cleanup */
  PyUniverseSpec_StateLock(universe_spec, -2);
  PyEval_RestoreThread(evaluator->tstate_save);
  PyTrajectory_OutputFinish(output, i, 0, 1, data_descriptors);
  free(scratch);
  free(pdata);
  free(data_descriptors);
  Py_DECREF(gradients);
  Py_INCREF(Py_None);
  return Py_None;

  /* Error return */
error:
  PyTrajectory_OutputFinish(output, i, 1, 1, data_descriptors);
error2:
  free(scratch);
  free(pdata);
  free(data_descriptors);
  Py_DECREF(gradients);
  return NULL;
}

#include "bTD_INTER.h"  
#include "bTD_INTER2.h"  
#include "bTD_INTER3.h"  
#include "bTD_SPLIT.h"
#include "bTD_DEBUG.h"
#include "bTD_CLIENT.h"








/* Trajectory functions */

static int
getMassesAndVelocities(PyTrajectoryVariable *dynamic_data,
		       PyArrayObject **masses, PyArrayObject **velocities)
{
  PyTrajectoryVariable *var = dynamic_data;
  int found = 0;
  while (var->name != NULL) {
    if (strcmp(var->name, "masses") == 0) {
      *masses = var->value.array;
      found++;
    }
    if (strcmp(var->name, "velocities") == 0) {
      *velocities = var->value.array;
      found++;
    }
    var++;
  }
  if (found != 2) {
    PyErr_SetString(PyExc_ValueError,
		    "trajectory function needs masses and velocities");
    return 0;
  }
  return 1;
}

static int
getDegreesOfFreedom(PyTrajectoryVariable *dynamic_data)
{
  int natoms = -1, ndf = -1;
  PyTrajectoryVariable *var = dynamic_data;
  while (var->name != NULL) {
    if (strcmp(var->name, "degrees_of_freedom") == 0)
      ndf = *var->value.ip;
    if (strcmp(var->name, "configuration") == 0)
      natoms = var->value.array->dimensions[0];
    var++;
  }
  if (ndf >= 0)
    return ndf;
  else
    return 3*natoms;
}

static double *
getScalar(PyTrajectoryVariable *dynamic_data, char *name)
{
  PyTrajectoryVariable *var = dynamic_data;
  double *variable = NULL;
  while (var->name != NULL) {
    if (strcmp(var->name, name) == 0) {
      variable = var->value.dp;
    }
    var++;
  }
  return variable;
}

static PyArrayObject *
getConfiguration(PyTrajectoryVariable *dynamic_data)
{
  PyTrajectoryVariable *var = dynamic_data;
  PyArrayObject *configuration = NULL;
  while (var->name != NULL) {
    if (strcmp(var->name, "configuration") == 0) {
      configuration = var->value.array;
    }
    var++;
  }
  return configuration;
}

/* Temperature scaling option */

static int
scaleVelocities(PyTrajectoryVariable *dynamic_data, PyObject *parameters,
		int step, void **scratch, PyObject *universe)
{
  PyArrayObject *p_array = (PyArrayObject *)parameters;
  double temperature = ((double *)p_array->data)[0];
  double window = ((double *)p_array->data)[1];
  typedef struct {PyArrayObject *masses;
                  PyArrayObject *velocities;
                  double *thermostat;
                  int df;} workspace;
  workspace *ws = (workspace *)*scratch;

  if (step == -1) {
    ws = (workspace *)malloc(sizeof(workspace));
    *scratch = (void *)ws;
    if (ws == NULL) {
      PyErr_NoMemory();
      return 0;
    }
    if (!getMassesAndVelocities(dynamic_data, &ws->masses, &ws->velocities))
      return 0;
    ws->thermostat = getScalar(dynamic_data, "thermostat_coordinate");
    ws->df = getDegreesOfFreedom(dynamic_data);
  }

  else if (step == -2) {
    free(ws);
  }

  else {  /* Scale velocities */
    vector3 *v = (vector3 *)ws->velocities->data;
    int atoms = ws->velocities->dimensions[0];
    double *m = (double *)ws->masses->data;
    double k_energy = 0.;
    double t;
    int j;
    for (j = 0; j < atoms; j++)
      k_energy += m[j]*(v[j][0]*v[j][0]+v[j][1]*v[j][1]+v[j][2]*v[j][2]);
    t = k_energy*temperature_factor/ws->df;
    if (t > 0. && fabs(t-temperature) > window) {
      double f = sqrt(temperature/t);
      for (j = 0; j < atoms; j++) {
	v[j][0] *= f;
	v[j][1] *= f;
	v[j][2] *= f;
      }
    }
    k_energy = 0.;
    for (j = 0; j < atoms; j++)
      k_energy += m[j]*(v[j][0]*v[j][0]+v[j][1]*v[j][1]+v[j][2]*v[j][2]);
    t = k_energy*temperature_factor/ws->df;
    if (ws->thermostat != NULL) {
      ws->thermostat[0] = 0.;
      ws->thermostat[1] = 0.;
    }
  }
  return 1;
}

static int
heat(PyTrajectoryVariable *dynamic_data, PyObject *parameters,
     int step, void **scratch, PyObject *universe)
{
  PyArrayObject *p_array = (PyArrayObject *)parameters;
  double temp1 = ((double *)p_array->data)[0];
  double temp2 = ((double *)p_array->data)[1];
  double gradient = ((double *)p_array->data)[2];
  typedef struct {PyArrayObject *masses;
                  PyArrayObject *velocities;
                  double *thermostat;
                  double *time;
                  int df;} workspace;
  workspace *ws = (workspace *)*scratch;
  if (step == -1) {  /* Initialization */
    ws = (workspace *)malloc(sizeof(workspace));
    *scratch = (void *)ws;
    if (ws == NULL) {
      PyErr_NoMemory();
      return 0;
    }
    if (!getMassesAndVelocities(dynamic_data, &ws->masses, &ws->velocities))
      return 0;
    ws->thermostat = getScalar(dynamic_data, "thermostat_coordinate");
    if (ws->thermostat != NULL) {
      PyErr_SetString(PyExc_ValueError,
		      "heating not allowed with thermostat");
      return 0;
    }
    ws->df = getDegreesOfFreedom(dynamic_data);
    ws->time = getScalar(dynamic_data, "time");
  }

  else if (step == -2) {  /* Clean up */
    free(ws);
  }

  else {  /* Scale velocities */
    vector3 *v = (vector3 *)ws->velocities->data;
    int atoms = ws->velocities->dimensions[0];
    double *m = (double *)ws->masses->data;
    double temperature = temp1 + (*ws->time) * gradient;
    double k_energy = 0.;
    double t, f;
    int j;
    if ((gradient > 0. && temperature > temp2)
	|| (gradient < 0. && temperature < temp2))
      temperature = temp2;
    for (j = 0; j < atoms; j++)
      k_energy += m[j]*(v[j][0]*v[j][0]+v[j][1]*v[j][1]+v[j][2]*v[j][2]);
    t = k_energy*temperature_factor/ws->df;
    if (t > 0.){
      f = sqrt(temperature/t);
      for (j = 0; j < atoms; j++) {
	v[j][0] *= f;
	v[j][1] *= f;
	v[j][2] *= f;
      }
    }
  }
  return 1;
}

/* Resetting the barostat */

static int
resetBarostat(PyTrajectoryVariable *dynamic_data, PyObject *parameters,
	      int step, void **scratch, PyObject *universe)
{
  typedef struct {double *barostat;} workspace;
  workspace *ws = (workspace *)*scratch;

  if (step == -1) {
    ws = (workspace *)malloc(sizeof(workspace));
    *scratch = (void *)ws;
    if (ws == NULL) {
      PyErr_NoMemory();
      return 0;
    }
    ws->barostat = getScalar(dynamic_data, "barostat_coordinate");
    if (ws->barostat == NULL) {
      PyErr_SetString(PyExc_ValueError,
		      "no barostat to reset");
      return 0;
    }
  }

  else if (step == -2) {
    free(ws);
  }

  else {  /* Reset barostat */
    *ws->barostat = 0.;
  }
  return 1;
}

/* Elimination of global translation */

static int
removeTranslation(PyTrajectoryVariable *dynamic_data, PyObject *parameters,
		  int step, void **scratch, PyObject *universe)
{
  typedef struct {PyArrayObject *masses;
                  PyArrayObject *velocities;} workspace;
  workspace *ws = (workspace *)*scratch;

  if (step == -1) {
    ws = (workspace *)malloc(sizeof(workspace));
    *scratch = (void *)ws;
    if (ws == NULL) {
      PyErr_NoMemory();
      return 0;
    }
    if (!getMassesAndVelocities(dynamic_data, &ws->masses, &ws->velocities))
      return 0;
  }
  else if (step == -2) {
    free(ws);
  }
  else {  /* Remove translation */
    vector3 *v = (vector3 *)ws->velocities->data;
    int atoms = ws->velocities->dimensions[0];
    double *m = (double *)ws->masses->data;
    double total_mass, momentum;
    int i, j;
    total_mass = 0.;
    for (j = 0; j < atoms; j++)
      total_mass += m[j];
    for (i = 0; i < 3; i++) {
      momentum = 0.;
      for (j = 0; j < atoms; j++)
	momentum += m[j]*v[j][i];
      momentum /= total_mass;
      for (j = 0; j < atoms; j++)
	v[j][i] -= momentum;
    }
  }
  return 1;
}

/* Eliminitation of global rotation */

static void
solve_3x3(tensor3 A, vector3 B, vector3 X)
{
  double a = A[0][0];
  double b = A[1][1];
  double c = A[2][2];
  double d = A[0][1];
  double e = A[0][2];
  double f = A[1][2];
  double o = B[0];
  double p = B[1];
  double q = B[2];

  double af_de = a*f-d*e;
  double aq_eo = a*q-e*o;
  double ab_dd = a*b-d*d;
  double ac_ee = a*c-e*e;

  double z = (af_de*(a*p-d*o)-ab_dd*aq_eo) / (af_de*af_de-ab_dd*ac_ee);
  double y = (aq_eo - z*ac_ee)/af_de;
  double x = (o - d*y - e*z)/a;

  X[0] = x;
  X[1] = y;
  X[2] = z;
}

static int
removeRotation(PyTrajectoryVariable *dynamic_data, PyObject *parameters,
	       int step, void **scratch, PyObject *universe)
{
  typedef struct {PyArrayObject *configuration;
                  PyArrayObject *masses;
                  PyArrayObject *velocities;} workspace;
  workspace *ws = (workspace *)*scratch;

  if (step == -1) {  /* Initialization */
    ws = (workspace *)malloc(sizeof(workspace));
    *scratch = (void *)ws;
    if (ws == NULL) {
      PyErr_NoMemory();
      return 0;
    }
    if (!getMassesAndVelocities(dynamic_data, &ws->masses, &ws->velocities))
      return 0;
    ws->configuration = getConfiguration(dynamic_data);
    if (ws->configuration == NULL) {
      PyErr_SetString(PyExc_ValueError,
		      "rotation remover needs configuration");
      return 0;
    }
  }
  else if (step == -2) {  /* Clean up */
    free(ws);
  }
  else {  /* Remove rotation */
    vector3 *v = (vector3 *)ws->velocities->data;
    vector3 *x = (vector3 *)ws->configuration->data;
    double *m = (double *)ws->masses->data;
    int atoms = ws->masses->dimensions[0];

    double total_mass, trace;
    vector3 cm, l, o;
    tensor3 inertia;
    int i, j, k;

    total_mass = 0.;
    cm[0] = cm[1] = cm[2] = 0.;
    for (i = 0; i < atoms; i++) {
      total_mass += m[i];
      cm[0] += m[i]*x[i][0];
      cm[1] += m[i]*x[i][1];
      cm[2] += m[i]*x[i][2];
    }
    cm[0] /= total_mass; cm[1] /= total_mass; cm[2] /= total_mass;

    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
	inertia[j][k] = 0.;
    l[0] = l[1] = l[2] = 0.;
    for (i = 0; i < atoms; i++) {
      vector3 temp, l1;
      tensor3 i1;
      temp[0] = x[i][0]-cm[0];
      temp[1] = x[i][1]-cm[1];
      temp[2] = x[i][2]-cm[2];
      cross(l1, temp, v[i]);
      vector_scale(l1, m[i]);
      l[0] += l1[0]; l[1] += l1[1]; l[2] += l1[2];
      tensor_product(i1, temp, temp, m[i]);
      for (j = 0; j < 3; j++)
	for (k = 0; k < 3; k++)
	  inertia[j][k] -= i1[j][k];
    }
    trace = inertia[0][0] + inertia[1][1] + inertia[2][2];
    inertia[0][0] -= trace;
    inertia[1][1] -= trace;
    inertia[2][2] -= trace;

    solve_3x3(inertia, l, o);

    for (i = 0; i < atoms; i++) {
      vector3 temp, v1;
      temp[0] = x[i][0]-cm[0];
      temp[1] = x[i][1]-cm[1];
      temp[2] = x[i][2]-cm[2];
      cross(v1, o, temp);
      v[i][0] -= v1[0];
      v[i][1] -= v1[1];
      v[i][2] -= v1[2];
    }

  }
  return 1;
}

/*
 * Enforce distance constraints by calling SHAKE
 */
static PyObject *
enforceConstraints(PyObject *dummy, PyObject *args)
{
  PyUniverseSpecObject *universe_spec;
  PyArrayObject *configuration;
  PyArrayObject *masses;
  PyArrayObject *constraints, *constraint_distances_squared, *c_blocks;

  vector3 *x;
  double *m;
  int n_const, n_const_blocks;
  long *const_pairs, *const_blocks;
  double *const_dist;
  vector3 *const_vect = NULL;
  int i;

  if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!",
			&PyUniverseSpec_Type, &universe_spec,
			&PyArray_Type, &configuration,
			&PyArray_Type, &masses,
			&PyArray_Type, &constraints,
			&PyArray_Type, &constraint_distances_squared,
			&PyArray_Type, &c_blocks))
    return NULL;

  n_const = constraints->dimensions[0];
  n_const_blocks = c_blocks->dimensions[0]-1;
  x = (vector3 *)configuration->data;
  m = (double *)masses->data;
  const_pairs = (long *)constraints->data;
  const_dist = (double *)constraint_distances_squared->data;
  const_blocks = (long *)c_blocks->data;

  const_vect = (vector3 *)malloc(n_const*sizeof(vector3));
  if (const_vect == NULL) {
    PyErr_NoMemory();
    return NULL;
  }

  for (i = 0; i < n_const; i++)
    universe_spec->distance_function(const_vect[i],
				     x[const_pairs[2*i]],
				     x[const_pairs[2*i+1]],
				     universe_spec->geometry_data);
  for (i = 0; i < n_const_blocks; i++)
    shake(const_pairs, const_blocks[i], const_blocks[i+1],
	  x, m, const_vect, const_dist, universe_spec->distance_function,
	  universe_spec->geometry_data);

  free(const_vect);
  Py_INCREF(Py_None);
  return Py_None;
}

/*
 * Project velocities onto constraint surface
 */
static PyObject *
projectVelocities(PyObject *dummy, PyObject *args)
{
  PyUniverseSpecObject *universe_spec;
  PyArrayObject *configuration;
  PyArrayObject *velocities;
  PyArrayObject *masses;
  PyArrayObject *constraints, *constraint_distances_squared, *c_blocks;

  int atoms;
  vector3 *x, *v;
  double *m;
  int n_const;
  long *const_pairs;
  double *const_dist;
  vector3 *const_vect = NULL;
  vector3 *vc = NULL;
  projection_data *pdata = NULL;
  int i;

  if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!",
			&PyUniverseSpec_Type, &universe_spec,
			&PyArray_Type, &configuration,
			&PyArray_Type, &velocities,
			&PyArray_Type, &masses,
			&PyArray_Type, &constraints,
			&PyArray_Type, &constraint_distances_squared,
			&PyArray_Type, &c_blocks))
    return NULL;

  atoms = configuration->dimensions[0];
  n_const = constraints->dimensions[0];
  x = (vector3 *)configuration->data;
  v = (vector3 *)velocities->data;
  m = (double *)masses->data;
  const_pairs = (long *)constraints->data;
  const_dist = (double *)constraint_distances_squared->data;

  pdata = (projection_data *)malloc(n_const*sizeof(projection_data));
  const_vect = (vector3 *)malloc(n_const*sizeof(vector3));
  vc = (vector3 *)malloc(atoms*sizeof(vector3));
  if (pdata == NULL || const_vect == NULL || vc == NULL) {
    free(pdata);
    free(const_vect);
    free(vc);
    PyErr_NoMemory();
    return NULL;
  }

  for (i = 0; i < n_const; i++) {
    universe_spec->distance_function(const_vect[i],
				     x[const_pairs[2*i]],
				     x[const_pairs[2*i+1]],
				     universe_spec->geometry_data);
    pdata[i].multiplier[MU1] = 0.;
    pdata[i].diag = const_dist[i] *
	(1./m[const_pairs[2*i]] + 1./m[const_pairs[2*i+1]]);
  }
  project(n_const, const_pairs, const_dist, const_vect,
	  pdata, MU1, m, v, vc, atoms);
  for (i = 0; i < atoms; i++) {
    v[i][0] -= vc[i][0];
    v[i][1] -= vc[i][1];
    v[i][2] -= vc[i][2];
  }

  free(pdata);
  free(const_vect);
  free(vc);
  Py_INCREF(Py_None);
  return Py_None;
}


/*
 * List of functions defined in the module
 */

static PyMethodDef dynamics_methods[] = {
  {"integrateRKM", integrateRKM, 1},
  {"integrateRKM_EVOL", integrateRKM_EVOL, 1},
  {"integrateRKM_INTER", integrateRKM_INTER, 1},
  {"integrateRKM_INTER2", integrateRKM_INTER2, 1},
  {"integrateRKM_INTER3", integrateRKM_INTER3, 1},
  {"integrateRKM_SPLIT", integrateRKM_SPLIT, 1},
  {"integrateRKM_DEBUG", integrateRKM_DEBUG, 1},
  {"integrateRKM_client", integrateRKM_client, 1},
  {"enforceConstraints", enforceConstraints, 1},
  {"projectVelocities", projectVelocities, 1},
  {NULL, NULL}		/* sentinel */
};

/* Initialization function for the module */

DL_EXPORT(void)
initbTD_dynamics(void)
{
  PyObject *m, *dict;
  PyObject *units;

  /* Create the module and add the functions */
  m = Py_InitModule("bTD_dynamics", dynamics_methods);
  dict = PyModule_GetDict(m);

  /* Import the array module */
#ifdef import_array
  import_array();
#endif

  /* Import MMTK modules */
  import_MMTK_universe();
  import_MMTK_forcefield();
  import_MMTK_trajectory();

  /* Get the temperature conversion factor from Units */
  units = PyImport_ImportModule("MMTK.Units");
  if (units != NULL) {
    PyObject *module_dict = PyModule_GetDict(units);
    PyObject *factor = PyDict_GetItemString(module_dict, "k_B");
    kB = PyFloat_AsDouble(factor);
    temperature_factor = 1./kB;
  }

  /* Add addresses of C functions to module dictionary */
  PyDict_SetItemString(dict, "scaleVelocities",
		       PyCObject_FromVoidPtr((void *)scaleVelocities, NULL));
  PyDict_SetItemString(dict, "heat",
		       PyCObject_FromVoidPtr((void *)heat, NULL));
  PyDict_SetItemString(dict, "resetBarostat",
		       PyCObject_FromVoidPtr((void *)resetBarostat, NULL));
  PyDict_SetItemString(dict, "removeTranslation",
		       PyCObject_FromVoidPtr((void *)removeTranslation, NULL));
  PyDict_SetItemString(dict, "removeRotation",
		       PyCObject_FromVoidPtr((void *)removeRotation, NULL));

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module bTD_dynamics");
}
