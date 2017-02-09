#include "MMTK/universe.h"
#include "MMTK/forcefield.h"
#include "MMTK/forcefield_private.h"

/* This function does the actual energy (and gradient) calculation.
   Everything else is just bookkeeping. */
static void
ef_evaluator(PyFFEnergyTermObject *self,
	     PyFFEvaluatorObject *eval,
	     energy_spec *input,
	     energy_data *energy)
     /* The four parameters are pointers to structures that are
	defined in MMTK/forcefield.h.
	PyFFEnergyTermObject: All data relevant to this particular
                              energy term.
        PyFFEvaluatorObject:  Data referring to the global energy
                              evaluation process, e.g. parallelization
                              options. Not used here.
        energy_spec:          Input parameters for this routine, i.e.
                              atom positions and parallelization parameters.
        energy_data:          Storage for the results (energy terms,
                              gradients, second derivatives).
     */
{
  // Input variables
  vector3 *coordinates = (vector3 *)input->coordinates->data;
  int natoms = input->coordinates->dimensions[0];
  vector3 *g;
  
  double strength = self->param[0];
  double k = self->param[1];
  int nyz = (int) self->param[2];
  
  vector3 hCorner;
  hCorner[0] = self->param[3];
  hCorner[1] = self->param[4];
  hCorner[2] = self->param[5];
  
  PyArrayObject *spacing_array = (PyArrayObject *)self->data[3];
  double* spacing = (double *)spacing_array->data;
  PyArrayObject *counts_array = (PyArrayObject *)self->data[4];
  long* counts = (long *)counts_array->data;
  PyArrayObject *vals_array = (PyArrayObject *)self->data[5];
  double* vals = (double *)vals_array->data;
  PyArrayObject *scaling_factor_array = (PyArrayObject *)self->data[6];
  double* scaling_factor = (double *)scaling_factor_array->data;
  
  /* energy_terms is an array because each routine could compute
     several terms that should logically be kept apart. For example,
     a single routine calculates Lennard-Jones and electrostatic interactions
     in a single iteration over the nonbonded list. The separation of
     terms is only done for the benefit of user code (universe.energyTerms())
     returns the value of each term separately), the total energy is
     always the sum of all terms. Here we have only one energy term,
     which is initialized to zero.
     Note that the virial is also stored in the array energy_terms,
     at the index self->virial_index. However, there is only one virial
     term for the whole system, it is not added up term by term. Therefore
     we don't set it to zero, we just add to it in the loop. */
  // energy->energy_terms[self->index] = 0.;

  // Variables for output
  double gridEnergy = 0.;

  /* Add the gradient contribution to the global gradient array.
     It would be a serious error to use '=' instead of '+=' here,
     in that case all previously calculated forces would be erased.
     If energy_gradients is NULL, then the calling routine does not
     want gradients, and didn't provide storage for them.
     Second derivatives are not calculated because they are zero. */
  if (energy->gradients != NULL)
    g = (vector3 *)((PyArrayObject*)energy->gradients)->data;
  
  // Variables for processing
  int i, ix, iy, iz, ind;
  double vmmm, vmmp, vmpm, vmpp, vpmm, vpmp, vppm, vppp;
  double vmm, vmp, vpm, vpp, vm, vp;
  double fx, fy, fz, ax, ay, az;
  double dvdx, dvdy, dvdz;
  double dev;
  double interpolated, prefactor;
  
  for (ind = 0; ind < natoms; ind++) {
    if ((coordinates[ind][0]>0.) && (coordinates[ind][1]>0.)
        && (coordinates[ind][2]>0.) &&
        (coordinates[ind][0]<hCorner[0]) && (coordinates[ind][1]<hCorner[1])
        && (coordinates[ind][2]<hCorner[2]))
    {
      // Index within the grid
      ix = (int) (coordinates[ind][0]/spacing[0]);
      iy = (int) (coordinates[ind][1]/spacing[1]);
      iz = (int) (coordinates[ind][2]/spacing[2]);
      
      i = ix*nyz + iy*counts[2] + iz;
      
      // Corners of the box surrounding the point
      vmmm = vals[i];
      vmmp = vals[i+1];
      vmpm = vals[i+counts[2]];
      vmpp = vals[i+counts[2]+1];

      vpmm = vals[i+nyz];
      vpmp = vals[i+nyz+1];
      vppm = vals[i+nyz+counts[2]];
      vppp = vals[i+nyz+counts[2]+1];
      
      // Fraction within the box
      fx = (coordinates[ind][0] - (ix*spacing[0]))/spacing[0];
      fy = (coordinates[ind][1] - (iy*spacing[1]))/spacing[1];
      fz = (coordinates[ind][2] - (iz*spacing[2]))/spacing[2];
      
      // Fraction ahead
      ax = 1 - fx;
      ay = 1 - fy;
      az = 1 - fz;

      // Trilinear interpolation for energy
      vmm = az*vmmm + fz*vmmp;
      vmp = az*vmpm + fz*vmpp;
      vpm = az*vpmm + fz*vpmp;
      vpp = az*vppm + fz*vppp;
      
      vm = ay*vmm + fy*vmp;
      vp = ay*vpm + fy*vpp;
      
      interpolated = (ax*vm + fx*vp);
      gridEnergy += scaling_factor[ind]*interpolated*interpolated*interpolated*interpolated;
      if (energy->gradients != NULL) {
        // x coordinate
        dvdx = -vm + vp;
        // y coordinate
        dvdy = (-vmm + vmp)*ax + (-vpm + vpp)*fx;
        // z coordinate
        dvdz = ((-vmmm + vmmp)*ay + (-vmpm + vmpp)*fy)*ax +
               ((-vpmm + vpmp)*ay + (-vppm + vppp)*fy)*fx;
      
        prefactor = strength*scaling_factor[ind]*4.*interpolated*interpolated*interpolated;
        g[ind][0] += prefactor*dvdx/spacing[0];
        g[ind][1] += prefactor*dvdy/spacing[1];
        g[ind][2] += prefactor*dvdz/spacing[2];
      }
    }
    else {
      for (i = 0; i<3; i++) {
        if (coordinates[ind][i]<i) {
          gridEnergy += k*coordinates[ind][i]*coordinates[ind][i]/2.;
          if (energy->gradients != NULL)
            g[ind][i] += k*coordinates[ind][i];
        }
        else {
          if (coordinates[ind][i]>hCorner[i]) {
            dev = (coordinates[ind][i]-hCorner[i]);
            gridEnergy += k*dev*dev/2.;
            if (energy->gradients != NULL)
              g[ind][i] += k*dev;
          }
        }
      }
    }
    energy->energy_terms[self->index] = gridEnergy*strength;
  }
}

/* A utility function that allocates memory for a copy of a string */
static char *
allocstring(char *string)
{
  char *memory = (char *)malloc(strlen(string)+1);
  if (memory != NULL)
    strcpy(memory, string);
  return memory;
}

/* The next function is meant to be called from Python. It creates the
   energy term object at the C level and stores all the parameters in
   there in a form that is convient to access for the C routine above.
   This is the routine that is imported into and called by the Python
   module. */
static PyObject *
TrilinearOneFourthGridTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self;
  PyArrayObject *spacing;
  PyArrayObject *counts;
  PyArrayObject *vals;
  double strength;
  PyArrayObject *scaling_factor;
  char *name;
  
  /* Create a new energy term object and return if the creation fails. */
  self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;
  /* Convert the parameters to C data types. */
  if (!PyArg_ParseTuple(args, "O!O!O!O!dO!s",
			&PyUniverseSpec_Type, &self->universe_spec,
      &PyArray_Type, &spacing,
      &PyArray_Type, &counts,
      &PyArray_Type, &vals,
      &strength,
      &PyArray_Type, &scaling_factor,
			&name))
    return NULL;
  
  /* We keep a reference to the universe_spec in the newly created
     energy term object, so we have to increase the reference count. */
  Py_INCREF(self->universe_spec);
  /* A pointer to the evaluation routine. */
  self->eval_func = ef_evaluator;
  /* The name of the energy term object. */
  self->evaluator_name = "trilinear_one_fourth_grid";
  /* The names of the individual energy terms - just one here. */
  self->term_names[0] = allocstring(name);
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->nterms = 1;

  long* counts_v = (long* )counts->data;
  double* spacing_v = (double* )spacing->data;

//  int ind;
//
//  double* vals_v = (double *)vals->data;
//
//  printf("Counts:\n");
//  for (ind = 0; ind < 3; ind ++)
//    printf("count %d = %d\n", ind, counts_v[ind]);
//  fflush(stdout);
//
//  printf("Spacing:\n");
//  for (ind = 0; ind < 3; ind ++)
//    printf("spacing %d = %f\n", ind, spacing_v[ind]);
//  fflush(stdout);
//
//  printf("Values:\n");
//  for (ind = 0; ind < 10; ind ++)
//    printf("vals %d = %f\n", ind, vals_v[ind]);
//  fflush(stdout);
  
  /* self->param is a storage area for parameters. Note that there
     are only 40 slots (double) there. */
  self->param[0] = strength;
  self->param[1] = 10000.; // k, the spring constant in kJ/mol nm**2
  self->param[2] = counts_v[1]*counts_v[2]; // nyz
  self->param[3] = spacing_v[0]*(counts_v[0]-1); // hCorner in x
  self->param[4] = spacing_v[1]*(counts_v[1]-1); // hCorner in y
  self->param[5] = spacing_v[2]*(counts_v[2]-1); // hCorner in z
  
//  printf("Params:\n");
//  for (ind = 0; ind < 6; ind ++)
//    printf("param %d = %f\n", ind, self->param[ind]);
//  fflush(stdout);
  
  /* self->data is the other storage area for parameters. There are
     40 Python object slots there */
  self->data[3] = (PyObject *)spacing;
  Py_INCREF(spacing);
  self->data[4] = (PyObject *)counts;
  Py_INCREF(counts);
  self->data[5] = (PyObject *)vals;
  Py_INCREF(vals);
  self->data[6] = (PyObject *)scaling_factor;
  Py_INCREF(scaling_factor);
  
  /* Return the energy term object. */
  return (PyObject *)self;
}

/* This is a list of all Python-callable functions defined in this
   module. Each list entry consists of the name of the function object
   in the module, the C routine that implements it, and a "1" signalling
   new-style parameter passing conventions (only veterans care about the
   alternatives). The list is terminated by a NULL entry. */
static PyMethodDef functions[] = {
  {"TrilinearOneFourthGridTerm", TrilinearOneFourthGridTerm, 1},
  {NULL, NULL}		/* sentinel */
};


/* The initialization function for the module. This is the only function
   that must be publicly visible, everything else should be declared
   static to prevent name clashes with other modules. The name of this
   function must be "init" followed by the module name. */
DL_EXPORT(void)
initMMTK_trilinear_one_fourth_grid(void)
{
  PyObject *m;

  /* Create the module and add the functions. */
  m = Py_InitModule("MMTK_trilinear_one_fourth_grid", functions);

  /* Import the array module. */
#ifdef import_array
  import_array();
#endif

  /* Import MMTK modules. */
  import_MMTK_universe();
  import_MMTK_forcefield();

  /* Check for errors. */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_trilinear_one_fourth_grid");
}
