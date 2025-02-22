/* C routines for OBC.py */

#include "MMTK/universe.h"
#include "MMTK/forcefield.h"
#include "MMTK/forcefield_private.h"
#include "ObcWrapper.h"

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
  vector3* coordinates = (vector3 *)input->coordinates->data;
  vector3* g;

  int natoms = input->coordinates->dimensions[0];
  double Igrid[natoms];
  double fractionalDesolvationToIgrid = self->param[6];

  /* Interpolate the fractional desolvation grid */
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

  // Variables for processing
  int i, ix, iy, iz, ind;
  double vmmm, vmmp, vmpm, vmpp, vpmm, vpmp, vppm, vppp;
  double vmm, vmp, vpm, vpp, vm, vp;
  double fx, fy, fz, ax, ay, az;

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
      
      // (ax*vm + fx*vp) is the interpolation of fractional desolvation
      // fractionalDesolvationToIgrid converts it to Igrid
      Igrid[ind] = fractionalDesolvationToIgrid*(ax*vm + fx*vp);
    }
  }
  
  struct ObcParameters* obcParameters = (struct ObcParameters*)self->data[6];
  struct ReferenceObc* obc = (struct ReferenceObc*)self->data[7];
  
  if (energy->gradients != NULL) {
    g = (vector3 *)((PyArrayObject*)energy->gradients)->data;
    energy->energy_terms[self->index] =
      computeBornEnergyForces(obc, obcParameters, Igrid, coordinates, g);
  }
  else {
    energy->energy_terms[self->index] =
      computeBornEnergy(obc, obcParameters, (const double*)Igrid, coordinates);
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
   module, OBC.py. */
static PyObject *
OBCDesolvTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self;
  int numParticles;
  PyArrayObject *charges;
  PyArrayObject *atomicRadii;
  PyArrayObject *scaleFactors;
  PyArrayObject *spacing;
  PyArrayObject *counts;
  PyArrayObject *vals;
  double strength;
  // Integration range for the Coulomb integral
  double r_min;
  double r_max;

  /* Create a new energy term object and return if the creation fails. */
  self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;
  /* Convert the parameters to C data types. */
  if (!PyArg_ParseTuple(args, "O!idO!O!O!O!O!O!dd",
			&PyUniverseSpec_Type, &self->universe_spec,
      &numParticles, &strength,
			&PyArray_Type, &charges,
      &PyArray_Type, &atomicRadii,
      &PyArray_Type, &scaleFactors,
      &PyArray_Type, &spacing,
      &PyArray_Type, &counts,
      &PyArray_Type, &vals,
      &r_min, &r_max))
    return NULL;
  /* We keep a reference to the universe_spec in the newly created
     energy term object, so we have to increase the reference count. */
  Py_INCREF(self->universe_spec);
  /* A pointer to the evaluation routine. */
  self->eval_func = ef_evaluator;
  /* The name of the energy term object. */
  self->evaluator_name = "OBC_desolv";
  /* The names of the individual energy terms - just one here. */
  self->term_names[0] = allocstring("OBC_desolv");
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->nterms = 1;
  
  struct ObcParameters* obcParameters = newObcParameters(
    numParticles, strength, (double *)charges->data,
    (double *)atomicRadii->data, (double *)scaleFactors->data);
  struct ReferenceObc* obc = newReferenceObc(obcParameters);

  long* counts_v = (long* )counts->data;
  double* spacing_v = (double* )spacing->data;

  /* self->param is a storage area for parameters. Note that there
     are only 40 slots (double) there. */
  self->param[0] = strength;
  self->param[1] = (double) numParticles;
  self->param[2] = counts_v[1]*counts_v[2]; // nyz
  self->param[3] = spacing_v[0]*(counts_v[0]-1); // hCorner in x
  self->param[4] = spacing_v[1]*(counts_v[1]-1); // hCorner in y
  self->param[5] = spacing_v[2]*(counts_v[2]-1); // hCorner in z
  // Multiplicative factor for fractional desolvation to get Igrid
  // TODO: Understand where the factor of 4*pi comes from
  self->param[6] = 4*3.14159265359*(1/r_min - 1/r_max);
  
  /* self->data is the other storage area for parameters. There are
     40 Python object slots there */
  self->data[0] = (PyObject *)charges;
  Py_INCREF(charges);
  self->data[1] = (PyObject *)atomicRadii;
  Py_INCREF(atomicRadii);
  self->data[2] = (PyObject *)scaleFactors;
  Py_INCREF(scaleFactors);
  self->data[3] = (PyObject *)spacing;
  Py_INCREF(spacing);
  self->data[4] = (PyObject *)counts;
  Py_INCREF(counts);
  self->data[5] = (PyObject *)vals;
  Py_INCREF(vals);
  self->data[6] = (PyObject *)obcParameters;
  Py_INCREF(obcParameters); // Seems to increment the number of particles
  setNumberOfAtoms(obcParameters, numParticles);
  self->data[7] = (PyObject *)obc;
  Py_INCREF(obc);
  
  /* Return the energy term object. */
  return (PyObject *)self;
}

/* This is a list of all Python-callable functions defined in this
   module. Each list entry consists of the name of the function object
   in the module, the C routine that implements it, and a "1" signalling
   new-style parameter passing conventions (only veterans care about the
   alternatives). The list is terminated by a NULL entry. */
static PyMethodDef functions[] = {
  {"OBCDesolvTerm", OBCDesolvTerm, 1},
  {NULL, NULL}		/* sentinel */
};


/* The initialization function for the module. This is the only function
   that must be publicly visible, everything else should be declared
   static to prevent name clashes with other modules. The name of this
   function must be "init" followed by the module name. */
DL_EXPORT(void)
initMMTK_OBC_desolv(void)
{
  PyObject *m;

  /* Create the module and add the functions. */
  m = Py_InitModule("MMTK_OBC_desolv", functions);

  /* Import the array module. */
#ifdef import_array
  import_array();
#endif

  /* Import MMTK modules. */
  import_MMTK_universe();
  import_MMTK_forcefield();

  /* Check for errors. */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_OBC_desolv");
}
