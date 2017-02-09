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
  
  struct ObcParameters* obcParameters = (struct ObcParameters*)self->data[6];
  struct ReferenceObc* obc = (struct ReferenceObc*)self->data[7];
  
  if (energy->gradients != NULL) {
    g = (vector3 *)((PyArrayObject*)energy->gradients)->data;
    energy->energy_terms[self->index] =
      computeBornEnergyForces(obc, obcParameters, coordinates, g);
  }
  else {
    energy->energy_terms[self->index] =
      computeBornEnergy(obc, obcParameters, coordinates);
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
OBCTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self;
  int numParticles;
  PyArrayObject *charges;
  PyArrayObject *atomicRadii;
  PyArrayObject *scaleFactors;
  double strength;

  /* Create a new energy term object and return if the creation fails. */
  self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;
  /* Convert the parameters to C data types. */
  if (!PyArg_ParseTuple(args, "O!idO!O!O!",
			&PyUniverseSpec_Type, &self->universe_spec,
      &numParticles, &strength,
			&PyArray_Type, &charges,
      &PyArray_Type, &atomicRadii,
      &PyArray_Type, &scaleFactors))
    return NULL;
  /* We keep a reference to the universe_spec in the newly created
     energy term object, so we have to increase the reference count. */
  Py_INCREF(self->universe_spec);
  /* A pointer to the evaluation routine. */
  self->eval_func = ef_evaluator;
  /* The name of the energy term object. */
  self->evaluator_name = "OBC";
  /* The names of the individual energy terms - just one here. */
  self->term_names[0] = allocstring("OBC");
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->nterms = 1;
  
  struct ObcParameters* obcParameters = newObcParameters(
    numParticles, strength, (double *)charges->data,
    (double *)atomicRadii->data, (double *)scaleFactors->data);
  struct ReferenceObc* obc = newReferenceObc(obcParameters);

  /* self->param is a storage area for parameters. Note that there
     are only 40 slots (double) there. */
  self->param[0] = strength;
  self->param[1] = (double) numParticles;
  
  /* self->data is the other storage area for parameters. There are
     40 Python object slots there */
  self->data[0] = (PyObject *)charges;
  Py_INCREF(charges);
  self->data[1] = (PyObject *)atomicRadii;
  Py_INCREF(atomicRadii);
  self->data[2] = (PyObject *)scaleFactors;
  Py_INCREF(scaleFactors);
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
  {"OBCTerm", OBCTerm, 1},
  {NULL, NULL}		/* sentinel */
};


/* The initialization function for the module. This is the only function
   that must be publicly visible, everything else should be declared
   static to prevent name clashes with other modules. The name of this
   function must be "init" followed by the module name. */
DL_EXPORT(void)
initMMTK_OBC(void)
{
  PyObject *m;

  /* Create the module and add the functions. */
  m = Py_InitModule("MMTK_OBC", functions);

  /* Import the array module. */
#ifdef import_array
  import_array();
#endif

  /* Import MMTK modules. */
  import_MMTK_universe();
  import_MMTK_forcefield();

  /* Check for errors. */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_OBC");
}
