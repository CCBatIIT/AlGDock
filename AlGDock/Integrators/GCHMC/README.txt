GCHMC integrator requires:
  - a .mol2 structure file from which it can read
the connectivity and coordinates
  - a .rb file which specifies the rings, the pi bonds 
and the rigid bodies as a Python list. Pi bonds, if specified 
will be treated as rigid bodies. While the pi bonds or the 
rigid bodies are optional, the rings need to be specified.

The two options for dynamics are 
  IC - "Internal Coordinates" and
  TD -  "Torsional Dynamics" 
(TODO: to be passed as an argument to the integrator's constructor)

For ring assesment:
/share/apps/amber/14/AmberTools/bin/antechamber -i <ligand>.mol2 -fi mol2 -o <ligand>.sdf -fo sdf
python /home/lspirido/py_progs/32rigid_bodies.py <ligand>.sdf > <ligand>.rb


