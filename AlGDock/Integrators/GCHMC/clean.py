import os, sys
import glob

gmolDIR = "gmolmodel"
if not os.path.exists(gmolDIR):
  print "gmolmodel directory not found. Exit ..."
  sys.exit(1)

print "Cleaning objects... ",
objFNs = glob.glob("gmolmodel/*.o")
for objFN in objFNs:
  os.remove(os.path.abspath(objFN))
print "Done."

print "Cleaning libraries... ",
libgmolFNs  = glob.glob("libgmolmodel*.so*")
libboostFNs = glob.glob("libboost*.so*")
for libgmolFN in libgmolFNs:
  os.remove(os.path.abspath(libgmolFN))
for libboostFN in libboostFNs:
  os.remove(os.path.abspath(libboostFN))
print "Done."

print "Cleaning extensions... ",
extGCFNs = glob.glob("GCHMC*.so*")
for extGCFN in extGCFNs:
  os.remove(os.path.abspath(extGCFN))
print "Done."

print "Cleaning pdbs... ",
pdbFNs = glob.glob("pdbs/sb*.pdb")
for pdbFN in pdbFNs:
  os.remove(os.path.abspath(pdbFN))
if os.path.exists("pdbs"):
  os.rmdir("pdbs")
print "Done."

print "Cleaning logs... ",
logFNs  = [os.path.abspath(f) for f in glob.glob("gmolmodel/*.log")]
logFNs += [os.path.abspath(f) for f in glob.glob("*.log")]
for logFN in logFNs:
  os.remove(logFN)
print "Done."

print "Cleaning bjam files... ",
bjamFNs = ["Jamroot"]
for bjamFN in bjamFNs:
  os.remove(bjamFN)
print "Done."


