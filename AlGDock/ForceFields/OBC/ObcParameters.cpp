
/* Portions copyright (c) 2006-2009 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <math.h>
#include <sstream>
#include <string.h>

//  #include "openmm/OpenMMException.h"
#include "ObcParameters.h"

using std::vector;
//  using namespace OpenMM;

/**---------------------------------------------------------------------------------------

    ObcParameters constructor 

    @param numberOfAtoms       number of atoms
    @param obcType             OBC type (Eq. 7 or 8 in paper)

    --------------------------------------------------------------------------------------- */

ObcParameters::ObcParameters(int numberOfAtoms, ObcParameters::ObcType obcType) :
                              _numberOfAtoms(numberOfAtoms),
                              _solventDielectric(78.3),
                              _soluteDielectric(1.0),
                              _electricConstant(-0.5*ONE_4PI_EPS0),
                              _probeRadius(0.14),
                              _pi4Asolv(28.3919551),
                              _dielectricOffset(0.009),
                              _obcType(obcType),
                              _cutoff(false),
                              _periodic(false) {

    _atomicRadii.resize(numberOfAtoms);
    _scaledRadiusFactors.resize(numberOfAtoms);
    _partialCharges.resize(numberOfAtoms);
    setObcTypeParameters(obcType);

}

/**---------------------------------------------------------------------------------------

    ObcParameters destructor 

    --------------------------------------------------------------------------------------- */

ObcParameters::~ObcParameters() {
}

/**---------------------------------------------------------------------------------------

   Get number of atoms

   @return number of atoms

   --------------------------------------------------------------------------------------- */

int ObcParameters::getNumberOfAtoms() const {
    return _numberOfAtoms;
}

/**---------------------------------------------------------------------------------------

   Set number of atoms

   @param numberOfAtoms number of atoms

   --------------------------------------------------------------------------------------- */

void ObcParameters::setNumberOfAtoms(int numberOfAtoms) {
    _numberOfAtoms = numberOfAtoms;
}

/**---------------------------------------------------------------------------------------

   Get strength

   @return strength

   --------------------------------------------------------------------------------------- */

double ObcParameters::getStrength() const {
    return _strength;
}

/**---------------------------------------------------------------------------------------

   Set number of atoms

   @param numberOfAtoms number of atoms

   --------------------------------------------------------------------------------------- */

void ObcParameters::setStrength(double strength) {
    _strength = strength;
}


/**---------------------------------------------------------------------------------------

    Get OBC type

    @return OBC type

    --------------------------------------------------------------------------------------- */

ObcParameters::ObcType ObcParameters::getObcType() const {
    return _obcType;
}

/**---------------------------------------------------------------------------------------

    Set OBC type specific parameters

    @param obcType OBC type (ObcTypeI or ObcTypeII -- Eq. 7 or 8)

    --------------------------------------------------------------------------------------- */

void ObcParameters::setObcTypeParameters(ObcParameters::ObcType obcType) {
    if (obcType == ObcTypeI) {
        _alphaObc   = 0.8f;
        _betaObc    = 0.0f;
        _gammaObc   = 2.91f;
    } else {
        _alphaObc   = 1.0f;
        _betaObc    = 0.8f;
        _gammaObc   = 4.85f;
    }
    _obcType = obcType;
}

/**---------------------------------------------------------------------------------------

    Get dielectricOffset

    @return _dielectricOffset

    --------------------------------------------------------------------------------------- */

double ObcParameters::getDielectricOffset() const {
    return _dielectricOffset;
}

/**---------------------------------------------------------------------------------------

    Get alpha OBC (Eqs. 6 & 7) in Proteins paper

    @return alphaObc

    --------------------------------------------------------------------------------------- */

double ObcParameters::getAlphaObc() const {
    return _alphaObc;
}

/**---------------------------------------------------------------------------------------

    Get beta OBC (Eqs. 6 & 7) in Proteins paper

    @return betaObc

    --------------------------------------------------------------------------------------- */

double ObcParameters::getBetaObc() const {
    return _betaObc;
}

/**---------------------------------------------------------------------------------------

    Get gamma OBC (Eqs. 6 & 7) in Proteins paper

    @return gammaObc

    --------------------------------------------------------------------------------------- */

double ObcParameters::getGammaObc() const {
    return _gammaObc;
}

/**---------------------------------------------------------------------------------------

   Get solvent dielectric

   @return solvent dielectric

   --------------------------------------------------------------------------------------- */

double ObcParameters::getSolventDielectric() const {
    return _solventDielectric;
}

/**---------------------------------------------------------------------------------------

   Set solvent dielectric

   @param solventDielectric solvent dielectric

   --------------------------------------------------------------------------------------- */

void ObcParameters::setSolventDielectric(double solventDielectric) {
    _solventDielectric = solventDielectric;
}
/**---------------------------------------------------------------------------------------

   Get solute dielectric

   @return soluteDielectric

   --------------------------------------------------------------------------------------- */

double ObcParameters::getSoluteDielectric() const {
    return _soluteDielectric;
}

/**---------------------------------------------------------------------------------------

   Set solute dielectric

   @param soluteDielectric solute dielectric

   --------------------------------------------------------------------------------------- */

void ObcParameters::setSoluteDielectric(double soluteDielectric) {
    _soluteDielectric = soluteDielectric;
}

/**---------------------------------------------------------------------------------------

   Get electric constant 

   @return electricConstant

   --------------------------------------------------------------------------------------- */

double ObcParameters::getElectricConstant() const {
    return _electricConstant;
}

/**---------------------------------------------------------------------------------------

   Get probe radius 

   @return probeRadius

   --------------------------------------------------------------------------------------- */

double ObcParameters::getProbeRadius() const {
    return _probeRadius;
}

/**---------------------------------------------------------------------------------------

   Set probe radius  

   @param probeRadius   probe radius

   --------------------------------------------------------------------------------------- */

void ObcParameters::setProbeRadius(double probeRadius) {
    _probeRadius = probeRadius;
}

/**---------------------------------------------------------------------------------------

   Get pi*4*Asolv:  used in ACE approximation for nonpolar term  
         ((double) M_PI)*4.0f*0.0049*1000.0; (Still) 
         ((double) M_PI)*4.0f*0.0054*1000.0; (OBC) 

   @return pi4Asolv

   --------------------------------------------------------------------------------------- */

double ObcParameters::getPi4Asolv() const {
    return _pi4Asolv;
}

void ObcParameters::setPi4Asolv(double pi4Asolv) {
    _pi4Asolv = pi4Asolv;
}

/**---------------------------------------------------------------------------------------

    Get AtomicRadii array

    @return array of atomic radii

    --------------------------------------------------------------------------------------- */

const vector<double>& ObcParameters::getAtomicRadii() const {
    return _atomicRadii;
}

/**---------------------------------------------------------------------------------------

    Set AtomicRadii array

    @param atomicRadii vector of atomic radii

    --------------------------------------------------------------------------------------- */

void ObcParameters::setAtomicRadii(const vector<double>& atomicRadii) {

    if (atomicRadii.size() == _atomicRadii.size()) {
        for (unsigned int ii = 0; ii < atomicRadii.size(); ii++) {
            _atomicRadii[ii] = atomicRadii[ii];
        }   
    } else {
        std::stringstream msg;
        msg << "ObcParameters: input size for atomic radii does not agree w/ current size: input=";
        msg << atomicRadii.size();
        msg << " current size=" << _atomicRadii.size();
    }

}

/**---------------------------------------------------------------------------------------

    Return OBC scale factors

    @return array 

    --------------------------------------------------------------------------------------- */

const vector<double>& ObcParameters::getScaledRadiusFactors() const {
    return _scaledRadiusFactors;
}

/**---------------------------------------------------------------------------------------

    Set OBC scale factors

    @param scaledRadiusFactors  scaledRadiusFactors

    --------------------------------------------------------------------------------------- */

void ObcParameters::setScaledRadiusFactors(const vector<double>& scaledRadiusFactors) {

    if (scaledRadiusFactors.size() == _scaledRadiusFactors.size()) {
        for (unsigned int ii = 0; ii < scaledRadiusFactors.size(); ii++) {
            _scaledRadiusFactors[ii] = scaledRadiusFactors[ii];
        }   
    } else {
        std::stringstream msg;
        msg << "ObcParameters: input size for scaled radius factors does not agree w/ current size: input=";
        msg << scaledRadiusFactors.size();
        msg << " current size=" << _scaledRadiusFactors.size();
//        throw OpenMM::OpenMMException(msg.str());
    }

}

/**---------------------------------------------------------------------------------------

    Get AtomicRadii array

    @return array of atomic radii

    --------------------------------------------------------------------------------------- */

const vector<double>& ObcParameters::getPartialCharges() const {
    return _partialCharges;
}

/**---------------------------------------------------------------------------------------

    Set AtomicRadii array

    @param atomicRadii vector of atomic radii

    --------------------------------------------------------------------------------------- */

void ObcParameters::setPartialCharges(const vector<double>& partialCharges) {

    if (partialCharges.size() == _partialCharges.size()) {
        for (unsigned int ii = 0; ii < partialCharges.size(); ii++) {
            _partialCharges[ii] = partialCharges[ii];
        }   
    } else {
        std::stringstream msg;
        msg << "ObcParameters: input size for partial charges does not agree w/ current size: input=";
        msg << partialCharges.size();
        msg << " current size=" << _partialCharges.size();
//        throw OpenMM::OpenMMException(msg.str());
    }

}

/**---------------------------------------------------------------------------------------

      Set the force to use a cutoff.

      @param distance            the cutoff distance

      --------------------------------------------------------------------------------------- */

void ObcParameters::setUseCutoff(double distance) {

     _cutoff         = true;
     _cutoffDistance = distance;
}

/**---------------------------------------------------------------------------------------

      Get whether to use a cutoff.

      --------------------------------------------------------------------------------------- */

bool ObcParameters::getUseCutoff() const {
     return _cutoff;
}

/**---------------------------------------------------------------------------------------

      Get the cutoff distance.

      --------------------------------------------------------------------------------------- */

double ObcParameters::getCutoffDistance() const {
     return _cutoffDistance;
}

/**---------------------------------------------------------------------------------------

      Set the force to use periodic boundary conditions.  This requires that a cutoff has
      also been set, and the smallest side of the periodic box is at least twice the cutoff
      distance.

      @param vectors    the vectors defining the periodic box

      --------------------------------------------------------------------------------------- */

void ObcParameters::setPeriodic(OpenMM::RealVec* vectors) {

    assert(_cutoff);

    assert(vectors[0][0] >= 2.0*_cutoffDistance);
    assert(vectors[1][1] >= 2.0*_cutoffDistance);
    assert(vectors[2][2] >= 2.0*_cutoffDistance);

    _periodic           = true;
    _periodicBoxVectors[0] = vectors[0];
    _periodicBoxVectors[1] = vectors[1];
    _periodicBoxVectors[2] = vectors[2];
}

/**---------------------------------------------------------------------------------------

      Get whether to use periodic boundary conditions.

      --------------------------------------------------------------------------------------- */

bool ObcParameters::getPeriodic() {
     return _periodic;
}

/**---------------------------------------------------------------------------------------

      Get the periodic box dimension

      --------------------------------------------------------------------------------------- */

const OpenMM::RealVec* ObcParameters::getPeriodicBox() {
     return _periodicBoxVectors;
}

//  extern "C" double call_C_f(C* p, int i) // wrapper function
//  {
//      return p->f(i);
//  }
