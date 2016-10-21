
/* Portions copyright (c) 2006-2013 Stanford University and Simbios.
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

#include <cstring>
#include <sstream>

//  #include "SimTKOpenMMUtilities.h"
#include "ReferenceForce.h"

#include <cstdio>

using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceForce constructor

   --------------------------------------------------------------------------------------- */

ReferenceForce::ReferenceForce() {
}

/**---------------------------------------------------------------------------------------

   ReferenceForce destructor

   --------------------------------------------------------------------------------------- */

ReferenceForce::~ReferenceForce() {
}

/**---------------------------------------------------------------------------------------

   Given two coordinates on a periodic lattice, return the difference between them.

   --------------------------------------------------------------------------------------- */

double ReferenceForce::periodicDifference(double val1, double val2, double period) {
    double diff = val1-val2;
    double base = (double) (floor(diff/period+0.5)*period);
    return diff-base;

}

//void ReferenceForce::getDeltaR(const RealVec& atomCoordinatesI, const RealVec& atomCoordinatesJ,
//                               double* deltaR) {
void ReferenceForce::getDeltaR(const vector3& atomCoordinatesI, const vector3& atomCoordinatesJ,
                               double* deltaR) {

   deltaR[XIndex]    = atomCoordinatesJ[0] - atomCoordinatesI[0];
   deltaR[YIndex]    = atomCoordinatesJ[1] - atomCoordinatesI[1];
   deltaR[ZIndex]    = atomCoordinatesJ[2] - atomCoordinatesI[2];

   deltaR[R2Index]   = DOT3(deltaR, deltaR);
   deltaR[RIndex]    = (double) SQRT(deltaR[R2Index]);
}


// Removed getDeltaRPeriodic

//  double* ReferenceForce::getVariablePointer(Lepton::CompiledExpression& expression, const std::string& name) {
//      if (expression.getVariables().find(name) == expression.getVariables().end())
//          return NULL;
//      return &expression.getVariableReference(name);
//  }

void ReferenceForce::setVariable(double* pointer, double value) {
    if (pointer != NULL)
        *pointer = value;
}
