/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "SimTKmath.h"
#include <vector>
#include <exception>

using namespace SimTK;
using namespace std;

const Real TOL = 1e-10;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

template <class T>
void assertEqual(T val1, T val2) {
    ASSERT(abs(val1-val2) < TOL);
}

template <int N>
void assertEqual(Vec<N> val1, Vec<N> val2) {
    for (int i = 0; i < N; ++i)
        ASSERT(abs(val1[i]-val2[i]) < TOL);
}

void testHalfSpace() {
    ContactGeometry::HalfSpace hs;
    
    // Check intersections with various rays.
    
    Real distance;
    UnitVec3 normal;
    ASSERT(!hs.intersectsRay(Vec3(4, 0, 0), UnitVec3(1, 0, 0), distance, normal));
    ASSERT(!hs.intersectsRay(Vec3(-1, 0, 0), UnitVec3(-1, 1, 0), distance, normal));
    ASSERT(!hs.intersectsRay(Vec3(-1, 0, 0), UnitVec3(0, 1, 0), distance, normal));
    ASSERT(hs.intersectsRay(Vec3(-1, 0, 0), UnitVec3(1, 0, 0), distance, normal));
    assertEqual(1.0, distance);
    assertEqual(Vec3(-1, 0, 0), normal);
    ASSERT(hs.intersectsRay(Vec3(-2, 15, 37), UnitVec3(1, 0, 0), distance, normal));
    assertEqual(2.0, distance);
    assertEqual(Vec3(-1, 0, 0), normal);
    ASSERT(hs.intersectsRay(Vec3(-3, 1, 2), UnitVec3(1, 1, 1), distance, normal));
    assertEqual(3*Sqrt3, distance);
    assertEqual(Vec3(-1, 0, 0), normal);
    ASSERT(hs.intersectsRay(Vec3(2, 0, 0), UnitVec3(-1, 0, 1), distance, normal));
    assertEqual(2*Sqrt2, distance);
    assertEqual(Vec3(-1, 0, 0), normal);
    
    // Test finding the nearest point.
    
    Random::Gaussian random(0, 3);
    for (int i = 0; i < 100; i++) {
        Vec3 pos(random.getValue(), random.getValue(), random.getValue());
        bool inside;
        UnitVec3 normal;
        Vec3 nearest = hs.findNearestPoint(pos, inside, normal);
        assertEqual(nearest, Vec3(0, pos[1], pos[2]));
        ASSERT(inside == (pos[0] >= 0));
        assertEqual(normal, Vec3(-1, 0, 0));
    }
}

void testSphere() {
    // Create a sphere.
    
    Real radius = 3.5;
    ContactGeometry::Sphere sphere(radius);
    assert(sphere.getRadius() == radius);
    
    // Check intersections with various rays.
    
    Real distance;
    UnitVec3 normal;
    ASSERT(!sphere.intersectsRay(Vec3(4, 0, 0), UnitVec3(1, 0, 0), distance, normal));
    ASSERT(sphere.intersectsRay(Vec3(2, 0, 0), UnitVec3(1, 0, 0), distance, normal));
    assertEqual(1.5, distance);
    assertEqual(Vec3(1, 0, 0), normal);
    ASSERT(sphere.intersectsRay(Vec3(4, 0, 0), UnitVec3(-1, 0, 0), distance, normal));
    assertEqual(0.5, distance);
    assertEqual(Vec3(1, 0, 0), normal);
    ASSERT(sphere.intersectsRay(Vec3(2, 0, 0), UnitVec3(-1, 0, 0), distance, normal));
    assertEqual(5.5, distance);
    assertEqual(Vec3(-1, 0, 0), normal);
    ASSERT(sphere.intersectsRay(Vec3(0, 0, 0), UnitVec3(1, 1, 1), distance, normal));
    assertEqual(3.5, distance);
    assertEqual(Vec3(1.0/Sqrt3), normal);

    // Test finding the nearest point.
    
    Random::Gaussian random(0, 3);
    for (int i = 0; i < 100; i++) {
        Vec3 pos(random.getValue(), random.getValue(), random.getValue());
        bool inside;
        UnitVec3 normal;
        Vec3 nearest = sphere.findNearestPoint(pos, inside, normal);
        assertEqual(nearest, pos.normalize()*radius);
        ASSERT(inside == (pos.norm() <= radius));
        assertEqual(normal, pos.normalize());
    }
}

void testEllipsoid() {
    // Create a ellipsoid.

    Vec3 radii(1.5, 2.2, 3.1);
    ContactGeometry::Ellipsoid ellipsoid(radii);
    assert(ellipsoid.getRadii() == radii);

    // Check intersections with various rays.

    Real distance;
    UnitVec3 normal;
    ASSERT(!ellipsoid.intersectsRay(Vec3(4, 0, 0), UnitVec3(1, 0, 0), distance, normal));
    ASSERT(ellipsoid.intersectsRay(Vec3(1, 0, 0), UnitVec3(1, 0, 0), distance, normal));
    assertEqual(0.5, distance);
    assertEqual(Vec3(1, 0, 0), normal);
    ASSERT(ellipsoid.intersectsRay(Vec3(4, 0, 0), UnitVec3(-1, 0, 0), distance, normal));
    assertEqual(2.5, distance);
    assertEqual(Vec3(1, 0, 0), normal);
    ASSERT(ellipsoid.intersectsRay(Vec3(0, -5, 0), UnitVec3(0, 1, 0), distance, normal));
    assertEqual(2.8, distance);
    assertEqual(Vec3(0, -1, 0), normal);
    ASSERT(ellipsoid.intersectsRay(Vec3(0, 0, 0), UnitVec3(1, 1, 1), distance, normal));
    assertEqual(sqrt(3/(1/(radii[0]*radii[0])+1/(radii[1]*radii[1])+1/(radii[2]*radii[2]))), distance);
    assertEqual(UnitVec3(1/(radii[0]*radii[0]), 1/(radii[1]*radii[1]), 1/(radii[2]*radii[2])), normal);

    // Test finding the nearest point.

    Random::Gaussian random(0, 2);
    for (int i = 0; i < 100; i++) {
        Vec3 pos(random.getValue(), random.getValue(), random.getValue());
        bool inside;
        UnitVec3 normal;
        Vec3 nearest = ellipsoid.findNearestPoint(pos, inside, normal);
        assertEqual(nearest[0]*nearest[0]/(radii[0]*radii[0])+nearest[1]*nearest[1]/(radii[1]*radii[1])+nearest[2]*nearest[2]/(radii[2]*radii[2]), 1.0);
        Real projectedRadius = pos[0]*pos[0]/(radii[0]*radii[0])+pos[1]*pos[1]/(radii[1]*radii[1])+pos[2]*pos[2]/(radii[2]*radii[2]);
        ASSERT(inside == (projectedRadius < 1.0));
        Vec3 projectedPoint = pos/sqrt(projectedRadius);
        ASSERT((nearest-pos).normSqr() < (projectedPoint-pos).normSqr());
        assertEqual(normal, UnitVec3(nearest[0]/(radii[0]*radii[0]), nearest[1]/(radii[1]*radii[1]), nearest[2]/(radii[2]*radii[2])));
    }
}

int main() {
    try {
        testHalfSpace();
        testSphere();
        testEllipsoid();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
