#ifndef SimTK_SIMMATH_CONTACT_TRACKER_SUBSYSTEM_H_
#define SimTK_SIMMATH_CONTACT_TRACKER_SUBSYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman, Peter Eastman                                    *
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

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Contact.h"

namespace SimTK {

//==============================================================================
//                              CONTACT TRACKER
//==============================================================================
/** A ContactTracker implements an algorithm for detecting overlaps or 
potential overlaps between pairs of ContactGeometry objects, and managing
Contact objects that track individual contacts as they evolve through
time. This class is used internally by ContractTrackerSubsystem and there 
usually is no reason to access it directly. The exception is if you are 
defining a new ContactGeometry subclass. In that case, you will also need to 
define one or more ContactTrackers to detect collisions with your new geometry
type, then register it with the ContactTrackerSubsystem. 

The result of a ContactTracker when applied to a pair of contact
surfaces, is either a determination that the surfaces are not in contact,
or a Contact object describing their contact interaction. There are different
types of these Contact objects (for example, PointContact, LineContact, 
MeshContact) and the same algorithm may result in different kinds of Contact 
under different circumstances. At each evaluation, the caller passes in the 
previous Contact object, if any, that was associated with two ContactSurfaces,
then receives an update from the algorithm.

Note that ContactTrackers that manage dissimilar geometry type pairs expect 
the two types in a particular order, e.g. (halfspace,sphere) rather than
(sphere,halfspace) but are used for all contacts involving that pair of 
types. It is up to the ContactTrackerSubsystem to ensure that the contact
surfaces are presented in the correct order regardless of how they are 
encountered. The Contact objects that are created and managed by trackers
always have their (surface1,surface2) pairs in the order required by the
tracker that handles those types. **/
class SimTK_SIMMATH_EXPORT ContactTracker {
public:
class HalfSpaceSphere;
class HalfSpaceEllipsoid;
class SphereSphere;
class HalfSpaceTriangleMesh;
class SphereTriangleMesh;
class TriangleMeshTriangleMesh;
class ConvexImplicitPair;
class GeneralImplicitPair;

/** Base class constructor for use by the concrete classes. **/
ContactTracker(ContactGeometryTypeId typeOfSurface1,
               ContactGeometryTypeId typeOfSurface2)
:   m_surfaceTypes(typeOfSurface1, typeOfSurface2)
{
}

/** Return the pair of contact geometry type ids handled by this tracker,
in the order that they must be presented to the tracker's methods. **/
const std::pair<ContactGeometryTypeId,ContactGeometryTypeId>&
getContactGeometryTypeIds() const {return m_surfaceTypes;}

virtual ~ContactTracker() {}

/** The ContactTrackerSubsystem will invoke this method for any pair of
contact surfaces that is already being tracked, or for which the static broad 
phase analysis indicated that they might be in contact now. Only position 
information is available. Note that the arguments and Contact object surfaces
must be ordered by geometry type id as required by this tracker. **/
virtual bool trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& surface1,
    const Transform& X_GS2, 
    const ContactGeometry& surface2,
    Real                   cutoff,
    Contact&               currentStatus) const = 0;

/** The ContactTrackerSubsystem will invoke this method for any tracked pair of
contact surfaces that is still not in contact after trackContact() looked at
it, or any untracked pair for which the dynamic broad phase indicated that 
they might be in contact within the interval of interest. Position, velocity, 
and acceleration information may be used. Ordering must be correct as 
discussed for trackContact(). **/
virtual bool predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const = 0;

/** At the beginning of a simulation we will have no past information to
help disambiguate tricky contact situations. This method may use current
position and velocity information in heuristics for guessing the contact
status between the indicated pair of surfaces. Ordering must be correct as 
discussed for trackContact(). **/
virtual bool initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const = 0;

/** Given two shapes for which implicit functions are known, and a rough-guess
contact point for each shape, refine those contact points to obtain the nearest
pair that satisfies contact conditions to a requested accuracy. For separated
objects, these will be the points of closest approach between the surfaces; for
contacting objects these are the points of maximum penetration.

In implicit form there are six unknowns (spatial coordinates of the contact 
points). The six contact conditions we use are:
  - 2 equations: Point P is on the surface of shape A and point Q is on the 
    surface of shape B (that is, the implicit functions both return zero).
  - 2 equations: The normal vector at point P is perpendicular to the
    tangent plane at point Q.
  - 2 equations: The separation vector (P-Q) is zero or perpendicular to the
    tangent plane at P.

Note that these equations could be satisfied by incorrect points that have the
opposite normals because the perpendicularity conditions can't distinguish n
from -n. We are depending on having an initial guess that is good enough so
that we find the correct solution by going downhill from there. Don't try to
use this if you don't have a reasonably good guess already. 

@returns \c true if the requested accuracy is achieved but returns its best
attempt at the refined points regardless. **/
static bool refineImplicitPair
   (const ContactGeometry& shapeA, Vec3& pointP,    // in/out
    const ContactGeometry& shapeB, Vec3& pointQ,    // in/out
    const Transform& X_AB, Real accuracyRequested,
    Real& accuracyAchieved, int& numIterations);

/** Use Minkowski Portal Refinement (XenoCollide method by G. Snethen) to
generate a reasonably good starting estimate of the contact points between
two implicit shapes that are in contact. MPR cannot find those points if the
surfaces are separated. Returns \c true if a contact was found; otherwise
the returned points are meaningless. **/
static bool estimateImplicitPairContactUsingMPR
   (const ContactGeometry& shapeA, const ContactGeometry& shapeB, 
    const Transform& X_AB, 
    Vec3& pointP, Vec3& pointQ, int& numIterations);

static Vec6 findImplicitPairError
   (const ContactGeometry& shapeA, const Vec3& pointP,
    const ContactGeometry& shapeB, const Vec3& pointQ,
    const Transform& X_AB);

static Mat66 calcImplicitPairJacobian
   (const ContactGeometry& shapeA, const Vec3& pointP,
    const ContactGeometry& shapeB, const Vec3& pointQ,
    const Transform& X_AB, const Vec6& err0);

//--------------------------------------------------------------------------
                                private:
// This tracker should be called only for surfaces of these two types,
// in this order.
std::pair<ContactGeometryTypeId,ContactGeometryTypeId> m_surfaceTypes;
};



//==============================================================================
//                     HALFSPACE-SPHERE CONTACT TRACKER
//==============================================================================
/** This ContactTracker handles contacts between a ContactGeometry::HalfSpace
and a ContactGeometry::Sphere, in that order. **/
class SimTK_SIMMATH_EXPORT ContactTracker::HalfSpaceSphere 
:   public ContactTracker {
public:
HalfSpaceSphere() 
:   ContactTracker(ContactGeometry::HalfSpace::classTypeId(),
                   ContactGeometry::Sphere::classTypeId()) {}

virtual ~HalfSpaceSphere() {}

virtual bool trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& surface1,
    const Transform& X_GS2, 
    const ContactGeometry& surface2,
    Real                   cutoff,
    Contact&               currentStatus) const;

virtual bool predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const;

virtual bool initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const;
};



//==============================================================================
//                     HALFSPACE-ELLIPSOID CONTACT TRACKER
//==============================================================================
/** This ContactTracker handles contacts between a ContactGeometry::HalfSpace
and a ContactGeometry::Ellipsoid, in that order. **/
class SimTK_SIMMATH_EXPORT ContactTracker::HalfSpaceEllipsoid 
:   public ContactTracker {
public:
HalfSpaceEllipsoid() 
:   ContactTracker(ContactGeometry::HalfSpace::classTypeId(),
                   ContactGeometry::Ellipsoid::classTypeId()) {}

virtual ~HalfSpaceEllipsoid() {}

virtual bool trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& surface1,
    const Transform& X_GS2, 
    const ContactGeometry& surface2,
    Real                   cutoff,
    Contact&               currentStatus) const;

virtual bool predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const;

virtual bool initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const;
};



//==============================================================================
//                       SPHERE-SPHERE CONTACT TRACKER
//==============================================================================
/** This ContactTracker handles contacts between two ContactGeometry::Sphere
objects. **/
class SimTK_SIMMATH_EXPORT ContactTracker::SphereSphere 
:   public ContactTracker {
public:
SphereSphere() 
:   ContactTracker(ContactGeometry::Sphere::classTypeId(),
                   ContactGeometry::Sphere::classTypeId()) {}

virtual ~SphereSphere() {}

virtual bool trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& surface1,
    const Transform& X_GS2, 
    const ContactGeometry& surface2,
    Real                   cutoff,
    Contact&               currentStatus) const;

virtual bool predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const;

virtual bool initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const;
};



//==============================================================================
//                 HALFSPACE-TRIANGLE MESH CONTACT TRACKER
//==============================================================================
/** This ContactTracker handles contacts between a ContactGeometry::HalfSpace
and a ContactGeometry::TriangleMesh, in that order. **/
class SimTK_SIMMATH_EXPORT ContactTracker::HalfSpaceTriangleMesh
:   public ContactTracker {
public:
HalfSpaceTriangleMesh() 
:   ContactTracker(ContactGeometry::HalfSpace::classTypeId(),
                   ContactGeometry::TriangleMesh::classTypeId()) {}

virtual ~HalfSpaceTriangleMesh() {}

virtual bool trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& surface1,    // the half space
    const Transform& X_GS2, 
    const ContactGeometry& surface2,    // the mesh
    Real                   cutoff,
    Contact&               currentStatus) const;

virtual bool predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const;

virtual bool initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const;

private:
void processBox(const ContactGeometry::TriangleMesh&              mesh, 
                const ContactGeometry::TriangleMesh::OBBTreeNode& node, 
                const Transform& X_HM, const UnitVec3& hsNormal_M, 
                Real hsFaceHeight_M, std::set<int>& insideFaces) const;
void addAllTriangles(const ContactGeometry::TriangleMesh::OBBTreeNode& node, 
                     std::set<int>& insideFaces) const; 
};



//==============================================================================
//                 SPHERE - TRIANGLE MESH CONTACT TRACKER
//==============================================================================
/** This ContactTracker handles contacts between a ContactGeometry::Sphere
and a ContactGeometry::TriangleMesh, in that order. **/
class SimTK_SIMMATH_EXPORT ContactTracker::SphereTriangleMesh
:   public ContactTracker {
public:
SphereTriangleMesh() 
:   ContactTracker(ContactGeometry::Sphere::classTypeId(),
                   ContactGeometry::TriangleMesh::classTypeId()) {}

virtual ~SphereTriangleMesh() {}

virtual bool trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& surface1,    // the sphere
    const Transform& X_GS2, 
    const ContactGeometry& surface2,    // the mesh
    Real                   cutoff,
    Contact&               currentStatus) const;

virtual bool predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const;

virtual bool initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const;

private:
void processBox
   (const ContactGeometry::TriangleMesh&              mesh, 
    const ContactGeometry::TriangleMesh::OBBTreeNode& node, 
    const Vec3& center_M, Real radius2,   
    std::set<int>& insideFaces) const ;
};



//==============================================================================
//             TRIANGLE MESH - TRIANGLE MESH CONTACT TRACKER
//==============================================================================
/** This ContactTracker handles contacts between two 
ContactGeometry::TriangleMesh surfaces. **/
class SimTK_SIMMATH_EXPORT ContactTracker::TriangleMeshTriangleMesh
:   public ContactTracker {
public:
TriangleMeshTriangleMesh() 
:   ContactTracker(ContactGeometry::TriangleMesh::classTypeId(),
                   ContactGeometry::TriangleMesh::classTypeId()) {}

virtual ~TriangleMeshTriangleMesh() {}

virtual bool trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& surface1,    // mesh1
    const Transform& X_GS2, 
    const ContactGeometry& surface2,    // mesh2
    Real                   cutoff,
    Contact&               currentStatus) const;

virtual bool predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const;

virtual bool initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const;

private:
void findIntersectingFaces
   (const ContactGeometry::TriangleMesh&                mesh1, 
    const ContactGeometry::TriangleMesh&                mesh2,
    const ContactGeometry::TriangleMesh::OBBTreeNode&   node1, 
    const ContactGeometry::TriangleMesh::OBBTreeNode&   node2, 
    const OrientedBoundingBox&                          node2Bounds_M1,
    const Transform&                                    X_M1M2, 
    std::set<int>&                                      insideFaces1, 
    std::set<int>&                                      insideFaces2) const; 
void findBuriedFaces
   (const ContactGeometry::TriangleMesh&    mesh,
    const ContactGeometry::TriangleMesh&    otherMesh,
    const Transform&                        X_OM, 
    std::set<int>&                          insideFaces) const;
void tagFaces(const ContactGeometry::TriangleMesh&   mesh, 
              Array_<int>&                           faceType,
              std::set<int>&                         triangles, 
              int                                    index,
              int                                    depth) const;
};


//==============================================================================
//               CONVEX IMPLICIT SURFACE PAIR CONTACT TRACKER
//==============================================================================
/** This ContactTracker handles contacts between two smooth, convex objects
by using their implicit functions. Create one of these for each possible
pair that you want handled this way. **/
class SimTK_SIMMATH_EXPORT ContactTracker::ConvexImplicitPair 
:   public ContactTracker {
public:
ConvexImplicitPair(ContactGeometryTypeId type1, ContactGeometryTypeId type2) 
:   ContactTracker(type1, type2) {}

virtual ~ConvexImplicitPair() {}

virtual bool trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& surface1,
    const Transform& X_GS2, 
    const ContactGeometry& surface2,
    Real                   cutoff,
    Contact&               currentStatus) const;

virtual bool predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const;

virtual bool initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const;
};


//==============================================================================
//                GENERAL IMPLICIT SURFACE PAIR CONTACT TRACKER
//==============================================================================
/** This ContactTracker handles contacts between two arbitrary smooth surfaces
by using their implicit functions, with no shape restrictions. Each surface
must provide a bounding hierarchy with "safe" leaf objects, meaning that
interactions between a leaf of each surface yield at most one solution.

Create one of these for each possible pair that you want handled this way. **/
class SimTK_SIMMATH_EXPORT ContactTracker::GeneralImplicitPair 
:   public ContactTracker {
public:
GeneralImplicitPair(ContactGeometryTypeId type1, ContactGeometryTypeId type2) 
:   ContactTracker(type1, type2) {}

virtual ~GeneralImplicitPair() {}

virtual bool trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& surface1,
    const Transform& X_GS2, 
    const ContactGeometry& surface2,
    Real                   cutoff,
    Contact&               currentStatus) const;

virtual bool predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const;

virtual bool initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const;
};

} // namespace SimTK

#endif // SimTK_SIMMATH_CONTACT_TRACKER_SUBSYSTEM_H_
