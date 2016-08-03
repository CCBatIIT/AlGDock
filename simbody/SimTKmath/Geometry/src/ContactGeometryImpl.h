#ifndef SimTK_SIMMATH_CONTACT_GEOMETRY_IMPL_H_
#define SimTK_SIMMATH_CONTACT_GEOMETRY_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
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


#include "simmath/internal/Geo.h"
#include "simmath/internal/Geo_Sphere.h"
#include "simmath/internal/OBBTree.h"
#include "simmath/internal/ContactGeometry.h"

#include <limits>

namespace SimTK {

//==============================================================================
//                             CONTACT GEOMETRY IMPL
//==============================================================================
class SimTK_SIMMATH_EXPORT ContactGeometryImpl {
public:
    ContactGeometryImpl() : myHandle(0) {}
    virtual ~ContactGeometryImpl() {
        clearMyHandle();
    }

    /* Create a new ContactGeometryTypeId and return this unique integer 
    (thread safe). Each distinct type of ContactGeometry should use this to
    initialize a static variable for that concrete class. */
    static ContactGeometryTypeId  createNewContactGeometryTypeId()
    {   static AtomicInteger nextAvailableId = 1;
        return ContactGeometryTypeId(nextAvailableId++); }

    virtual ContactGeometryImpl*  clone() const = 0;
    virtual ContactGeometryTypeId getTypeId() const = 0;

    virtual Vec3 findNearestPoint(const Vec3& position, bool& inside, 
                                  UnitVec3& normal) const = 0;
    virtual bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                               Real& distance, UnitVec3& normal) const = 0;

    virtual void getBoundingSphere(Vec3& center, Real& radius) const = 0;


    virtual bool isSmooth() const = 0;
    virtual bool isConvex() const = 0;
    virtual bool isFinite() const = 0;

    // Smooth surfaces only.
    virtual void calcCurvature(const Vec3& point, Vec2& curvature, 
                       Rotation& orientation) const
    {   SimTK_THROW2(Exception::UnimplementedVirtualMethod, 
        "ContactGeometryImpl", "calcCurvature()"); }

    // Smooth surfaces only.
    virtual const Function& getImplicitFunction() const
    {   SimTK_THROW2(Exception::UnimplementedVirtualMethod, 
        "ContactGeometryImpl", "getImplicitFunction()"); }

    // Convex surfaces only.
    virtual Vec3 calcSupportPoint(UnitVec3 direction) const
    {   SimTK_THROW2(Exception::UnimplementedVirtualMethod, 
        "ContactGeometryImpl", "calcSupportPoint()"); }

    const OBBTree& getOBBTree() const {return obbTree;}

    ContactGeometry* getMyHandle() {return myHandle;}
    void setMyHandle(ContactGeometry& h) {myHandle = &h;}
    void clearMyHandle() {myHandle = 0;}
protected:
    ContactGeometry*        myHandle;
    OBBTree                 obbTree;

};



//==============================================================================
//                             HALF SPACE IMPL
//==============================================================================
class HalfSpaceImplicitFunction : public Function {
public:
    HalfSpaceImplicitFunction() : ownerp(0) {}
    HalfSpaceImplicitFunction(const ContactGeometry::HalfSpace::Impl& owner) 
    :   ownerp(&owner) {}
    void setOwner(const ContactGeometry::HalfSpace::Impl& owner) {ownerp=&owner;}

    // Value is positive for x>0.
    Real calcValue(const Vector& P) const {return P[0];}
    // First derivative w.r.t. x is 1, all else is zero.
    Real calcDerivative(const Array_<int>& components, 
                        const Vector& P) const
    {   if (components.empty()) return calcValue(P);
        if (components.size()==1 && components[0]==0) return 1;
        return 0; }

    int getArgumentSize() const {return 3;}
    int getMaxDerivativeOrder() const
    {   return std::numeric_limits<int>::max(); }
private:
    const ContactGeometry::HalfSpace::Impl* ownerp; // just a ref.; don't delete
};


class ContactGeometry::HalfSpace::Impl : public ContactGeometryImpl {
public:
    Impl() : ContactGeometryImpl() {
    }
    ContactGeometryImpl* clone() const {
        return new Impl();
    }

    ContactGeometryTypeId getTypeId() const {return classTypeId();}

    Vec3 findNearestPoint(const Vec3& position, bool& inside, 
                          UnitVec3& normal) const;
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                       Real& distance, UnitVec3& normal) const;
    void getBoundingSphere(Vec3& center, Real& radius) const;

    bool isSmooth() const {return true;}
    bool isConvex() const {return false;}
    bool isFinite() const {return false;}

    // Curvature is zero everywhere. Since the half plane occupies x>0 in
    // its own frame, the surface normal is -x, and -x,y,-z forms a right 
    // handed set.
    void calcCurvature(const Vec3& point, Vec2& curvature, 
                       Rotation& orientation) const
    {   curvature = 0;
        orientation.setRotationFromUnitVecsTrustMe
            (UnitVec3(-XAxis), UnitVec3(YAxis), UnitVec3(-ZAxis));
    }

    const Function& getImplicitFunction() const {return function;}

    static ContactGeometryTypeId classTypeId() {
        static const ContactGeometryTypeId id = 
            createNewContactGeometryTypeId();
        return id;
    }
private:
    HalfSpaceImplicitFunction function;
};



//==============================================================================
//                                SPHERE IMPL
//==============================================================================
class SphereImplicitFunction : public Function {
public:
    SphereImplicitFunction() : ownerp(0) {}
    SphereImplicitFunction(const ContactGeometry::Sphere::Impl& owner) 
    :   ownerp(&owner) {}
    void setOwner(const ContactGeometry::Sphere::Impl& owner) {ownerp=&owner;}
    Real calcValue(const Vector& x) const;
    Real calcDerivative(const Array_<int>& derivComponents, 
                        const Vector& x) const;
    int getArgumentSize() const {return 3;}
    int getMaxDerivativeOrder() const
    {   return std::numeric_limits<int>::max(); }
private:
    const ContactGeometry::Sphere::Impl* ownerp; // just a reference; don't delete
};

class ContactGeometry::Sphere::Impl : public ContactGeometryImpl {
public:
    explicit Impl(Real radius) : radius(radius) {
        function.setOwner(*this);
        createOBBTree(); 
    }

    ContactGeometryImpl* clone() const {
        return new Impl(radius);
    }
    Real getRadius() const {
        return radius;
    }
    void setRadius(Real r) {
        radius = r;
    }

    ContactGeometryTypeId getTypeId() const {return classTypeId();}

    Vec3 findNearestPoint(const Vec3& position, bool& inside, 
                          UnitVec3& normal) const;
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                       Real& distance, UnitVec3& normal) const;
    void getBoundingSphere(Vec3& center, Real& radius) const;

    bool isSmooth() const {return true;}
    bool isConvex() const {return true;}
    bool isFinite() const {return true;}

    Vec3 calcSupportPoint(UnitVec3 direction) const {
        return radius*direction;
    }
    void calcCurvature(const Vec3& point, Vec2& curvature, 
                       Rotation& orientation) const;
    const Function& getImplicitFunction() const {
        return function;
    }

    static ContactGeometryTypeId classTypeId() {
        static const ContactGeometryTypeId id = 
            createNewContactGeometryTypeId();
        return id;
    }
private:
    void createOBBTree();

    Real                    radius;
    SphereImplicitFunction  function;
};



//==============================================================================
//                                ELLIPSOID IMPL
//==============================================================================
class EllipsoidImplicitFunction : public Function {
public:
    EllipsoidImplicitFunction() : ownerp(0) {}
    EllipsoidImplicitFunction(const ContactGeometry::Ellipsoid::Impl& owner) 
    :   ownerp(&owner) {}
    void setOwner(const ContactGeometry::Ellipsoid::Impl& owner) {ownerp=&owner;}
    Real calcValue(const Vector& x) const;
    Real calcDerivative(const Array_<int>& derivComponents, 
                        const Vector& x) const;
    int getArgumentSize() const {return 3;}
    int getMaxDerivativeOrder() const
    {   return std::numeric_limits<int>::max(); }
private:
    const ContactGeometry::Ellipsoid::Impl* ownerp;// just a ref.; don't delete
};

class ContactGeometry::Ellipsoid::Impl : public ContactGeometryImpl  {
public:
    explicit Impl(const Vec3& radii)
    :   radii(radii),
        curvatures(Vec3(1/radii[0],1/radii[1],1/radii[2])) 
    {   function.setOwner(*this);
        createOBBTree(); }

    ContactGeometryImpl* clone() const {return new Impl(radii);}
    const Vec3& getRadii() const {return radii;}
    void setRadii(const Vec3& r) 
    {   radii = r; curvatures = Vec3(1/r[0],1/r[1],1/r[2]); }

    const Vec3& getCurvatures() const {return curvatures;}

    // See below.
    inline Vec3 findPointWithThisUnitNormal(const UnitVec3& n) const;
    inline Vec3 findPointInSameDirection(const Vec3& Q) const;
    inline UnitVec3 findUnitNormalAtPoint(const Vec3& Q) const;

    // Cost is findParaboloidAtPointWithNormal + about 40 flops.
    void findParaboloidAtPoint(const Vec3& Q, Transform& X_EP, Vec2& k) const
    {   findParaboloidAtPointWithNormal(Q,findUnitNormalAtPoint(Q),X_EP,k); }

    void findParaboloidAtPointWithNormal(const Vec3& Q, const UnitVec3& n,
        Transform& X_EP, Vec2& k) const;

    ContactGeometryTypeId getTypeId() const {return classTypeId();}


    Vec3 findNearestPoint(const Vec3& position, bool& inside, 
                          UnitVec3& normal) const;
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                       Real& distance, UnitVec3& normal) const;
    void getBoundingSphere(Vec3& center, Real& radius) const;

    bool isSmooth() const {return true;}
    bool isConvex() const {return true;}
    bool isFinite() const {return true;}

    // The point furthest in this direction is the unique point whose outward
    // normal is this direction.
    Vec3 calcSupportPoint(UnitVec3 direction) const {
        return findPointWithThisUnitNormal(direction);
    }
    void calcCurvature(const Vec3& point, Vec2& curvature, 
                       Rotation& orientation) const;
    const Function& getImplicitFunction() const {
        return function;
    }

    static ContactGeometryTypeId classTypeId() {
        static const ContactGeometryTypeId id = 
            createNewContactGeometryTypeId();
        return id;
    }
private:
    void createOBBTree();


    Vec3 radii;
    // The curvatures are calculated whenever the radii are set.
    Vec3 curvatures; // (1/radii[0], 1/radii[1], 1/radii[2])
    EllipsoidImplicitFunction function;
};

// Given an ellipsoid and a unit normal direction, find the unique point on the
// ellipsoid whose outward normal matches. The unnormalized normal at a point
// p=[x y z] is the gradient of the ellipsoid's implicit equation there:
// (1)      n(p) = grad(f(p)) = 2*[x/a^2, y/b^2, z/c^2]
// If we had that, we'd have 
// (2)      p=[n[0]*a^2, n[1]*b^2, n[2]*c^2]/2,
// but instead we're given the normalized normal that has been divided
// by the length of n:
// (3)      nn = n/|n| = s * [x/a^2, y/b^2, z/c^2]
// where s = 2/|n|. 
// 
// We can solve for s using the fact that x,y,z must lie on the ellipsoid so
// |x/a,y/b,z/c|=1. Construct the vector
//          v = [nn[0]*a, nn[1]*b, nn[2]*c] = s*[x/a, y/b, z/c]
// Now we have |v|=s. So n/2 = nn/|v| and we can use equation (2) to solve
// for p. Cost is about 40 flops.
inline Vec3 ContactGeometry::Ellipsoid::Impl::
findPointWithThisUnitNormal(const UnitVec3& nn) const {
    const Real& a=radii[0]; const Real& b=radii[1]; const Real& c=radii[2];
    const Vec3 v  = Vec3(nn[0]*a, nn[1]*b, nn[2]*c);
    const Vec3 p  = Vec3( v[0]*a,  v[1]*b,  v[2]*c) / v.norm();
    return p;
}

// Given a point Q=(x,y,z) measured from ellipse center O, find the intersection 
// of the ray d=Q-O with the ellipse surface. This just requires scaling the 
// direction vector d by a factor s so that f(s*d)=0, that is, 
//       s*|x/a y/b z/c|=1  => s = 1/|x/a y/b z/c|
// Cost is about 40 flops.
inline Vec3 ContactGeometry::Ellipsoid::Impl::
findPointInSameDirection(const Vec3& Q) const {
    Real s = 1/Vec3(Q[0]*curvatures[0], 
                    Q[1]*curvatures[1], 
                    Q[2]*curvatures[2]).norm();
    return s*Q;
}

// The implicit equation of the ellipsoid surface is f(x,y,z)=0 where
// f(x,y,z) = (ka x)^2 + (kb y)^2 (kc z)^2 - 1. Points p inside the ellipsoid
// have f(p)<0, outside f(p)>0. f defines a field in space; its positive
// gradient [df/dx df/dy df/dz] points outward. So, given an ellipsoid with 
// principal curvatures ka,kb,kc and a point Q allegedly on the ellipsoid, the
// outward normal (unnormalized) n at that point is
//    n(p) = grad(f(p)) = 2*[ka^2 x, kb^2 y, kc^2 z]
// so the unit norm we're interested in is nn=n/|n| (the "2" drops out).
// If Q is not on the ellipsoid this is equivalent to scaling the ray Q-O
// until it hits the ellipsoid surface at Q'=s*Q, and then reporting the outward
// normal at Q' instead.
// Cost is about 40 flops.
inline UnitVec3 ContactGeometry::Ellipsoid::Impl::
findUnitNormalAtPoint(const Vec3& Q) const {
    const Vec3 kk(square(curvatures[0]), square(curvatures[1]), 
                  square(curvatures[2]));
    return UnitVec3(kk[0]*Q[0], kk[1]*Q[1], kk[2]*Q[2]);
}



//==============================================================================
//                          SMOOTH HEIGHT MAP IMPL
//==============================================================================
class SmoothHeightMapImplicitFunction : public Function {
public:
    SmoothHeightMapImplicitFunction() : ownerp(0) {}
    SmoothHeightMapImplicitFunction
       (const ContactGeometry::SmoothHeightMap::Impl& owner) 
    :   ownerp(&owner) {}
    void setOwner(const ContactGeometry::SmoothHeightMap::Impl& owner) 
    {   ownerp=&owner; }
    Real calcValue(const Vector& x) const;
    Real calcDerivative(const Array_<int>& derivComponents, 
                        const Vector& x) const;
    int getArgumentSize() const {return 3;}
    int getMaxDerivativeOrder() const
    {   return std::numeric_limits<int>::max(); }
private:
    // just a reference; don't delete
    const ContactGeometry::SmoothHeightMap::Impl*   ownerp; 
};



class ContactGeometry::SmoothHeightMap::Impl : public ContactGeometryImpl {
public:
    explicit Impl(const BicubicSurface& surface);

    ContactGeometryImpl* clone() const {
        return new Impl(surface);
    }

    const BicubicSurface& getBicubicSurface() const {return surface;}
    BicubicSurface::PatchHint& updHint() const {return hint;}

    ContactGeometryTypeId getTypeId() const {return classTypeId();}

    Vec3 findNearestPoint(const Vec3& position, bool& inside, 
                          UnitVec3& normal) const;

    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                       Real& distance, UnitVec3& normal) const;

    void getBoundingSphere(Vec3& center, Real& radius) const {
        center = boundingSphere.getCenter();
        radius = boundingSphere.getRadius();
    }

    bool isSmooth() const {return true;}
    bool isConvex() const {return false;}
    bool isFinite() const {return true;}

    Vec3 calcSupportPoint(UnitVec3 direction) const {
        assert(false);
        return Vec3(NaN);
    }

    // We ignore the z coordinate here and just return the curvature of
    // the unique point at (x,y).
    void calcCurvature(const Vec3& point, Vec2& curvature, 
                       Rotation& orientation) const {
        Transform X_SP;
        surface.calcParaboloid(Vec2(point[0],point[1]), hint, X_SP, curvature);
        orientation = X_SP.R();
    }

    const Function& getImplicitFunction() const {return implicitFunction;}

    static ContactGeometryTypeId classTypeId() {
        static const ContactGeometryTypeId id = 
            createNewContactGeometryTypeId();
        return id;
    }
private:
    void createBoundingVolumes();
    // The given OBBNode is assigned this range of patches. If there is
    // more than one patch in the range, it will dole those out to its
    // children recursively until the leaves each have responsibility for
    // one patch. Then we'll deal with that patch, which may have to be
    // subdivided into submission.
    void splitPatches(int x0,int y0, int nx, int ny, 
                      OBBNode& node, int depth,
                      Array_<const Vec3*>* parentControlPoints=0) const;

    // The supplied OBBNode has responsibility for the given subpatch, which
    // may need further subdivision.
    void assignPatch(const Geo::BicubicBezierPatch& patch,
                     OBBNode& node, int depth, 
                     Array_<const Vec3*>* parentControlPoints=0) const;


    BicubicSurface                      surface;
    mutable BicubicSurface::PatchHint   hint;
    Geo::Sphere                         boundingSphere;
    SmoothHeightMapImplicitFunction     implicitFunction;
};




//==============================================================================
//                            OBB TREE NODE IMPL
//==============================================================================
class OBBTreeNodeImpl {
public:
    OBBTreeNodeImpl() : child1(NULL), child2(NULL) {
    }
    OBBTreeNodeImpl(const OBBTreeNodeImpl& copy);
    ~OBBTreeNodeImpl();
    OrientedBoundingBox bounds;
    OBBTreeNodeImpl* child1;
    OBBTreeNodeImpl* child2;
    Array_<int> triangles;
    int numTriangles;
    Vec3 findNearestPoint(const ContactGeometry::TriangleMesh::Impl& mesh, 
                          const Vec3& position, Real cutoff2, Real& distance2, 
                          int& face, Vec2& uv) const;
    bool intersectsRay(const ContactGeometry::TriangleMesh::Impl& mesh, 
                       const Vec3& origin, const UnitVec3& direction, 
                       Real& distance, int& face, Vec2& uv) const;
};



//==============================================================================
//                            TRIANGLE MESH IMPL
//==============================================================================
class ContactGeometry::TriangleMesh::Impl : public ContactGeometryImpl {
public:
    class Edge;
    class Face;
    class Vertex;

    Impl(const ArrayViewConst_<Vec3>& vertexPositions, 
         const ArrayViewConst_<int>& faceIndices, bool smooth);
    Impl(const PolygonalMesh& mesh, bool smooth);
    ContactGeometryImpl* clone() const {
        return new Impl(*this);
    }

    ContactGeometryTypeId getTypeId() const {return classTypeId();}

    Vec3     findPoint(int face, const Vec2& uv) const;
    Vec3     findCentroid(int face) const;
    UnitVec3 findNormalAtPoint(int face, const Vec2& uv) const;
    Vec3 findNearestPoint(const Vec3& position, bool& inside, 
                          UnitVec3& normal) const;
    Vec3 findNearestPoint(const Vec3& position, bool& inside, int& face, 
                          Vec2& uv) const;
    Vec3 findNearestPointToFace(const Vec3& position, int face, Vec2& uv) const;
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                       Real& distance, UnitVec3& normal) const;
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                       Real& distance, int& face, Vec2& uv) const;
    void getBoundingSphere(Vec3& center, Real& radius) const;

    bool isSmooth() const {return false;}
    bool isConvex() const {return false;}
    bool isFinite() const {return true;}

    void createPolygonalMesh(PolygonalMesh& mesh) const;

    static ContactGeometryTypeId classTypeId() {
        static const ContactGeometryTypeId id = 
            createNewContactGeometryTypeId();
        return id;
    }
private:
    void init(const Array_<Vec3>& vertexPositions, const Array_<int>& faceIndices);
    void createObbTree(OBBTreeNodeImpl& node, const Array_<int>& faceIndices);
    void splitObbAxis(const Array_<int>& parentIndices, 
                      Array_<int>& child1Indices, 
                      Array_<int>& child2Indices, int axis);
    void findBoundingSphere(Vec3* point[], int p, int b, 
                            Vec3& center, Real& radius);
    friend class ContactGeometry::TriangleMesh;
    friend class OBBTreeNodeImpl;

    Array_<Edge>    edges;
    Array_<Face>    faces;
    Array_<Vertex>  vertices;
    Vec3            boundingSphereCenter;
    Real            boundingSphereRadius;
    OBBTreeNodeImpl obb;
    bool            smooth;
};



//==============================================================================
//                          TriangleMeshImpl EDGE
//==============================================================================
class ContactGeometry::TriangleMesh::Impl::Edge {
public:
    Edge(int vert1, int vert2, int face1, int face2) {
        vertices[0] = vert1;
        vertices[1] = vert2;
        faces[0] = face1;
        faces[1] = face2;
    }
    int     vertices[2];
    int     faces[2];
};



//==============================================================================
//                           TriangleMeshImpl FACE
//==============================================================================
class ContactGeometry::TriangleMesh::Impl::Face {
public:
    Face(int vert1, int vert2, int vert3, 
         const Vec3& normal, Real area) 
    :   normal(normal), area(area) {
        vertices[0] = vert1;
        vertices[1] = vert2;
        vertices[2] = vert3;
    }
    int         vertices[3];
    int         edges[3];
    UnitVec3    normal;
    Real        area;
};



//==============================================================================
//                          TriangleMeshImpl VERTEX
//==============================================================================
class ContactGeometry::TriangleMesh::Impl::Vertex {
public:
    Vertex(Vec3 pos) : pos(pos), firstEdge(-1) {
    }
    Vec3        pos;
    UnitVec3    normal;
    int         firstEdge;
};

} // namespace SimTK

#endif // SimTK_SIMMATH_CONTACT_GEOMETRY_IMPL_H_
