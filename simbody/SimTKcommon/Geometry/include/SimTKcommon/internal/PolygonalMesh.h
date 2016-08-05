#ifndef SimTK_SimTKCOMMON_POLYGONAL_MESH_H_
#define SimTK_SimTKCOMMON_POLYGONAL_MESH_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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

#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/PrivateImplementation.h"

namespace SimTK {

class PolygonalMesh;
class PolygonalMeshImpl;

// We only want the template instantiation to occur once. This symbol is defined
// in the SimTKcommon compilation unit that defines the PolygonalMesh class but 
// should not be defined any other time.
#ifndef SimTK_SIMTKCOMMON_DEFINING_POLYGONALMESH
    extern template class PIMPLHandle<PolygonalMesh, PolygonalMeshImpl, true>;
#endif

/** This class provides a description of a mesh made of polygonal faces. Its 
primary purpose is for loading geometry from files, which can then be used for
visualization or collision detection. For example, the following lines load 
a mesh from a Wavefront OBJ file, then create a DecorativeMesh from it.
@code
    PolygonalMesh mesh;
    std::ifstream file;
    file.open("teapot.obj");
    mesh.loadObjFile(file);
    file.close();
    DecorativeMesh decoration(mesh);
@endcode 
You can also read a polygon mesh from a VTK PolyData (.vtp) file.

We expect this to be a large object so give it shared (reference) semantics; that is,
the copy constructor and copy assignment default to shallow copies (both handles
will refer to the same data). If you want to make a deep (non-shared) copy of a
PolygonalMesh, use the copyAssign() method provided by the PIMPLHandle base class.
**/
class SimTK_SimTKCOMMON_EXPORT PolygonalMesh 
:   public PIMPLHandle<PolygonalMesh, PolygonalMeshImpl, true> {
public:
    /**
     * Create a PolygonalMesh, which initially contains no vertices or faces.
     */
    PolygonalMesh() {}

    /**
     * Get the number of faces in the mesh.
     */
    int getNumFaces() const;
    /**
     * Get the number of vertices in the mesh.
     */
    int getNumVertices() const;
    /**
     * Get the position of a vertex in the mesh.
     *
     * @param vertex  the index of the vertex to get
     * @return the position of the specified vertex
     */
    const Vec3& getVertexPosition(int vertex) const;
    /**
     * Get the number of vertices that make up a face.
     *
     * @param face    the index of the face
     */
    int getNumVerticesForFace(int face) const;
    /**
     * Get the index of one of the vertices of a face.
     *
     * @param face    the index of the face
     * @param vertex  the index of the vertex within the face 
     *                (from 0, 1, or 2 for a triangular face, etc.)
     * @return the index of the specified vertex
     */
    int getFaceVertex(int face, int vertex) const;
    /**
     * Add a vertex to the mesh.
     *
     * @param position    the position of the vertex to add
     * @return the index of the newly added vertex
     */
    int addVertex(const Vec3& position);
    /**
     * Add a face to the mesh.
     *
     * @param vertices    indices of the vertices which make up the new face
     * @return the index of the newly added face
     */
    int addFace(const Array_<int>& vertices);
    /**
     * Scale a mesh by multiplying every vertex by a fixed value.
     */
    void scaleMesh(Real scale);
    /**
     * Transform a mesh by applying a Transform to every vertex.
     */
    void transformMesh(const Transform& transform);
    /**
     * Load a Wavefront OBJ file, adding the vertices and faces it contains
     * to this mesh.
     *
     * @param file    an input stream from which to load the file contents
     */
    void loadObjFile(std::istream& file);
    /**
     * Load a VTK PolyData (.vtp) file, adding the vertices and faces it 
     * contains to this mesh.
     *
     * @param pathname    the name of a .vtp file
     */
    void loadVtpFile(const String& pathname);

private:
    explicit PolygonalMesh(PolygonalMeshImpl* impl) : HandleBase(impl) {}
    void initializeHandleIfEmpty();
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_POLYGONAL_MESH_H_
