#include "IndexedMesh.h"

IndexedMesh::IndexedMesh()
{
    //ctor
}

IndexedMesh::~IndexedMesh()
{
    //dtor
}

IndexedMesh::IndexedMesh(const IndexedMesh& other)
{
    //copy ctor
    this->polys = other.polys;
    this->hasColor = other.hasColor;
    this->color = other.color;
    this->hasMaterialId = other.hasMaterialId;
    this->materialId = other.materialId;
}

IndexedMesh& IndexedMesh::operator=(const IndexedMesh& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator

    this->polys = rhs.polys;
    this->hasColor = rhs.hasColor;
    this->color = rhs.color;
    this->hasMaterialId = rhs.hasMaterialId;
    this->materialId = rhs.materialId;

    return *this;
}


void IndexedMesh::AddMesh (const IndexedMesh& other)
{
    int baseIndex = polys.faces.size();

    polys.vertices.insert( polys.vertices.end(), other.polys.vertices.begin(), other.polys.vertices.end());

    for (size_t i=0; i < other.polys.faces.size(); i ++) {
        IndexedFace face = other.polys.faces[i];

        for (size_t j=0; j < face.size(); j++)
            face[j] += baseIndex;
        polys.faces.push_back(face);
    }

}
