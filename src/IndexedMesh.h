
#pragma once

#include "GeometryUtils.h"


struct imIndexedPolygons {
	std::vector<Vector3d> vertices;
	std::vector<IndexedFace> faces;
};

class IndexedMesh
{
    public:
        IndexedMesh();
        virtual ~IndexedMesh();
        IndexedMesh(const IndexedMesh& other);
        IndexedMesh& operator=(const IndexedMesh& other);

        void AddMesh (const IndexedMesh& other);

        imIndexedPolygons polys;

        bool hasColor;
        Color4f color;

        bool hasMaterialId;
        int materialId;

    private:
};
