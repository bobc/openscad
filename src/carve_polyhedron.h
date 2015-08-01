#pragma once

#include "Geometry.h"

#include "memory.h"
#include <string>
#include "linalg.h"

#include "CSGIF.h"

class Carve_volume
{
public:
    Carve_volume (carve::mesh::MeshSet<3> *poly)
    {
        this->poly = poly;
    }

    Carve_volume (const Carve_volume &src)
    {
        this->poly = src.poly->clone();
        this->hasColor = src.hasColor;
        this->color[0] = src.color[0];
        this->color[1] = src.color[1];
        this->color[2] = src.color[2];
        this->color[3] = src.color[3];
    }

    Carve_volume () {
        hasColor = false;
        color[0] = color[1] = color[2] = 0.5f;
    }

    ~Carve_volume () {
    }

    // attributes
    // color
    bool hasColor;
    Color4f color;

    carve::mesh::MeshSet<3> *poly;
};

class CSGIF_polyhedron : public Geometry
{
public:
	CSGIF_polyhedron(Carve_volume *p = NULL);

	CSGIF_polyhedron(const CSGIF_polyhedron &src);
	~CSGIF_polyhedron() {}

	virtual size_t memsize() const;
	// FIXME: Implement, but we probably want a high-resolution BBox..
	virtual BoundingBox getBoundingBox() const { assert(false && "not implemented"); }
	virtual std::string dump() const;
	virtual unsigned int getDimension() const { return 3; }
    // Empty means it is a geometric node which has zero area/volume
	virtual bool isEmpty() const;
	virtual Geometry *copy() const { return new CSGIF_polyhedron(*this); }

	virtual void setColor(const Color4f &c) {
	    color = c; hasColor = true;
	    for (size_t i=0; i < this->Volumes.size(); i++) {
            Volumes[i].hasColor = true;
            Volumes[i].color = c;
	    }
    };

	void reset() {
		Volumes.clear();
	}

	CSGIF_polyhedron &operator+=(const CSGIF_polyhedron &other);
	CSGIF_polyhedron &operator*=(const CSGIF_polyhedron &other);
	CSGIF_polyhedron &operator-=(const CSGIF_polyhedron &other);
	CSGIF_polyhedron &minkowski(const CSGIF_polyhedron &other);

	void transform( const Transform3d &matrix );
	void resize(Vector3d newsize, const Eigen::Matrix<bool,3,1> &autosize);

	//shared_ptr<CSGIF_poly3> poly;   //vector poly_group
	std::vector<Carve_volume> Volumes;

};
