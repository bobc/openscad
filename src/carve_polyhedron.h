#pragma once

#include "Geometry.h"

#include "memory.h"
#include <string>
#include "linalg.h"

#include "CSGIF.h"

class CSGIF_poly3
{
public:
    CSGIF_poly3 (carve::mesh::MeshSet<3> *poly)
    { this->poly = poly; }

    CSGIF_poly3 (const CSGIF_poly3 &src)
    {
        this->poly = src.poly->clone();
    }

    ~CSGIF_poly3 () {}

    carve::mesh::MeshSet<3> *poly;

    // attributes
    // color
    bool has_color;
    Color4f color;

    //std::vector<carve::mesh::MeshSet<3> > *polys;
};

class CSGIF_polyhedron : public Geometry
{
public:
	CSGIF_polyhedron(CSGIF_poly3 *p = NULL);

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

	void reset() {
		poly.reset();
	}
	CSGIF_polyhedron &operator+=(const CSGIF_polyhedron &other);
	CSGIF_polyhedron &operator*=(const CSGIF_polyhedron &other);
	CSGIF_polyhedron &operator-=(const CSGIF_polyhedron &other);
	CSGIF_polyhedron &minkowski(const CSGIF_polyhedron &other);

	void transform( const Transform3d &matrix );
	void resize(Vector3d newsize, const Eigen::Matrix<bool,3,1> &autosize);

	shared_ptr<CSGIF_poly3> poly;   //vector poly_group
};
