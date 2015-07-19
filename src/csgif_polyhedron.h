#pragma once

#include "Geometry.h"
//#include "cgal.h"
#include "memory.h"
#include <string>
#include "linalg.h"

#include <carve/csg.hpp>
#include <carve/poly.hpp>
#include <carve/geom.hpp>

class CSGIF_poly3
{
public:

    CSGIF_poly3 (carve::mesh::MeshSet<3> *poly)
    { this->poly = poly; }

    CSGIF_poly3 (const CSGIF_poly3 &src)
    {
        // todo: copy mesh?
        this->poly = src.poly;
    }

    ~CSGIF_poly3 () {}

    carve::mesh::MeshSet<3> *poly;
};

class CGAL_Nef_polyhedron : public Geometry
{
public:
//	CGAL_Nef_polyhedron(CGAL_Nef_polyhedron3 *p = NULL);

	CGAL_Nef_polyhedron(CSGIF_poly3 *p = NULL);

	CGAL_Nef_polyhedron(const CGAL_Nef_polyhedron &src);
	~CGAL_Nef_polyhedron() {}

	virtual size_t memsize() const;
	// FIXME: Implement, but we probably want a high-resolution BBox..
	virtual BoundingBox getBoundingBox() const { assert(false && "not implemented"); }
	virtual std::string dump() const;
	virtual unsigned int getDimension() const { return 3; }
  // Empty means it is a geometric node which has zero area/volume
	virtual bool isEmpty() const;
	virtual Geometry *copy() const { return new CGAL_Nef_polyhedron(*this); }

	void reset() {
		//p3.reset();
	}
	CGAL_Nef_polyhedron &operator+=(const CGAL_Nef_polyhedron &other);
	CGAL_Nef_polyhedron &operator*=(const CGAL_Nef_polyhedron &other);
	CGAL_Nef_polyhedron &operator-=(const CGAL_Nef_polyhedron &other);
	CGAL_Nef_polyhedron &minkowski(const CGAL_Nef_polyhedron &other);
// FIXME: Deprecated by CGALUtils::createPolySetFromNefPolyhedron3
//	class PolySet *convertToPolyset() const;
	void transform( const Transform3d &matrix );
	void resize(Vector3d newsize, const Eigen::Matrix<bool,3,1> &autosize);

	shared_ptr<CSGIF_poly3> poly;

};
