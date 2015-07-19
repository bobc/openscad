#pragma once

#include "polyset.h"
#include "csgif_polyhedron.h"
#include "enums.h"
#include "Geometry.h"

namespace csgif_utils{


	CGAL_Nef_polyhedron *createCsgPolyhedronFromGeometry(const class Geometry &geom);

//	bool createPolySetFromNefPolyhedron3(const CGAL_Nef_polyhedron3 &N, PolySet &ps);
	bool createPolySetFromCsgPolyhedron (const CGAL_Nef_polyhedron &N, PolySet &ps);

    CGAL_Nef_polyhedron *applyOperator(const Geometry::ChildList &children, OpenSCADOperator op);
}


