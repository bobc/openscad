#pragma once

#include "polyset.h"
#include "enums.h"
#include "Geometry.h"

namespace CSGIF_Utils{


	CSGIF_polyhedron *createCsgPolyhedronFromGeometry(const class Geometry &geom);

//	bool createPolySetFromNefPolyhedron3(const CGAL_Nef_polyhedron3 &N, PolySet &ps);
	bool createPolySetFromCsgPolyhedron (const CSGIF_polyhedron &N, PolySet &ps);

    CSGIF_polyhedron *applyOperator(const Geometry::ChildList &children, OpenSCADOperator op);

   	Polygon2d *project(const CSGIF_polyhedron &N, bool cut);

}


