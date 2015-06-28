#include "csgif_utils.h"

#include "printutils.h"

namespace csgif_utils{


    CGAL_Nef_polyhedron *createCsgPolyhedronFromGeometry(const class Geometry &geom)
    {
        PRINT ("createCsgPoly");
        return NULL;
    }

/*!
	Applies op to all children and returns the result.
	The child list should be guaranteed to contain non-NULL 3D or empty Geometry objects
*/
	CGAL_Nef_polyhedron *applyOperator(const Geometry::ChildList &children, OpenSCADOperator op)
	{
		CGAL_Nef_polyhedron *N = NULL;
		PRINT ("applyOperator");
        return N;
    }

}

