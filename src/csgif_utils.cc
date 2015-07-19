#include "csgif_utils.h"
#include "csgif_polyhedron.h"

#include "printutils.h"
#include "Polygon2d.h"
#include "polyset-utils.h"

#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>

/*!
	Triangulates this polygon2d and returns a 2D PolySet.
*/
PolySet *Polygon2d::tessellate() const
{
	PRINTDB("Polygon2d::tessellate(): %d outlines", this->outlines().size());
	PolySet *polyset = new PolySet(*this);

    // TODO: Stub

	return polyset;
}


static CGAL_Nef_polyhedron *createCsgPolyhedronFromPolySet(const PolySet &ps)
{
    PRINT ("createCsgPolyhedronFromPolySet");

	if (ps.isEmpty()) return new CGAL_Nef_polyhedron();
	assert(ps.getDimension() == 3);

	// Since is_convex doesn't work well with non-planar faces,
	// we tessellate the polyset before checking.
	PolySet psq(ps);
	psq.quantizeVertices();
	PolySet ps_tri(3, psq.convexValue());
	PolysetUtils::tessellate_faces(psq, ps_tri);

	if (ps_tri.is_convex())
    {
        PRINT ("is_convex");

        //std::vector<carve::mesh::MeshSet<3>::vertex_t> verts;
        std::vector<carve::mesh::MeshSet<3>::vertex_t *> corners;
        std::vector<carve::mesh::MeshSet<3>::face_t *> faces;
        //int j = 0;

		BOOST_FOREACH(const Polygon &poly, psq.polygons) {
		    // poly
			corners.clear();
			BOOST_FOREACH(const Vector3d &p, poly) {
				//verts.push_back( carve::mesh::MeshSet<3>::vertex_t(carve::geom::VECTOR(p[0], p[1], p[2])) );
				carve::mesh::MeshSet<3>::vertex_t *v = new carve::mesh::MeshSet<3>::vertex_t(carve::geom::VECTOR(p[0], p[1], p[2]));
                corners.push_back(v);
			}
            faces.push_back(new carve::mesh::MeshSet<3>::face_t(corners.begin(), corners.end()));
		}

        carve::mesh::MeshSet<3> *poly_mesh = new carve::mesh::MeshSet<3>(faces);

        CSGIF_poly3 *p3 = new CSGIF_poly3 (poly_mesh);

        return new CGAL_Nef_polyhedron(p3);
#if 0
		typedef CGAL::Epick K;
		// Collect point cloud
		// FIXME: Use unordered container (need hash)
		// NB! CGAL's convex_hull_3() doesn't like std::set iterators, so we use a list
		// instead.
		std::list<K::Point_3> points;
		BOOST_FOREACH(const Polygon &poly, psq.polygons) {
			BOOST_FOREACH(const Vector3d &p, poly) {
				points.push_back(vector_convert<K::Point_3>(p));
			}
		}

		if (points.size() <= 3) return new CGAL_Nef_polyhedron();;

		// Apply hull
		CGAL::Polyhedron_3<K> r;
		CGAL::convex_hull_3(points.begin(), points.end(), r);
		CGAL_Polyhedron r_exact;
		CGALUtils::copyPolyhedron(r, r_exact);
		return new CGAL_Nef_polyhedron(new CGAL_Nef_polyhedron3(r_exact));
#endif // 0
	}
    else
    {
        PRINT ("not_convex");
    }

	CSGIF_poly3 *N = NULL;
	bool plane_error = false;
#if 0
	//CGAL::Failure_behaviour old_behaviour = CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
	try {
		CGAL_Polyhedron P;
		bool err = CGALUtils::createPolyhedronFromPolySet(psq, P);
		 if (!err) {
		 	PRINTDB("Polyhedron is closed: %d", P.is_closed());
		 	PRINTDB("Polyhedron is valid: %d", P.is_valid(false, 0));
		 }

		if (!err) N = new CGAL_Nef_polyhedron3(P);
	}
	catch (const CGAL::Assertion_exception &e) {
		if (std::string(e.what()).find("Plane_constructor")!=std::string::npos &&
            std::string(e.what()).find("has_on")!=std::string::npos) {
				PRINT("PolySet has nonplanar faces. Attempting alternate construction");
				plane_error=true;
		} else {
			PRINTB("ERROR: CGAL error in CGAL_Nef_polyhedron3(): %s", e.what());
		}
	}
	if (plane_error) try {
			CGAL_Polyhedron P;
			bool err = CGALUtils::createPolyhedronFromPolySet(ps_tri, P);
            if (!err) {
                PRINTDB("Polyhedron is closed: %d", P.is_closed());
                PRINTDB("Polyhedron is valid: %d", P.is_valid(false, 0));
            }
			if (!err) N = new CGAL_Nef_polyhedron3(P);
		}
		catch (const CGAL::Assertion_exception &e) {
			PRINTB("ERROR: Alternate construction failed. CGAL error in CGAL_Nef_polyhedron3(): %s", e.what());
		}
	//CGAL::set_error_behaviour(old_behaviour);
#endif // 0

	return new CGAL_Nef_polyhedron(N);
}

static CGAL_Nef_polyhedron *createCsgPolyhedronFromPolygon2d(const Polygon2d &polygon)
{
	shared_ptr<PolySet> ps(polygon.tessellate());
	return createCsgPolyhedronFromPolySet(*ps);
}

namespace csgif_utils{


    CGAL_Nef_polyhedron *createCsgPolyhedronFromGeometry(const class Geometry &geom)
    {
        PRINT ("createCsgPolyhedronFromGeometry");

   		const PolySet *ps = dynamic_cast<const PolySet*>(&geom);
		if (ps) {
			return createCsgPolyhedronFromPolySet(*ps);
		}
		else {
			const Polygon2d *poly2d = dynamic_cast<const Polygon2d*>(&geom);
			if (poly2d) return createCsgPolyhedronFromPolygon2d(*poly2d);
		}
		assert(false && "createNefPolyhedronFromGeometry(): Unsupported geometry type");
		return NULL;
    }

#if 0
        std::vector<carve::mesh::MeshSet<3>::vertex_t> tet_verts;
        std::vector<carve::mesh::MeshSet<3>::face_t *> tet_faces;
        std::vector<carve::mesh::MeshSet<3>::vertex_t *> corners;

        tet_verts.push_back(carve::mesh::MeshSet<3>::vertex_t(carve::geom::VECTOR(0.0, 0.0, 0.0)));
        tet_verts.push_back(carve::mesh::MeshSet<3>::vertex_t(carve::geom::VECTOR(1.0, 0.0, 0.0)));
        tet_verts.push_back(carve::mesh::MeshSet<3>::vertex_t(carve::geom::VECTOR(0.0, 1.0, 0.0)));
        tet_verts.push_back(carve::mesh::MeshSet<3>::vertex_t(carve::geom::VECTOR(0.0, 0.0, 1.0)));

        corners.push_back(&tet_verts[0]);
        corners.push_back(&tet_verts[2]);
        corners.push_back(&tet_verts[1]);
        tet_faces.push_back(new carve::mesh::MeshSet<3>::face_t(corners.begin(), corners.end()));

        corners.clear();
        corners.push_back(&tet_verts[0]);
        corners.push_back(&tet_verts[1]);
        corners.push_back(&tet_verts[3]);
        tet_faces.push_back(new carve::mesh::MeshSet<3>::face_t(corners.begin(), corners.end()));

        corners.clear();
        corners.push_back(&tet_verts[0]);
        corners.push_back(&tet_verts[3]);
        corners.push_back(&tet_verts[2]);
        tet_faces.push_back(new carve::mesh::MeshSet<3>::face_t(corners.begin(), corners.end()));

        corners.clear();
        corners.push_back(&tet_verts[1]);
        corners.push_back(&tet_verts[2]);
        corners.push_back(&tet_verts[3]);
        tet_faces.push_back(new carve::mesh::MeshSet<3>::face_t(corners.begin(), corners.end()));

        carve::mesh::MeshSet<3> *tetrahedron = new carve::mesh::MeshSet<3>(tet_faces);

        CSGIF_poly3 *p3 = new CSGIF_poly3 (tetrahedron);
        CGAL_Nef_polyhedron *pp = new CGAL_Nef_polyhedron(p3);

        return pp;
    }
#endif

    bool createPolySetFromCsgPolyhedron (const CGAL_Nef_polyhedron &N, PolySet &ps)
    {
        PRINT ("createPolySetFromCsgPolyhedron");
        return false;
    }

/*!
	Applies op to all children and returns the result.
	The child list should be guaranteed to contain non-NULL 3D or empty Geometry objects
*/
	CGAL_Nef_polyhedron *applyOperator(const Geometry::ChildList &children, OpenSCADOperator op)
	{
		CGAL_Nef_polyhedron *N = NULL;
		PRINT ("applyOperator");

		try {

			BOOST_FOREACH(const Geometry::ChildItem &item, children) {
				const shared_ptr<const Geometry> &chgeom = item.second;
				shared_ptr<const CGAL_Nef_polyhedron> chN = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(chgeom);
				if (!chN) {
					const PolySet *chps = dynamic_cast<const PolySet*>(chgeom.get());
					if (chps) chN.reset(createCsgPolyhedronFromGeometry(*chps));
				}

				// Initialize N with first expected geometric object
				if (!N) {
					N = new CGAL_Nef_polyhedron(*chN);
					continue;
				}

				// Intersecting something with nothing results in nothing
				if (chN->isEmpty()) {
					if (op == OPENSCAD_INTERSECTION) *N = *chN;
					continue;
				}

				// empty op <something> => empty
				//! if (N->isEmpty()) continue;

				switch (op) {

                case OPENSCAD_UNION:
                    *N += *chN;
                    break;

				case OPENSCAD_INTERSECTION:
					*N *= *chN;
					break;
				case OPENSCAD_DIFFERENCE:
					*N -= *chN;
					break;
				case OPENSCAD_MINKOWSKI:
					N->minkowski(*chN);
					break;
				default:
					PRINTB("ERROR: Unsupported CGAL operator: %d", op);
				}

				//TODO
				//!item.first->progress_report();
			}
		}
        catch (const carve::exception &ex)
        {
            PRINTB ("ERROR: CSG exception %s", ex.str());
        }

        return N;
    }

}

