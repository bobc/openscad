
#include "printutils.h"
#include "Polygon2d.h"
#include "polyset-utils.h"

#include "CSGIF.h"
#include <carve/input.hpp>

#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>


extern void csgif_dump_meshset (carve::mesh::MeshSet<3> *meshset);

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


static carve::mesh::MeshSet<3> *makeCube(const carve::math::Matrix &transform)
{
  carve::input::PolyhedronData data;

  data.addVertex(transform * carve::geom::VECTOR(+1.0, +1.0, +1.0));
  data.addVertex(transform * carve::geom::VECTOR(-1.0, +1.0, +1.0));
  data.addVertex(transform * carve::geom::VECTOR(-1.0, -1.0, +1.0));
  data.addVertex(transform * carve::geom::VECTOR(+1.0, -1.0, +1.0));
  data.addVertex(transform * carve::geom::VECTOR(+1.0, +1.0, -1.0));
  data.addVertex(transform * carve::geom::VECTOR(-1.0, +1.0, -1.0));
  data.addVertex(transform * carve::geom::VECTOR(-1.0, -1.0, -1.0));
  data.addVertex(transform * carve::geom::VECTOR(+1.0, -1.0, -1.0));

  data.addFace(0, 1, 2, 3);
  data.addFace(7, 6, 5, 4);
  data.addFace(0, 4, 5, 1);
  data.addFace(1, 5, 6, 2);
  data.addFace(2, 6, 7, 3);
  data.addFace(3, 7, 4, 0);

  return new carve::mesh::MeshSet<3>(data.points, data.getFaceCount(), data.faceIndices);
}

static bool is_equal (double v1, double v2)
{
    return abs(v1-v2) < 0.000000001;
}

#define Abs(x)    ((x) < 0 ? -(x) : (x))
#define Max(a, b) ((a) > (b) ? (a) : (b))

#define EPSILON 0.00000000001

bool nearly_equal(double a, double b)
{
	double c = Abs(a);
	double d = Abs(b);

    if ((a<EPSILON) && (b<EPSILON))
        return true;

	d = Max(c, d);

    if (d != 0.0)
        d = Abs(a - b) / d;

	return d < EPSILON;
}

static CSGIF_polyhedron *createCsgPolyhedronFromPolySet(const PolySet &ps)
{
    PRINT ("createCsgPolyhedronFromPolySet");

	if (ps.isEmpty()) return new CSGIF_polyhedron();
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
#if 0
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
#endif // 0

//        carve::mesh::MeshSet<3> *poly_mesh = makeCube(carve::math::Matrix::SCALE(10.0, 10.0, 10.0));

#if 1
        carve::input::PolyhedronData data;

		std::vector<std::string> vertices;

		BOOST_FOREACH(const Polygon &poly, psq.polygons) {
		    // poly
			std::vector<int> face_verts;
			BOOST_FOREACH(const Vector3d &p, poly) {

                int idx=0;
                // std::cout << "    v" << idx << " " << p[0] << " " << p[1] << " " << p[2] << " " << "\n";

                // This might not be the most obvious method, but appears more reliable
                // than comparing doubles.
				std::stringstream stream;
				stream << p[0] << " " << p[1] << " " << p[2];
				std::string vs1 = stream.str();

                std::vector<std::string>::iterator it;
                it = std::find(vertices.begin(), vertices.end(), vs1);

                if (it != vertices.end())
                {
                    idx = std::distance(vertices.begin(), it);
                    face_verts.push_back (idx);
                }
                else
                {
                    data.addVertex(carve::geom::VECTOR(p[0], p[1], p[2]));
                    vertices.push_back(vs1);
                    face_verts.push_back (vertices.size()-1);
                }

#if 0
                // this produces lots of bad mesh
                for (; idx < data.points.size(); idx++)
                {
                    if ( is_equal(p[0], data.points[idx].x) &&
                         is_equal(p[1], data.points[idx].y) &&
                         is_equal(p[2], data.points[idx].z)
                        )
                    {
                        face_verts.push_back (idx);
                        break;
                    }
                }
                if (idx == data.points.size())
                {
                    face_verts.push_back (data.addVertex(carve::geom::VECTOR(p[0], p[1], p[2])));
                }
#endif
			}

            data.addFace(face_verts.begin(), face_verts.end());
		}

        carve::mesh::MeshSet<3> *poly_mesh = new carve::mesh::MeshSet<3>(data.points, data.getFaceCount(), data.faceIndices);
#endif
//
        if (poly_mesh->meshes.size() > 1)
            PRINT("WARNING: bad poly conversion?");

        Carve_volume *new_volume = new Carve_volume (poly_mesh);
        new_volume->hasColor = ps.hasColor;
        new_volume->color = ps.color;

        CSGIF_polyhedron *result = new CSGIF_polyhedron(new_volume);

        return result;

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

		if (points.size() <= 3) return new CSGIF_polyhedron();;

		// Apply hull
		CGAL::Polyhedron_3<K> r;
		CGAL::convex_hull_3(points.begin(), points.end(), r);
		CGAL_Polyhedron r_exact;
		CGALUtils::copyPolyhedron(r, r_exact);
		return new CSGIF_polyhedron(new CGAL_Nef_polyhedron3(r_exact));
#endif // 0
	}
    else
    {
        PRINT ("not_convex");
    }

	Carve_volume *volume = NULL;
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

	return new CSGIF_polyhedron(volume);
}

static CSGIF_polyhedron *createCsgPolyhedronFromPolygon2d(const Polygon2d &polygon)
{
	shared_ptr<PolySet> ps(polygon.tessellate());
	return createCsgPolyhedronFromPolySet(*ps);
}

namespace CSGIF_Utils{


    CSGIF_polyhedron *createCsgPolyhedronFromGeometry(const class Geometry &geom)
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

    bool createPolySetFromCsgPolyhedron (const CSGIF_polyhedron &N, PolySet &ps)
    {
        PRINT ("createPolySetFromCsgPolyhedron");
        return false;
    }

/*!
	Applies op to all children and returns the result.
	The child list should be guaranteed to contain non-NULL 3D or empty Geometry objects
*/
	CSGIF_polyhedron *applyOperator(const Geometry::ChildList &children, OpenSCADOperator op)
	{
		CSGIF_polyhedron *N = NULL;
		PRINT ("applyOperator");

		try {

			BOOST_FOREACH(const Geometry::ChildItem &item, children) {
				const shared_ptr<const Geometry> &chgeom = item.second;
				shared_ptr<const CSGIF_polyhedron> chN = dynamic_pointer_cast<const CSGIF_polyhedron>(chgeom);
				if (!chN) {
					const PolySet *chps = dynamic_cast<const PolySet*>(chgeom.get());
					if (chps) chN.reset(createCsgPolyhedronFromGeometry(*chps));
				}

				// Initialize N with first expected geometric object
				if (!N) {
					N = new CSGIF_polyhedron(*chN);
					continue;
				}

				// Intersecting something with nothing results in nothing
				if (chN->isEmpty()) {
					if (op == OPENSCAD_INTERSECTION) *N = *chN;
					continue;
				}

				// empty op <something> => empty
				if (N->isEmpty()) continue;

#if 0
                //todo: debug
                std::cout << "-- N (left) --" << "\n";
                N->dump();
                std::cout << "-- chN (right) --" << "\n";
                chN->dump();
#endif // 0
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

                //std::cout << "-- RESULT --" << "\n";
                //N->dump();

				//TODO : not sure what this is, gives error?
				//!item.first->progress_report();
			}
		}
        catch (const carve::exception &ex)
        {
            PRINTB ("ERROR: CSG exception %s", ex.str());
        }

        return N;
    }

    // from cgalutils-project.cc
	Polygon2d *project(const CSGIF_polyhedron &N, bool cut)
	{
		Polygon2d *poly = NULL;
		if (N.getDimension() != 3) return poly;

        // todo
        PRINT("ERROR: not implemented");

        return poly;
	}
}

