#include "GeometryEvaluator.h"
#include "traverser.h"
#include "Tree.h"
#include "Geometry.h"
#include "GeometryCache.h"
#include "CGALCache.h"
#include "Polygon2d.h"
#include "clipper-utils.h"
#include "module.h"
#include "polyset.h"

//#include "cgal.h"
//#include <CGAL/Cartesian.h>
//#include <CGAL/convex_hull_2.h>
//#include <CGAL/Point_2.h>

#include <boost/foreach.hpp>

//#include "carve/carve.hpp"

#include <carve/csg.hpp>
#include <carve/poly.hpp>
#include <carve/geom.hpp>

#include "csgif_polyhedron.h"
#include "csgif_utils.h"

static carve::mesh::MeshSet<3> *carve_test(void)
{
 //create a tetrahedron
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

  carve::mesh::MeshSet<3> tetrahedron(tet_faces);

  //create a triangle
  std::vector<carve::mesh::MeshSet<3>::vertex_t> tri_verts;
  std::vector<carve::mesh::MeshSet<3>::face_t *> tri_faces;

  //Vertices
  //crashes if last coordinate set to 1e-8, but ok for 1e-7
  tri_verts.push_back(carve::mesh::MeshSet<3>::vertex_t(carve::geom::VECTOR(-0.3, 0.0, 1e-8)));
  tri_verts.push_back(carve::mesh::MeshSet<3>::vertex_t(carve::geom::VECTOR(1.0, 0.0, 1.1e-8)));
  tri_verts.push_back(carve::mesh::MeshSet<3>::vertex_t(carve::geom::VECTOR(-0.3, 1.0, 1.1e-8)));

  //Face
  corners.clear();
  corners.push_back(&tri_verts[0]);
  corners.push_back(&tri_verts[2]);
  corners.push_back(&tri_verts[1]);
  tri_faces.push_back(new carve::mesh::MeshSet<3>::face_t(corners.begin(), corners.end()));

//  corners.clear();
//  corners.push_back(&tri_verts[0]);
//  corners.push_back(&tri_verts[1]);
//  corners.push_back(&tri_verts[2]);
//  tri_faces.push_back(new carve::mesh::MeshSet<3>::face_t(corners));

  carve::mesh::MeshSet<3> triangle(tri_faces);

  //cut triangle with tetrahedron.
  carve::mesh::MeshSet<3> *is_poly = carve::csg::CSG().compute(&tetrahedron,
                                                               &triangle,
                                                               carve::csg::CSG::INTERSECTION);

  // std::cout << "Tetrahedron is ... \n" << tetrahedron;
  // std::cout << "Triangle is ... \n" << triangle;
  // std::cout << "Intersection is ... \n" << *is_poly;

    return is_poly;
}


GeometryEvaluator::GeometryEvaluator(const class Tree &tree):
	tree(tree)
{
}


/*!
	Set allownef to false to force the result to _not_ be a Nef polyhedron
*/
shared_ptr<const Geometry> GeometryEvaluator::evaluateGeometry(const AbstractNode &node,
																															 bool allownef)
{
#if 0
        carve::mesh::MeshSet<3> *p = carve_test();
        CSGIF_poly3 *p3 = new CSGIF_poly3 (p);
        CGAL_Nef_polyhedron *pp = new CGAL_Nef_polyhedron(p3);

        shared_ptr<CGAL_Nef_polyhedron> N;

		// If not found in any caches, we need to evaluate the geometry
		if (N) {
			this->root = N;
		}

		//if (shared_ptr<const CGAL_Nef_polyhedron> N = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(this->root))
        this->root.reset(pp);
        return this->root;
#endif

	if (!GeometryCache::instance()->contains(this->tree.getIdString(node))) {
		shared_ptr<const CGAL_Nef_polyhedron> N;
		if (CGALCache::instance()->contains(this->tree.getIdString(node))) {
			N = CGALCache::instance()->get(this->tree.getIdString(node));
		}

		// If not found in any caches, we need to evaluate the geometry
		if (N) {
			this->root = N;
		}
    else {
			Traverser trav(*this, node, Traverser::PRE_AND_POSTFIX);
			trav.execute();
		}

		if (!allownef) {
			if (shared_ptr<const CGAL_Nef_polyhedron> N = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(this->root)) {
				PolySet *ps = new PolySet(3);
				ps->setConvexity(N->getConvexity());
				this->root.reset(ps);
                if (!N->isEmpty()) {
                    bool err = csgif_utils::createPolySetFromCsgPolyhedron(*N, *ps);
                    if (err) {
                        PRINT("ERROR: Nef->PolySet failed");
                    }
                }

				smartCacheInsert(node, this->root);
			}
		}
        return this->root;
	}
	return GeometryCache::instance()->get(this->tree.getIdString(node));
}

GeometryEvaluator::ResultObject GeometryEvaluator::applyToChildren(const AbstractNode &node, OpenSCADOperator op)
{
	unsigned int dim = 0;

	PRINT("operator");

	BOOST_FOREACH(const Geometry::ChildItem &item, this->visitedchildren[node.index()]) {
		if (!item.first->modinst->isBackground() && item.second) {
			if (!dim) dim = item.second->getDimension();
			else if (dim != item.second->getDimension()) {
				PRINT("WARNING: Mixing 2D and 3D objects is not supported.");
				break;
			}
		}
	}
    if (dim == 2) {
        Polygon2d *p2d = applyToChildren2D(node, op);
        assert(p2d);
        return ResultObject(p2d);
    }
    else if (dim == 3) return applyToChildren3D(node, op);
	return ResultObject();
}

/*!
	Applies the operator to all child nodes of the given node.

	May return NULL or any 3D Geometry object (can be either PolySet or CGAL_Nef_polyhedron)
*/
GeometryEvaluator::ResultObject GeometryEvaluator::applyToChildren3D(const AbstractNode &node, OpenSCADOperator op)
{
	Geometry::ChildList children = collectChildren3D(node);
	if (children.size() == 0) return ResultObject();

#if ENABLE_CGAL
	if (op == OPENSCAD_HULL) {
		PolySet *ps = new PolySet(3, true);

		if (CGALUtils::applyHull(children, *ps)) {
			return ps;
		}

		delete ps;
		return ResultObject();
	}
#endif
	// Only one child -> this is a noop
	if (children.size() == 1) return ResultObject(children.front().second);

#if ENABLE_CGAL
	if (op == OPENSCAD_MINKOWSKI) {
		Geometry::ChildList actualchildren;
		BOOST_FOREACH(const Geometry::ChildItem &item, children) {
			if (!item.second->isEmpty()) actualchildren.push_back(item);
		}
		if (actualchildren.empty()) return ResultObject();
		if (actualchildren.size() == 1) return ResultObject(actualchildren.front().second);
		return ResultObject(CGALUtils::applyMinkowski(actualchildren));
	}
#endif // ENABLE_CGAL

	CGAL_Nef_polyhedron *N = csgif_utils::applyOperator(children, op);
	// FIXME: Clarify when we can return NULL and what that means
	if (!N) N = new CGAL_Nef_polyhedron;
	return ResultObject(N);
}


/*!
	Apply 2D hull.

	May return an empty geometry but will not return NULL.
*/
Polygon2d *GeometryEvaluator::applyHull2D(const AbstractNode &node)
{
	std::vector<const Polygon2d *> children = collectChildren2D(node);
	Polygon2d *geometry = new Polygon2d();

#if ENABLE_CGAL
	typedef CGAL::Point_2<CGAL::Cartesian<double> > CGALPoint2;
	// Collect point cloud
	std::list<CGALPoint2> points;
	BOOST_FOREACH(const Polygon2d *p, children) {
		BOOST_FOREACH(const Outline2d &o, p->outlines()) {
			BOOST_FOREACH(const Vector2d &v, o.vertices) {
				points.push_back(CGALPoint2(v[0], v[1]));
			}
		}
	}
	if (points.size() > 0) {
		// Apply hull
		std::list<CGALPoint2> result;
		CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(result));

		// Construct Polygon2d
		Outline2d outline;
		BOOST_FOREACH(const CGALPoint2 &p, result) {
			outline.vertices.push_back(Vector2d(p[0], p[1]));
		}
		geometry->addOutline(outline);
	}
#endif // ENABLE_CGAL
	return geometry;
}

#if ENABLE_CGAL
Geometry *GeometryEvaluator::applyHull3D(const AbstractNode &node)
{
	Geometry::ChildList children = collectChildren3D(node);

	PolySet *P = new PolySet(3);
	if (CGALUtils::applyHull(children, *P)) {
		return P;
	}
	delete P;
	return NULL;
}
#endif // ENABLE_CGAL

Polygon2d *GeometryEvaluator::applyMinkowski2D(const AbstractNode &node)
{
	std::vector<const Polygon2d *> children = collectChildren2D(node);
	if (!children.empty()) {
		return ClipperUtils::applyMinkowski(children);
	}
	return NULL;
}

/*!
	Returns a list of Polygon2d children of the given node.
	May return empty Polygon2d object, but not NULL objects
*/
std::vector<const class Polygon2d *> GeometryEvaluator::collectChildren2D(const AbstractNode &node)
{
	std::vector<const Polygon2d *> children;
	BOOST_FOREACH(const Geometry::ChildItem &item, this->visitedchildren[node.index()]) {
		const AbstractNode *chnode = item.first;
		const shared_ptr<const Geometry> &chgeom = item.second;
		// FIXME: Don't use deep access to modinst members
		if (chnode->modinst->isBackground()) continue;

		// NB! We insert into the cache here to ensure that all children of
		// a node is a valid object. If we inserted as we created them, the
		// cache could have been modified before we reach this point due to a large
		// sibling object.
		smartCacheInsert(*chnode, chgeom);

		if (chgeom) {
			if (chgeom->getDimension() == 2) {
				const Polygon2d *polygons = dynamic_cast<const Polygon2d *>(chgeom.get());
				assert(polygons);
				children.push_back(polygons);
			}
			else {
				PRINT("WARNING: Ignoring 3D child object for 2D operation");
			}
		}
	}
	return children;
}



/*!
	Since we can generate both Nef and non-Nef geometry, we need to insert it into
	the appropriate cache.
	This method inserts the geometry into the appropriate cache if it's not already cached.
*/
void GeometryEvaluator::smartCacheInsert(const AbstractNode &node,
																				 const shared_ptr<const Geometry> &geom)
{
	const std::string &key = this->tree.getIdString(node);

	shared_ptr<const CGAL_Nef_polyhedron> N = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(geom);
	if (N) {
		if (!CGALCache::instance()->contains(key)) CGALCache::instance()->insert(key, N);
	}
	else {
		if (!GeometryCache::instance()->contains(key)) {
			if (!GeometryCache::instance()->insert(key, geom)) {
				PRINT("WARNING: GeometryEvaluator: Node didn't fit into cache");
			}
		}
	}
}

bool GeometryEvaluator::isSmartCached(const AbstractNode &node)
{
	const std::string &key = this->tree.getIdString(node);
	return (GeometryCache::instance()->contains(key) ||
					CGALCache::instance()->contains(key));
}

shared_ptr<const Geometry> GeometryEvaluator::smartCacheGet(const AbstractNode &node, bool preferNef)
{
	const std::string &key = this->tree.getIdString(node);
	shared_ptr<const Geometry> geom;
	bool hasgeom = GeometryCache::instance()->contains(key);
	bool hascgal = CGALCache::instance()->contains(key);
	if (hascgal && (preferNef || !hasgeom)) geom = CGALCache::instance()->get(key);
	else if (hasgeom) geom = GeometryCache::instance()->get(key);
	return geom;
}

/*!
	Returns a list of 3D Geometry children of the given node.
	May return empty geometries, but not NULL objects
*/
Geometry::ChildList GeometryEvaluator::collectChildren3D(const AbstractNode &node)
{
	Geometry::ChildList children;
	BOOST_FOREACH(const Geometry::ChildItem &item, this->visitedchildren[node.index()]) {
		const AbstractNode *chnode = item.first;
		const shared_ptr<const Geometry> &chgeom = item.second;
		// FIXME: Don't use deep access to modinst members
		if (chnode->modinst->isBackground()) continue;

		// NB! We insert into the cache here to ensure that all children of
		// a node is a valid object. If we inserted as we created them, the
		// cache could have been modified before we reach this point due to a large
		// sibling object.
		smartCacheInsert(*chnode, chgeom);

		if (chgeom) {
			if (chgeom->getDimension() == 2) {
				PRINT("WARNING: Ignoring 2D child object for 3D operation");
			}
			else if (chgeom->isEmpty() || chgeom->getDimension() == 3) {
				children.push_back(item);
			}
		}
	}
	return children;
}

/*!

*/
Polygon2d *GeometryEvaluator::applyToChildren2D(const AbstractNode &node, OpenSCADOperator op)
{
	if (op == OPENSCAD_MINKOWSKI) {
		return applyMinkowski2D(node);
	}
	else if (op == OPENSCAD_HULL) {
		return applyHull2D(node);
	}

	std::vector<const Polygon2d *> children = collectChildren2D(node);

	if (children.empty()) {
		return NULL;
	}

	if (children.size() == 1) {
		return new Polygon2d(*children[0]); // Copy
	}

	ClipperLib::ClipType clipType;
	switch (op) {
	case OPENSCAD_UNION:
		clipType = ClipperLib::ctUnion;
		break;
	case OPENSCAD_INTERSECTION:
		clipType = ClipperLib::ctIntersection;
		break;
	case OPENSCAD_DIFFERENCE:
		clipType = ClipperLib::ctDifference;
		break;
	default:
		PRINTB("Error: Unknown boolean operation %d", int(op));
		return NULL;
		break;
	}

	return ClipperUtils::apply(children, clipType);
}

/*!
	Adds ourself to out parent's list of traversed children.
	Call this for _every_ node which affects output during traversal.
	Usually, this should be called from the postfix stage, but for some nodes,
	we defer traversal letting other components (e.g. CGAL) render the subgraph,
	and we'll then call this from prefix and prune further traversal.

	The added geometry can be NULL if it wasn't possible to evaluate it.
*/
void GeometryEvaluator::addToParent(const State &state,
																		const AbstractNode &node,
																		const shared_ptr<const Geometry> &geom)
{
	this->visitedchildren.erase(node.index());
	if (state.parent()) {
		this->visitedchildren[state.parent()->index()].push_back(std::make_pair(&node, geom));
	}
	else {
		// Root node, insert into cache
		smartCacheInsert(node, geom);
		this->root = geom;
        assert(this->visitedchildren.empty());
	}
}

/*!
   Custom nodes are handled here => implicit union
*/
Response GeometryEvaluator::visit(State &state, const AbstractNode &node)
{
    PRINT ("implicit");

	if (state.isPrefix()) {
		if (isSmartCached(node)) return PruneTraversal;
		state.setPreferNef(true); // Improve quality of CSG by avoiding conversion loss
	}
	if (state.isPostfix()) {
		shared_ptr<const class Geometry> geom;
		if (!isSmartCached(node)) {
			geom = applyToChildren(node, OPENSCAD_UNION).constptr();
		}
		else {
			geom = smartCacheGet(node, state.preferNef());
		}
		addToParent(state, node, geom);
	}
	return ContinueTraversal;
}

Response GeometryEvaluator::visit(State &state, const OffsetNode &node)
{
	return AbortTraversal;
}

/*!
   RenderNodes just pass on convexity
*/
Response GeometryEvaluator::visit(State &state, const RenderNode &node)
{
	return AbortTraversal;
}

/*!
	Leaf nodes can create their own geometry, so let them do that

	input: None
	output: PolySet or Polygon2d
*/
Response GeometryEvaluator::visit(State &state, const LeafNode &node)
{
	PRINT("leaf");

    if (state.isPrefix()) {
		shared_ptr<const Geometry> geom;
		if (!isSmartCached(node)) {
			const Geometry *geometry = node.createGeometry();
            assert(geometry);
			if (const Polygon2d *polygon = dynamic_cast<const Polygon2d*>(geometry)) {
				if (!polygon->isSanitized()) {
					Polygon2d *p = ClipperUtils::sanitize(*polygon);
					delete geometry;
					geometry = p;
				}
			}
            geom.reset(geometry);
		}
		else geom = smartCacheGet(node, state.preferNef());
		addToParent(state, node, geom);
	}
	return PruneTraversal;
}

Response GeometryEvaluator::visit(State &state, const TextNode &node)
{
	return AbortTraversal;
}


/*!
	input: List of 2D or 3D objects (not mixed)
	output: Polygon2d or 3D PolySet
	operation:
	  o Perform csg op on children
 */
Response GeometryEvaluator::visit(State &state, const CsgNode &node)
{
	return AbortTraversal;
}

/*!
	input: List of 2D or 3D objects (not mixed)
	output: Polygon2d or 3D PolySet
	operation:
	  o Union all children
	  o Perform transform
 */
Response GeometryEvaluator::visit(State &state, const TransformNode &node)
{
	return AbortTraversal;
}


/*!
	input: List of 2D objects
	output: 3D PolySet
	operation:
	  o Union all children
	  o Perform extrude
 */
Response GeometryEvaluator::visit(State &state, const LinearExtrudeNode &node)
{
	return AbortTraversal;
}



/*!
	input: List of 2D objects
	output: 3D PolySet
	operation:
	  o Union all children
	  o Perform extrude
 */
Response GeometryEvaluator::visit(State &state, const RotateExtrudeNode &node)
{
	return AbortTraversal;
}

/*!
	Handles non-leaf PolyNodes; projection
*/
Response GeometryEvaluator::visit(State &state, const AbstractPolyNode &node)
{
	assert(false);
	return AbortTraversal;
}

/*!
	input: List of 3D objects
	output: Polygon2d
	operation:
	  o Union all children
		o Perform projection
 */
Response GeometryEvaluator::visit(State &state, const ProjectionNode &node)
{
	return AbortTraversal;
}

/*!
	input: List of 2D or 3D objects (not mixed)
	output: any Geometry
	operation:
	  o Perform cgal operation
 */
Response GeometryEvaluator::visit(State &state, const CgaladvNode &node)
{
	return AbortTraversal;
}

Response GeometryEvaluator::visit(State &state, const AbstractIntersectionNode &node)
{
	return AbortTraversal;
}


