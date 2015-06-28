//#include "CGAL_Nef_polyhedron.h"
//#include "cgalutils.h"
#include "csgif_polyhedron.h"
#include "printutils.h"
#include "polyset.h"
//#include "svg.h"

CGAL_Nef_polyhedron::CGAL_Nef_polyhedron(CSGIF_poly3 *p)
{
	if (p) poly.reset(p);
}

// Copy constructor
CGAL_Nef_polyhedron::CGAL_Nef_polyhedron(const CGAL_Nef_polyhedron &src)
{
//	if (src.poly)
//		this->poly.reset(new CGAL_Nef_polyhedron3(*src.poly));
}

CGAL_Nef_polyhedron& CGAL_Nef_polyhedron::operator+=(const CGAL_Nef_polyhedron &other)
{
//	(*this->poly) += (*other.poly);
	return *this;
}

CGAL_Nef_polyhedron& CGAL_Nef_polyhedron::operator*=(const CGAL_Nef_polyhedron &other)
{
//	(*this->poly) *= (*other.poly);
	return *this;
}

CGAL_Nef_polyhedron& CGAL_Nef_polyhedron::operator-=(const CGAL_Nef_polyhedron &other)
{
//	(*this->poly) -= (*other.poly);
	return *this;
}

CGAL_Nef_polyhedron &CGAL_Nef_polyhedron::minkowski(const CGAL_Nef_polyhedron &other)
{
	//(*this->poly) = CGAL::minkowski_sum_3(*this->poly, *other.poly);
	return *this;
}

size_t CGAL_Nef_polyhedron::memsize() const
{
	if (this->isEmpty()) return 0;

	size_t memsize = sizeof(CGAL_Nef_polyhedron);
	//memsize += this->poly->bytes();
	return memsize;
}

bool CGAL_Nef_polyhedron::isEmpty() const
{
	//return !this->poly || this->poly->is_empty();
	return !this->poly;
}


void CGAL_Nef_polyhedron::resize(Vector3d newsize,
																 const Eigen::Matrix<bool,3,1> &autosize)
{
}

std::string CGAL_Nef_polyhedron::dump() const
{
//	return OpenSCAD::dump_svg( *this->poly );
	return "";
}


void CGAL_Nef_polyhedron::transform( const Transform3d &matrix )
{

}
