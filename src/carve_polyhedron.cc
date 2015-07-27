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
	if (src.poly)
		this->poly.reset(new CSGIF_poly3(*src.poly));
}

CGAL_Nef_polyhedron& CGAL_Nef_polyhedron::operator+=(const CGAL_Nef_polyhedron &other)
{
    PRINT ("csgif union");
//	(*this->poly) += (*other.poly);

    carve::mesh::MeshSet<3> * mesh1 = this->poly->poly;
    carve::mesh::MeshSet<3> * mesh2 = other.poly->poly;

    carve::mesh::MeshSet<3> *result = NULL;

    result = carve::csg::CSG().compute(mesh1, mesh2, carve::csg::CSG::UNION, NULL, carve::csg::CSG::CLASSIFY_NORMAL);

    this->poly->poly = result;
    //this->poly.reset (new CSGIF_poly3 (result));

	return *this;
}

CGAL_Nef_polyhedron& CGAL_Nef_polyhedron::operator*=(const CGAL_Nef_polyhedron &other)
{
//	(*this->poly) *= (*other.poly);

    PRINT ("csgif intersection");

    carve::mesh::MeshSet<3> * mesh1 = this->poly->poly;
    carve::mesh::MeshSet<3> * mesh2 = other.poly->poly;

    carve::mesh::MeshSet<3> *result = NULL;

    result = carve::csg::CSG().compute(mesh1, mesh2, carve::csg::CSG::INTERSECTION, NULL, carve::csg::CSG::CLASSIFY_NORMAL);

    this->poly->poly = result;

	return *this;
}

CGAL_Nef_polyhedron& CGAL_Nef_polyhedron::operator-=(const CGAL_Nef_polyhedron &other)
{
    PRINT ("csgif difference");

    carve::mesh::MeshSet<3> * mesh1 = this->poly->poly;
    carve::mesh::MeshSet<3> * mesh2 = other.poly->poly;

    carve::mesh::MeshSet<3> *result = NULL;

    result = carve::csg::CSG().compute(mesh1, mesh2, carve::csg::CSG::A_MINUS_B, NULL, carve::csg::CSG::CLASSIFY_NORMAL);

    this->poly->poly = result;

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

void csgif_dump_meshset (carve::mesh::MeshSet<3> *meshset)
{
    std::cout << "-- Meshset --" <<  "\n";
    for (size_t mesh_idx = 0; mesh_idx < meshset->meshes.size(); ++mesh_idx)
    {
        std::cout << "Mesh " << mesh_idx << "\n";
        carve::mesh::Mesh<3> *mesh = meshset->meshes[mesh_idx];

        for (size_t face_idx = 0, num_faces = mesh->faces.size(); face_idx != num_faces; ++face_idx)
        {
            carve::mesh::Face<3> *face = mesh->faces[face_idx];
            std::cout << "  Face " << face_idx << "\n";
            carve::mesh::Edge<3> *e = face->edge;
            for (size_t vert_idx = 0, num_verts = face->nVertices(); vert_idx != num_verts; ++vert_idx)
            {
                std::cout << "    v" << vert_idx << " " << e->v1()->v[0] << " " << e->v1()->v[1] << " " << e->v1()->v[2] << " " << "\n";
                e = e->next;
            }
        }
    }
}

std::string CGAL_Nef_polyhedron::dump() const
{
//	return OpenSCAD::dump_svg( *this->poly );
    csgif_dump_meshset(poly->poly);

	return "";
}


void CGAL_Nef_polyhedron::transform( const Transform3d &matrix )
{
	if (!this->isEmpty()) {
		if (matrix.matrix().determinant() == 0) {
			PRINT("WARNING: Scaling a 3D object with 0 - removing object");
			this->reset();
		}
		else {
            PRINT ("transform");
            carve::math::Matrix transform_matrix (
                matrix(0,0), matrix(0,1), matrix(0,2), matrix(0,3),
                matrix(1,0), matrix(1,1), matrix(1,2), matrix(1,3),
                matrix(2,0), matrix(2,1), matrix(2,2), matrix(2,3),
                matrix(3,0), matrix(3,1), matrix(3,2), matrix(3,3)
                );
//            transform_matrix.m [0][0] = matrix(0,0);
//            transform_matrix.m [0][1] = matrix(0,1);
//            transform_matrix.m [0][2] = matrix(0,2);
//            transform_matrix.m [0][3] = matrix(0,3);
//
//            transform_matrix.m [1][0] = matrix(1,0);
//            transform_matrix.m [1][1] = matrix(1,1);
//            transform_matrix.m [1][2] = matrix(1,2);
//            transform_matrix.m [1][3] = matrix(1,3);
//
//            transform_matrix.m [2][0] = matrix(2,0);
//            transform_matrix.m [2][1] = matrix(2,1);
//            transform_matrix.m [2][2] = matrix(2,2);
//            transform_matrix.m [2][3] = matrix(2,3);
//
//            transform_matrix.m [3][3] = matrix(3,3);

            this->poly->poly->transform (carve::math::matrix_transformation(transform_matrix));

//			CGAL_Aff_transformation t(
//				matrix(0,0), matrix(0,1), matrix(0,2), matrix(0,3),
//				matrix(1,0), matrix(1,1), matrix(1,2), matrix(1,3),
//				matrix(2,0), matrix(2,1), matrix(2,2), matrix(2,3), matrix(3,3));
//			this->p3->transform(t);

		}
	}
}
