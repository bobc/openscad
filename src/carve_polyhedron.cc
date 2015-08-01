
#include "printutils.h"
#include "polyset.h"

#include "CSGIF.h"

CSGIF_polyhedron::CSGIF_polyhedron(Carve_volume *p)
{
	if (p)
	{
	    Volumes.push_back(*p);
	}
}

// Copy constructor
CSGIF_polyhedron::CSGIF_polyhedron(const CSGIF_polyhedron &src)
{
    //todo
	if (src.Volumes.size() != 0)
		this->Volumes = src.Volumes;

    hasColor = src.hasColor;
    color = src.color;
}

CSGIF_polyhedron& CSGIF_polyhedron::operator+=(const CSGIF_polyhedron &other)
{
    PRINT ("csgif union");
//	(*this->poly) += (*other.poly);
#if 0
    carve::mesh::MeshSet<3> * mesh1 = this->Volumes[0].poly;
    carve::mesh::MeshSet<3> * mesh2 = other.Volumes[0].poly;
    carve::mesh::MeshSet<3> *result = carve::csg::CSG().compute(mesh1, mesh2, carve::csg::CSG::UNION);
//
    Carve_volume result_volume (result);
    result_volume.hasColor = this->Volumes[0].hasColor;
    result_volume.color = this->Volumes[0].color;

    this->Volumes.clear();
    this->Volumes.push_back (result_volume);
#endif // 0
    CSGIF_polyhedron *result = new CSGIF_polyhedron();

    // a.clipto(b)
    for (size_t i=0; i < this->Volumes.size(); i++) {

        Carve_volume *new_vol = new Carve_volume (this->Volumes[i]);
        //
        if (new_vol->poly->meshes.size() > 0)
            result->Volumes.push_back (*new_vol);
    }

    // b.clipto(a)
    for (size_t j=0; j < other.Volumes.size(); j++) {
        Carve_volume *new_vol = new Carve_volume (other.Volumes[j]);

        for (size_t i=0; i < this->Volumes.size(); i++) {
            new_vol->poly = carve::csg::CSG().compute(new_vol->poly, this->Volumes[i].poly, carve::csg::CSG::A_MINUS_B);
        }
        if (new_vol->poly->meshes.size() > 0)
            result->Volumes.push_back (*new_vol);
    }

    this->Volumes = result->Volumes;

    return *this;
}

CSGIF_polyhedron& CSGIF_polyhedron::operator*=(const CSGIF_polyhedron &other)
{
    PRINT ("csgif intersection");
//	(*this->poly) *= (*other.poly);
#if 0
    carve::mesh::MeshSet<3> * mesh1 = this->poly->poly;
    carve::mesh::MeshSet<3> * mesh2 = other.poly->poly;
    carve::mesh::MeshSet<3> *result = NULL;
    result = carve::csg::CSG().compute(mesh1, mesh2, carve::csg::CSG::INTERSECTION, NULL, carve::csg::CSG::CLASSIFY_NORMAL);
    this->poly->poly = result;
#endif
    CSGIF_polyhedron *result = new CSGIF_polyhedron();

    for (size_t i=0; i < this->Volumes.size(); i++) {
        for (size_t j=0; j < other.Volumes.size(); j++) {
            Carve_volume *new_vol = new Carve_volume (this->Volumes[i]);

            new_vol->poly = carve::csg::CSG().compute(new_vol->poly, other.Volumes[j].poly, carve::csg::CSG::INTERSECTION);
            if (new_vol->poly->meshes.size() > 0)
                result->Volumes.push_back (*new_vol);
        }
    }

    this->Volumes = result->Volumes;

	return *this;
}

CSGIF_polyhedron& CSGIF_polyhedron::operator-=(const CSGIF_polyhedron &other)
{
    PRINT ("csgif difference");
#if 0
    carve::mesh::MeshSet<3> * mesh1 = this->poly->poly;
    carve::mesh::MeshSet<3> * mesh2 = other.poly->poly;
    carve::mesh::MeshSet<3> *result = NULL;
    result = carve::csg::CSG().compute(mesh1, mesh2, carve::csg::CSG::A_MINUS_B, NULL, carve::csg::CSG::CLASSIFY_NORMAL);
    this->poly->poly = result;
#endif
    CSGIF_polyhedron *result = new CSGIF_polyhedron();

    for (size_t i=0; i < this->Volumes.size(); i++) {
        Carve_volume *my_vol = new Carve_volume (this->Volumes[i]);

        for (size_t j=0; j < other.Volumes.size(); j++) {
            my_vol->poly = carve::csg::CSG().compute(my_vol->poly, other.Volumes[j].poly, carve::csg::CSG::A_MINUS_B);
        }

        result->Volumes.push_back (*my_vol);
    }

    this->Volumes = result->Volumes;

	return *this;
}

CSGIF_polyhedron &CSGIF_polyhedron::minkowski(const CSGIF_polyhedron &other)
{
	//(*this->poly) = CGAL::minkowski_sum_3(*this->poly, *other.poly);
	//todo?
	return *this;
}

size_t CSGIF_polyhedron::memsize() const
{
	if (this->isEmpty()) return 0;

//todo
	size_t memsize = sizeof(CSGIF_polyhedron);
	//memsize += this->poly->bytes();
	return memsize;
}

bool CSGIF_polyhedron::isEmpty() const
{
    //todo
	//return !this->poly || this->poly->is_empty();
	return this->Volumes.size() == 0;
}


void CSGIF_polyhedron::resize(Vector3d newsize,
			const Eigen::Matrix<bool,3,1> &autosize)
{
    //todo
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

std::string CSGIF_polyhedron::dump() const
{
//	return OpenSCAD::dump_svg( *this->poly );

    std::cout << "------------------------------------------------------\n";
    std::cout << "Num Volumes " << this->Volumes.size() << "\n";
    for (size_t i=0; i < this->Volumes.size(); i++) {
        std::cout << "- Volume " << i << " -" <<  "\n";
        std::cout << "  hasColor " << Volumes[i].hasColor <<  "\n";
        std::cout << "  color " << Volumes[i].color[0] << " " << Volumes[i].color[1] << " " << Volumes[i].color[2] << " " <<  "\n";
        csgif_dump_meshset(Volumes[i].poly);
    }
	return "";
}


void CSGIF_polyhedron::transform( const Transform3d &matrix )
{
	if (!this->isEmpty()) {
		if (matrix.matrix().determinant() == 0) {
			PRINT("WARNING: Scaling a 3D object with 0 - removing object");
			this->reset();
		}
		else {
            carve::math::Matrix transform_matrix (
                matrix(0,0), matrix(0,1), matrix(0,2), matrix(0,3),
                matrix(1,0), matrix(1,1), matrix(1,2), matrix(1,3),
                matrix(2,0), matrix(2,1), matrix(2,2), matrix(2,3),
                matrix(3,0), matrix(3,1), matrix(3,2), matrix(3,3)
                );

            for (size_t i=0; i < this->Volumes.size(); i++) {
                this->Volumes[i].poly->transform (carve::math::matrix_transformation(transform_matrix));
            }
		}
	}
}
