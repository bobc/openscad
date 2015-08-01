/*
 *  OpenSCAD (www.openscad.org)
 *  Copyright (C) 2009-2011 Clifford Wolf <clifford@clifford.at> and
 *                          Marius Kintel <marius@kintel.net>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  As a special exception, you have permission to link this program
 *  with the CGAL library and distribute executables, as long as you
 *  follow the requirements of the GNU GPL in regard to all of the
 *  software in the executable aside from CGAL.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "export.h"
#include "printutils.h"
#include "polyset.h"
#include "polyset-utils.h"
#include "dxfdata.h"

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/convenience.hpp>

#define QUOTE(x__) # x__
#define QUOTED(x__) QUOTE(x__)

#include "CSGIF.h"
#include "IndexedMesh.h"

#include <fstream>

struct triangle {
    std::string vs1;
    std::string vs2;
    std::string vs3;
};

static std::vector<Color4f> find_unique_colors (const CSGIF_polyhedron *P)
{
    std::vector<Color4f> result;

    for (size_t vol=0; vol < P->Volumes.size(); vol++) {
        if (P->Volumes[vol].hasColor) {

            if (std::find(result.begin(), result.end(), P->Volumes[vol].color) == result.end()) {
                result.push_back(P->Volumes[vol].color);
            }
        }
    }
    return result;
}

static IndexedMesh *get_shape (const Carve_volume &volume, bool select_color, Color4f color)
{
    IndexedMesh *mesh = new IndexedMesh();
    bool pick;

    pick = false;
    if (!select_color)
        pick = true;
    else if (volume.hasColor)
    {
        if (volume.color == color)
            pick = true;
    }

    if (pick)
    {
        // todo: get triangles

        carve::mesh::MeshSet<3> * poly = volume.poly;

        for (size_t i=0; i < poly->vertex_storage.size(); i++) {
                Vector3d vert (poly->vertex_storage[i].v[0],
                               poly->vertex_storage[i].v[1],
                               poly->vertex_storage[i].v[2] );
                mesh->polys.vertices.push_back(vert);
        }

        //
        for (carve::mesh::MeshSet<3>::face_iter i = poly->faceBegin(); i != poly->faceEnd(); ++i) {
            carve::mesh::MeshSet<3>::face_t *f = *i;

            carve::mesh::MeshSet<3>::face_t::edge_iter_t edge = f->begin();
            int v1,v2,v3;

            v1 = (*edge).vert - &poly->vertex_storage[0];
            edge++;
            v3 = (*edge).vert - &poly->vertex_storage[0];
            edge++;

            do {
                v2 = v3;
                v3 = (*edge).vert - &poly->vertex_storage[0];
                edge++;

                //
                IndexedFace face;
                face.push_back(v1);
                face.push_back(v2);
                face.push_back(v3);

                mesh->polys.faces.push_back(face);

            } while (edge != f->end());

            //
        }

    }
    return mesh;
}

static IndexedMesh *get_shape (const Carve_volume &volume)
{
    return get_shape (volume, false, Color4f(0,0,0));
}

static void amf_output_mesh (std::vector<IndexedMesh> amfVolumes, std::ostream &output)
{
	output << "  <mesh>\r\n";

    output << "   <vertices>\r\n";

    std::vector<int> base_index;
    int index = 0;
    base_index.push_back(index);

    for (size_t vol=0; vol < amfVolumes.size(); vol++) {

        index += amfVolumes[vol].polys.vertices.size();
        base_index.push_back(index);

        for (size_t i=0; i < amfVolumes[vol].polys.vertices.size(); i++) {
            output << "    <vertex><coordinates>\r\n";
            output << "     <x>" << amfVolumes[vol].polys.vertices[i][0] << "</x>\r\n";
            output << "     <y>" << amfVolumes[vol].polys.vertices[i][1] << "</y>\r\n";
            output << "     <z>" << amfVolumes[vol].polys.vertices[i][2] << "</z>\r\n";
            output << "    </coordinates></vertex>\r\n";
        }
    }
    output << "   </vertices>\r\n";

    // volumes...
    for (size_t vol=0; vol < amfVolumes.size(); vol++) {

        if (amfVolumes[vol].hasMaterialId) {
            output << "   <volume materialid=\"" << amfVolumes[vol].materialId << "\">\r\n";
        }
        else
        {
            output << "   <volume>\r\n";
            if (amfVolumes[vol].hasColor)
            {
                output << "   <color>" <<
                "<r>" << amfVolumes[vol].color[0] << "</r>" <<
                "<g>" << amfVolumes[vol].color[1] << "</g>" <<
                "<b>" << amfVolumes[vol].color[2] << "</b>" <<
                "</color>\r\n";
            }
        }

        for (size_t i=0; i < amfVolumes[vol].polys.faces.size(); i++) {
            output << "    <triangle>\r\n";
            output << "     <v1>" << base_index[vol] + amfVolumes[vol].polys.faces[i][0] << "</v1>\r\n";
            output << "     <v2>" << base_index[vol] + amfVolumes[vol].polys.faces[i][1] << "</v2>\r\n";
            output << "     <v3>" << base_index[vol] + amfVolumes[vol].polys.faces[i][2] << "</v3>\r\n";
            output << "    </triangle>\r\n";
        }

        output << "   </volume>\r\n";
    }

    output << "  </mesh>\r\n";
 }

/*!
	Saves the current 3D CSG polyhedron as STL to the given file.
	The file must be open.
 */
void export_stl(const CSGIF_polyhedron *root_N, std::ostream &output)
{
}

void export_off(const CSGIF_polyhedron *root_N, std::ostream &output)
{
}

static void export_stl (const IndexedMesh &mesh, std::ostream &output)
{
  	setlocale(LC_NUMERIC, "C"); // Ensure radix is . (not ,) in output

	output << "solid OpenSCAD_Model\n";

    BOOST_FOREACH (const IndexedFace &face, mesh.polys.faces) {
        output << "  facet normal ";

        Polygon p;
        p.push_back(mesh.polys.vertices[face[0]]);
        p.push_back(mesh.polys.vertices[face[1]]);
        p.push_back(mesh.polys.vertices[face[2]]);

        Vector3d normal = (p[1] - p[0]).cross(p[2] - p[0]);
        normal.normalize();
        if (is_finite(normal) && !is_nan(normal)) {
            output << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
        }
        else {
            output << "0 0 0\n";
        }
        output << "    outer loop\n";

        BOOST_FOREACH(const Vector3d &v, p) {
            output << "      vertex " << v[0] << " " << v[1] << " " << v[2] << "\n";
        }
        output << "    endloop\n";
        output << "  endfacet\n";
    }
  	output << "endsolid OpenSCAD_Model\n";
	setlocale(LC_NUMERIC, "");      // Set default locale
}

static void stl_output_mesh_file (const IndexedMesh &mesh, std::string name2open)
{
    std::string name2display = name2open;    //TODO

    std::ofstream fstream(name2open.c_str());
    if (!fstream.is_open()) {
        PRINTB(_("Can't open file \"%s\" for export"), name2display);
    } else {
        bool onerror = false;
        fstream.exceptions(std::ios::badbit|std::ios::failbit);
        try {
            PRINTB ("stl filename : %s", name2display.c_str());

            export_stl (mesh, fstream);
        } catch (std::ios::failure x) {
            onerror = true;
        }
        try { // make sure file closed - resources released
            fstream.close();
        } catch (std::ios::failure x) {
            onerror = true;
        }
        if (onerror) {
            PRINTB(_("ERROR: \"%s\" write error. (Disk full?)"), name2display);
        }
    }
}

void export_stl_files(const class CSGIF_polyhedron *root_N, const char *filename)
{
    // output an STL file for each unique color
    std::vector<Color4f> unique_colors = find_unique_colors (root_N);
    std::string base_filename = boost::filesystem::change_extension(std::string(filename), "").string();
    int fileNumber = 0;

    std::vector<IndexedMesh> amfVolumes;

    for (size_t color_index=0; color_index < unique_colors.size(); color_index++) {
        IndexedMesh mesh;

        mesh.hasColor = true;
        mesh.color = unique_colors [color_index];

        std::string output_name = base_filename;
        output_name.append ("_");
        output_name.append ( str (boost::format ("%1%") % fileNumber ) );
        output_name = boost::filesystem::change_extension(output_name, ".stl").string();
        fileNumber++;

        for (size_t vol=0; vol < root_N->Volumes.size(); vol++) {

            if (mesh.color == root_N->Volumes[vol].color)
            {
                IndexedMesh *mesh2 = get_shape (root_N->Volumes[vol], true, mesh.color);
                mesh.AddMesh (*mesh2);
            }
        }

        stl_output_mesh_file (mesh, output_name);
    }
}

/*!
    Saves the current 3D CSG polyhedron as AMF to the given file.
    The file must be open.
 */
void export_amf(const CSGIF_polyhedron *root_N, std::ostream &output)
{
    std::vector<Color4f> unique_colors = find_unique_colors (root_N);

    output << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\r\n"
			<< "<amf unit=\"millimeter\">\r\n"
			<< " <metadata type=\"producer\">OpenSCAD " << QUOTED(OPENSCAD_VERSION)
#ifdef OPENSCAD_COMMIT
			<< " (git " << QUOTED(OPENSCAD_COMMIT) << ")"
#endif
			<< "</metadata>\r\n" ;



//
    // TODO: handle "no color" volumes

#define MATERIAL_PER_COLOR

#ifdef BY_VOLUME
    // output a single object with one volume for each geom volume
    // each volume with have a color according to the geom volume properties
	output << " <object id=\"0\">\r\n";
    std::vector<IndexedMesh> amfVolumes;
    for (size_t vol=0; vol < root_N->Volumes.size(); vol++) {
        IndexedMesh *mesh = get_shape (root_N->Volumes[vol]);
        mesh->hasColor = root_N->Volumes[vol].hasColor;
        mesh->color = root_N->Volumes[vol].color;

        amfVolumes.push_back(*mesh);
    }
    amf_output_mesh (amfVolumes, output);
	output << " </object>\r\n";
#endif // BY_VOLUME

#ifdef MATERIAL_PER_COLOR
    // output a single object with one volume for each unique color
    // and a material id specifier
	output << " <object id=\"0\">\r\n";
    std::vector<IndexedMesh> amfVolumes;

    int materialId = 1;

    for (size_t color_index=0; color_index < unique_colors.size(); color_index++) {
        IndexedMesh mesh;

        mesh.hasMaterialId = true;
        mesh.materialId = materialId++;
        mesh.hasColor = true;
        mesh.color = unique_colors [color_index];

        for (size_t vol=0; vol < root_N->Volumes.size(); vol++) {

            if (mesh.color == root_N->Volumes[vol].color)
            {
                IndexedMesh *mesh2 = get_shape (root_N->Volumes[vol], true, mesh.color);
                mesh.AddMesh (*mesh2);
            }
        }
        amfVolumes.push_back(mesh);
    }
    amf_output_mesh (amfVolumes, output);
	output << " </object>\r\n";

    // materials
    for (size_t vol=0; vol < amfVolumes.size(); vol++) {
        output << " <material id=\"" << amfVolumes[vol].materialId << "\">\r\n";
        output << "  <color>" <<
            "<r>" << amfVolumes[vol].color[0] << "</r>" <<
            "<g>" << amfVolumes[vol].color[1] << "</g>" <<
            "<b>" << amfVolumes[vol].color[2] << "</b>" <<
            "</color>\r\n";
        output << " </material>\r\n";
    }


#endif // MATERIAL_PER_COLOR

    // todo : dispose
//
#if 0
    std::vector<int> base_index;
    int index = 0;
    base_index.push_back(index);

    for (size_t vol=0; vol < root_N->Volumes.size(); vol++) {

        carve::mesh::MeshSet<3> * poly = root_N->Volumes[vol].poly;
        index += poly->vertex_storage.size();
        base_index.push_back(index);

        for (size_t i =0; i < poly->vertex_storage.size(); i++) {
                output << "    <vertex><coordinates>\r\n";
                output << "     <x>" << poly->vertex_storage[i].v[0] << "</x>\r\n";
                output << "     <y>" << poly->vertex_storage[i].v[1] << "</y>\r\n";
                output << "     <z>" << poly->vertex_storage[i].v[2] << "</z>\r\n";
                output << "    </coordinates></vertex>\r\n";
            }
    }
    output << "   </vertices>\r\n";

    for (size_t vol=0; vol < root_N->Volumes.size(); vol++) {

        carve::mesh::MeshSet<3> * poly = root_N->Volumes[vol].poly;

        output << "   <volume>\r\n";

        if (root_N->Volumes[vol].hasColor)
            output << "   <color>" <<
                "<r>" << root_N->Volumes[vol].color[0] << "</r>" <<
                "<g>" << root_N->Volumes[vol].color[1] << "</g>" <<
                "<b>" << root_N->Volumes[vol].color[2] << "</b>" <<
                "</color>\r\n"
                ;

        for (carve::mesh::MeshSet<3>::face_iter i = poly->faceBegin(); i != poly->faceEnd(); ++i) {
            carve::mesh::MeshSet<3>::face_t *f = *i;

            carve::mesh::MeshSet<3>::face_t::edge_iter_t edge = f->begin();
            int v1,v2,v3;

            v1 = (*edge).vert - &poly->vertex_storage[0];
            edge++;
            v3 = (*edge).vert - &poly->vertex_storage[0];
            edge++;

            do {
                v2 = v3;
                v3 = (*edge).vert - &poly->vertex_storage[0];
                edge++;

                //
                output << "    <triangle>\r\n";
                output << "     <v1>" << v1 + base_index [vol] << "</v1>\r\n";
                output << "     <v2>" << v2 + base_index [vol]<< "</v2>\r\n";
                output << "     <v3>" << v3 + base_index [vol]<< "</v3>\r\n";
                output << "    </triangle>\r\n";
            } while (edge != f->end());

            //
        }
        output << "   </volume>\r\n";
    }

    output << "  </mesh>\r\n"
        << " </object>\r\n"
#endif // 0

    output << "</amf>\r\n";
}


