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
#include <boost/algorithm/string.hpp>

#define QUOTE(x__) # x__
#define QUOTED(x__) QUOTE(x__)

#include "CSGIF.h"

#include <fstream>

struct triangle {
    std::string vs1;
    std::string vs2;
    std::string vs3;
};


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

/*!
    Saves the current 3D CSG polyhedron as AMF to the given file.
    The file must be open.
 */
void export_amf(const CSGIF_polyhedron *root_N, std::ostream &output)
{

    output << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\r\n"
			<< "<amf unit=\"millimeter\">\r\n"
			<< " <metadata type=\"producer\">OpenSCAD " << QUOTED(OPENSCAD_VERSION)
#ifdef OPENSCAD_COMMIT
			<< " (git " << QUOTED(OPENSCAD_COMMIT) << ")"
#endif
			<< "</metadata>\r\n"
			<< " <object id=\"0\">\r\n"
			<< "  <mesh>\r\n";
    output << "   <vertices>\r\n";

    carve::mesh::MeshSet<3> * poly = root_N->poly->poly;

    for (size_t i =0; i < poly->vertex_storage.size(); i++) {
			output << "    <vertex><coordinates>\r\n";
			output << "     <x>" << poly->vertex_storage[i].v[0] << "</x>\r\n";
			output << "     <y>" << poly->vertex_storage[i].v[1] << "</y>\r\n";
			output << "     <z>" << poly->vertex_storage[i].v[2] << "</z>\r\n";
			output << "    </coordinates></vertex>\r\n";
		}
	output << "   </vertices>\r\n";

	output << "   <volume>\r\n";
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
            output << "     <v1>" << v1 << "</v1>\r\n";
            output << "     <v2>" << v2 << "</v2>\r\n";
            output << "     <v3>" << v3 << "</v3>\r\n";
            output << "    </triangle>\r\n";
        } while (edge != f->end());

        //
    }
    output << "   </volume>\r\n";

    output << "  </mesh>\r\n"
        << " </object>\r\n"
        << "</amf>\r\n";
}


