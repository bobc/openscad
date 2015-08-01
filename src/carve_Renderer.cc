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

#ifdef _MSC_VER
// Boost conflicts with MPFR under MSVC (google it)
#include <mpfr.h>
#endif

// dxfdata.h must come first for Eigen SIMD alignment issues
#include "dxfdata.h"
#include "polyset.h"
#include "polyset-utils.h"
#include "printutils.h"

#include "CSGIF.h"

//#include "Preferences.h"

CSGIF_Renderer::CSGIF_Renderer(shared_ptr<const class Geometry> geom)
{
	if (shared_ptr<const PolySet> ps = dynamic_pointer_cast<const PolySet>(geom)) {
		assert(ps->getDimension() == 3);
		// We need to tessellate here, in case the generated PolySet contains concave polygons
        // See testdata/scad/3D/features/polyhedron-concave-test.scad
		PolySet *ps_tri = new PolySet(3, ps->convexValue());
		ps_tri->setConvexity(ps->getConvexity());
		PolysetUtils::tessellate_faces(*ps, *ps_tri);
		this->polyset.reset(ps_tri);
	}
	else if (shared_ptr<const Polygon2d> poly = dynamic_pointer_cast<const Polygon2d>(geom)) {
		this->polyset.reset(poly->tessellate());
	}
	else if (shared_ptr<const CSGIF_polyhedron> new_N = dynamic_pointer_cast<const CSGIF_polyhedron>(geom)) {
		assert(new_N->getDimension() == 3);
		if (!new_N->isEmpty()) {
			this->N = new_N;
		}
	}
}

CSGIF_Renderer::~CSGIF_Renderer()
{
}

void CSGIF_Renderer::buildPolyhedron() const
{
	PRINTD("buildPolyhedron");
	PRINTD("buildPolyhedron() end");
}

// Overridden from Renderer
void CSGIF_Renderer::setColorScheme(const ColorScheme &cs)
{
	PRINTD("setColorScheme");
	Renderer::setColorScheme(cs);
#if 0
	this->polyhedron.reset(); // Mark as dirty
#endif // 0
	PRINTD("setColorScheme done");
}

//
// triangle drawing functions copied from polyset-gl.cc
//
static void draw_triangle(GLint *shaderinfo, const Vector3d &p0, const Vector3d &p1, const Vector3d &p2, double e0f, double e1f, double e2f, double z,
    bool mirror)
{
		glVertexAttrib3d(shaderinfo[3], e0f, e1f, e2f);
		glVertexAttrib3d(shaderinfo[4], p1[0], p1[1], p1[2] + z);
		glVertexAttrib3d(shaderinfo[5], p2[0], p2[1], p2[2] + z);
		glVertexAttrib3d(shaderinfo[6], 0.0, 1.0, 0.0);
		glVertex3d(p0[0], p0[1], p0[2] + z);
		if (!mirror) {
			glVertexAttrib3d(shaderinfo[3], e0f, e1f, e2f);
			glVertexAttrib3d(shaderinfo[4], p0[0], p0[1], p0[2] + z);
			glVertexAttrib3d(shaderinfo[5], p2[0], p2[1], p2[2] + z);
			glVertexAttrib3d(shaderinfo[6], 0.0, 0.0, 1.0);
			glVertex3d(p1[0], p1[1], p1[2] + z);
		}
		glVertexAttrib3d(shaderinfo[3], e0f, e1f, e2f);
		glVertexAttrib3d(shaderinfo[4], p0[0], p0[1], p0[2] + z);
		glVertexAttrib3d(shaderinfo[5], p1[0], p1[1], p1[2] + z);
		glVertexAttrib3d(shaderinfo[6], 1.0, 0.0, 0.0);
		glVertex3d(p2[0], p2[1], p2[2] + z);
		if (mirror) {
			glVertexAttrib3d(shaderinfo[3], e0f, e1f, e2f);
			glVertexAttrib3d(shaderinfo[4], p0[0], p0[1], p0[2] + z);
			glVertexAttrib3d(shaderinfo[5], p2[0], p2[1], p2[2] + z);
			glVertexAttrib3d(shaderinfo[6], 0.0, 0.0, 1.0);
			glVertex3d(p1[0], p1[1], p1[2] + z);
		}

}

static void draw_tri(const Vector3d &p0, const Vector3d &p1, const Vector3d &p2, double z, bool mirror)
{
		glVertex3d(p0[0], p0[1], p0[2] + z);
		if (!mirror)
			glVertex3d(p1[0], p1[1], p1[2] + z);
		glVertex3d(p2[0], p2[1], p2[2] + z);
		if (mirror)
			glVertex3d(p1[0], p1[1], p1[2] + z);

}

static void gl_draw_triangle(GLint *shaderinfo, const Vector3d &p0, const Vector3d &p1, const Vector3d &p2, bool e0, bool e1, bool e2, double z, bool mirrored)
{
	double ax = p1[0] - p0[0], bx = p1[0] - p2[0];
	double ay = p1[1] - p0[1], by = p1[1] - p2[1];
	double az = p1[2] - p0[2], bz = p1[2] - p2[2];
	double nx = ay*bz - az*by;
	double ny = az*bx - ax*bz;
	double nz = ax*by - ay*bx;
	double nl = sqrt(nx*nx + ny*ny + nz*nz);

	glNormal3d(nx / nl, ny / nl, nz / nl);
#ifdef ENABLE_OPENCSG
	if (shaderinfo) {
		double e0f = e0 ? 2.0 : -1.0;
		double e1f = e1 ? 2.0 : -1.0;
		double e2f = e2 ? 2.0 : -1.0;
		draw_triangle(shaderinfo, p0, p1, p2, e0f, e1f, e2f, z, mirrored);
	}
	else
#endif
	{
		draw_tri(p0, p1, p2, z, mirrored);
	}
}

double g_scale = 1.0;

static inline void glVertex(const carve::geom3d::Vector &v) {
  glVertex3f(g_scale * (v.x),
             g_scale * (v.y),
             g_scale * (v.z));
}

static void drawFaceNormal(carve::mesh::Face<3> *face, float r, float g, float b) {
  glDisable(GL_LIGHTING);
  glDepthMask(GL_FALSE);

  glLineWidth(2.0f);

  glEnable(GL_DEPTH_TEST);

  glBegin(GL_LINES);
  glColor4f(r,g,b, 1.0);
  glVertex(face->centroid());
  glColor4f(r,g,b, 1.0);
  glVertex(face->centroid() + 1 / g_scale * face->plane.N);
  glEnd();

  glDepthMask(GL_TRUE);
  glEnable(GL_LIGHTING);
}

static void draw_volume (const Carve_volume &volume, const ColorScheme *colorscheme, bool showfaces, bool showedges)
{
    carve::mesh::MeshSet<3> * poly = volume.poly;
    GLint *shaderinfo = NULL;
    bool mirrored = false;
#if 0
    for (size_t i =0; i < poly->vertex_storage.size(); i++) {
        output << "    <vertex><coordinates>\r\n";
        output << "     <x>" << poly->vertex_storage[i].v[0] << "</x>\r\n";
        output << "     <y>" << poly->vertex_storage[i].v[1] << "</y>\r\n";
        output << "     <z>" << poly->vertex_storage[i].v[2] << "</z>\r\n";
        output << "    </coordinates></vertex>\r\n";
    }
#endif // 0

    if (showfaces)
    {
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);

        if (volume.hasColor)
        {
            glColor3f (volume.color[0], volume.color[1], volume.color[2]);
        }
        else
        {
            Color4f tcolor = ColorMap::getColor(*colorscheme, CGAL_FACE_FRONT_COLOR);
            glColor3f (tcolor[0], tcolor[1], tcolor[2]);
        }

        glBegin(GL_TRIANGLES);

        for (carve::mesh::MeshSet<3>::face_iter i = poly->faceBegin(); i != poly->faceEnd(); ++i) {
            carve::mesh::MeshSet<3>::face_t *f = *i;

            carve::mesh::MeshSet<3>::face_t::edge_iter_t edge = f->begin();
            int v1,v2,v3;

            v1 = (*edge++).vert - &poly->vertex_storage[0];
            v3 = (*edge++).vert - &poly->vertex_storage[0];

            do {
                v2 = v3;
                v3 = (*edge++).vert - &poly->vertex_storage[0];

                Vector3d p0 (poly->vertex_storage[v1].v[0],
                               poly->vertex_storage[v1].v[1],
                               poly->vertex_storage[v1].v[2]
                               );

                Vector3d p1 (poly->vertex_storage[v2].v[0],
                               poly->vertex_storage[v2].v[1],
                               poly->vertex_storage[v2].v[2]
                               );

                Vector3d p2 (poly->vertex_storage[v3].v[0],
                               poly->vertex_storage[v3].v[1],
                               poly->vertex_storage[v3].v[2]
                               );

                gl_draw_triangle(shaderinfo, p0, p1, p2, true, true, true, 0, mirrored);

            } while (edge != f->end());
        }

        glEnd();
    }

    // normals
//#define SHOW_NORMALS
#ifdef SHOW_NORMALS
        for (carve::mesh::MeshSet<3>::face_iter i = poly->faceBegin(); i != poly->faceEnd(); ++i) {
            carve::mesh::MeshSet<3>::face_t *f = *i;
            drawFaceNormal(f, 0,0,0);
        }
#endif

    if (showedges)
    {
        glLineWidth(5);
        //glLineWidth(1);

        Color4f color = ColorMap::getColor(*colorscheme, CGAL_EDGE_FRONT_COLOR);
        glColor3f (color[0], color[1], color[2]);

#if 0
        for (carve::mesh::MeshSet<3>::face_iter i = poly->faceBegin(); i != poly->faceEnd(); ++i) {
            carve::mesh::MeshSet<3>::face_t *face = *i;

            glBegin(GL_LINE_LOOP);
            for (carve::mesh::MeshSet<3>::face_t::edge_iter_t edge = face->begin(); edge != face->end(); ++edge)
            {
                int v = (*edge).vert - &poly->vertex_storage[0];

                glVertex3d  (poly->vertex_storage[v].v[0],
                               poly->vertex_storage[v].v[1],
                               poly->vertex_storage[v].v[2]
                               );

                glVertex3d (edge->))
            }
            glEnd();
#endif // 0

#if 1
        for (size_t k = 0; k < poly->meshes.size(); ++k) {

            carve::mesh::Mesh<3> *mesh = poly->meshes[k];

            for (size_t i = 0, l = mesh->faces.size(); i != l; ++i) {
                carve::mesh::Face<3> *face = mesh->faces[i];

                glBegin(GL_LINES);
                carve::mesh::Edge<3> *e = face->edge;

                for (size_t j = 0, l2 = face->nVertices(); j != l2; ++j) {
                    glVertex(e->v1()->v);
                    glVertex(e->v2()->v);
                    e = e->next;
                }
                glEnd();
            }
#endif
        }
    }
}

void CSGIF_Renderer::draw(bool showfaces, bool showedges) const
{
	PRINTD("draw()");
	if (this->polyset) {
		PRINTD("draw() polyset");
		if (this->polyset->getDimension() == 2) {
			// Draw 2D polygons
			glDisable(GL_LIGHTING);
// FIXME:		const QColor &col = Preferences::inst()->color(Preferences::CGAL_FACE_2D_COLOR);
            glColor3f(0.0f, 0.75f, 0.60f);

			for (size_t i=0; i < this->polyset->polygons.size(); i++) {
				glBegin(GL_POLYGON);
				for (size_t j=0; j < this->polyset->polygons[i].size(); j++) {
					const Vector3d &p = this->polyset->polygons[i][j];
					glVertex3d(p[0], p[1], -0.1);
				}
				glEnd();
			}

			// Draw 2D edges
			glDisable(GL_DEPTH_TEST);

			glLineWidth(2);
// FIXME:		const QColor &col2 = Preferences::inst()->color(Preferences::CGAL_EDGE_2D_COLOR);
			glColor3f(1.0f, 0.0f, 0.0f);
			this->polyset->render_edges(CSGMODE_NONE);
			glEnable(GL_DEPTH_TEST);
		}
		else {
			// Draw 3D polygons
			const Color4f c(-1,-1,-1,-1);
			setColor(COLORMODE_MATERIAL, c.data(), NULL);

            if (this->polyset->hasColor)
            {
                glColor3f (this->polyset->color[0], this->polyset->color[1], this->polyset->color[2]);
            }
            else
            {
                Color4f tcolor = ColorMap::getColor(*colorscheme, CGAL_FACE_FRONT_COLOR);
                glColor3f (tcolor[0], tcolor[1], tcolor[2]);
            }

			this->polyset->render_surface(CSGMODE_NORMAL, Transform3d::Identity(), NULL);
		}
	}
	else {
        for (size_t i=0; i < N->Volumes.size(); i++) {
            draw_volume (N->Volumes[i], colorscheme, showfaces, showedges);
        }
	}
	PRINTD("draw() end");
}



BoundingBox CSGIF_Renderer::getBoundingBox() const
{
	BoundingBox bbox;

	if (this->polyset) {
		bbox = this->polyset->getBoundingBox();
	}
	else {
        if (N)
        {
            carve::geom3d::AABB aabb = N->Volumes[0].poly->getAABB();
            carve::geom3d::Vector min, max;
            min = aabb.min();
            max = aabb.max();

			bbox = BoundingBox(
				Vector3d(min.x, min.y, min.z) ,
				Vector3d(max[0], max[1], max[2]) );
        }
	}
	return bbox;
}
