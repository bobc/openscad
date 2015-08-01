#pragma once

#include <iostream>
#include "Tree.h"
#include "Camera.h"
#include "memory.h"

enum FileFormat {
	OPENSCAD_STL,
	OPENSCAD_OFF,
	OPENSCAD_AMF,
	OPENSCAD_DXF,
	OPENSCAD_SVG
};

#ifdef ENABLE_CSGIF
void export_stl(const class CSGIF_polyhedron *root_N, std::ostream &output);
void export_off(const CSGIF_polyhedron *root_N, std::ostream &output);
void export_amf(const class CSGIF_polyhedron *root_N, std::ostream &output);
void export_png(const CSGIF_polyhedron *root_N, Camera &c, std::ostream &output);

void export_stl_files(const class CSGIF_polyhedron *root_N, const char *filename);

#endif // ENABLE_CSGIF

// void exportFile(const class Geometry *root_geom, std::ostream &output, FileFormat format);
void exportFileByName(const class Geometry *root_geom, FileFormat format,
	const char *name2open, const char *name2display);
void export_png(shared_ptr<const class Geometry> root_geom, Camera &c, std::ostream &output);

void export_stl(const class PolySet &ps, std::ostream &output);
void export_off(const class PolySet &ps, std::ostream &output);
void export_amf(const class PolySet &ps, std::ostream &output);
void export_dxf(const class Polygon2d &poly, std::ostream &output);
void export_svg(const class Polygon2d &poly, std::ostream &output);

void export_png_with_opencsg(Tree &tree, Camera &c, std::ostream &output);
void export_png_with_throwntogether(Tree &tree, Camera &c, std::ostream &output);


#ifdef DEBUG
void export_stl(const class PolySet &ps, std::ostream &output);
#endif
