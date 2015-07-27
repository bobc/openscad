/*
 *  OpenSCAD (www.openscad.org)
 *  Copyright (C) 2015 Bob Cousins <bobcousins42@googlemail.com>
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

#pragma once

#ifdef ENABLE_CSGIF

#include "GeometryEvaluator.h"

#include "CSGIF_Cache.h"

#endif // ENABLE_CSGIF


#ifdef ENABLE_CGAL

#define ENABLE_MINKOWSKI
#define ENABLE_HULL

#include "cgal.h"

#include <CGAL/convex_hull_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Cartesian.h>

#include "CGAL_Nef_polyhedron.h"
#include "cgalutils.h"
#include "CGALRenderer.h"

#endif // ENABLE_CGAL

#ifdef ENABLE_CARVE

#include <carve/carve.hpp>
#include <carve/csg.hpp>
#include <carve/poly.hpp>
#include <carve/geom.hpp>

#include "carve_polyhedron.h"
#include "carve_utils.h"
#include "carve_Renderer.h"

#endif // ENABLE_CARVE
