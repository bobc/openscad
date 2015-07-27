#pragma once

#include "renderer.h"

#include "CSGIF.h"


class CSGIF_Renderer : public Renderer
{
public:
	CSGIF_Renderer(shared_ptr<const class Geometry> geom);
	~CSGIF_Renderer();
	virtual void draw(bool showfaces, bool showedges) const;
	virtual void setColorScheme(const ColorScheme &cs);
	virtual BoundingBox getBoundingBox() const;

private:
	//shared_ptr<class CGAL_OGL_Polyhedron> getPolyhedron() const;
	void buildPolyhedron() const;

//	mutable shared_ptr<class CGAL_OGL_Polyhedron> polyhedron;

	shared_ptr<const CSGIF_polyhedron> N;
	shared_ptr<const class PolySet> polyset;
};
