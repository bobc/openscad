#pragma once

#include <stddef.h>
#include <string>
#include <list>

#include "linalg.h"
#include "memory.h"

class Geometry
{
public:
  typedef std::pair<const class AbstractNode *, shared_ptr<const Geometry> > ChildItem;
  typedef std::list<ChildItem> ChildList;

	Geometry() : convexity(1) { hasColor = false; color[0] = color[1] = color[2] = 0.5f;}
	virtual ~Geometry() {}

	virtual size_t memsize() const = 0;
	virtual BoundingBox getBoundingBox() const = 0;
	virtual std::string dump() const = 0;
	virtual unsigned int getDimension() const = 0;
	virtual bool isEmpty() const = 0;
	virtual Geometry *copy() const = 0;

	unsigned int getConvexity() const { return convexity; }
	void setConvexity(int c) { this->convexity = c; }

	virtual void setColor(const Color4f &c) { color = c; hasColor = true; }

    bool hasColor;
    Color4f color;
protected:
	int convexity;

};
