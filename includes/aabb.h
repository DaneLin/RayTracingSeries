#pragma once
#ifndef AABB_H
#define AABB_H

#include "rtweekend.h"

class aabb
{
public:
	interval x, y, z;

	aabb () {}

	aabb(const point3& a, const point3& b)
	{
		x = interval(fmin(a[0], b[0]), fmax(a[0], b[0]));
		y = interval(fmin(a[1], b[1]), fmax(a[1], b[1]));
		z = interval(fmin(a[2], b[2]), fmax(a[2], b[2]));
	}

	aabb(const aabb& box0, const aabb& box1)
	{
		x = interval(box0.x, box1.x);
		y = interval(box0.y, box1.y);
		z = interval(box0.z, box1.z);
	}

	aabb(const interval& _x, const interval& _y, const interval& _z)
	{
		x = _x;
		y = _y;
		z = _z;
	}

	aabb pad()
	{
		// return an aabb that has no side narrower than some delta, padding if necessary
		double delta = 0.0001;
		interval new_x = (x.size() >= delta) ? x : x.expand(delta);
		interval new_y = (y.size() >= delta) ? y : y.expand(delta);
		interval new_z = (z.size() >= delta) ? z : z.expand(delta);

		return aabb(new_x, new_y, new_z);
	}

	const interval& axis(int n) const
	{
		if (n == 1) return y;
		if (n == 2) return z;
		return x;
	}

	/*bool hit(const ray& r, interval ray_t) const
	{
		for (int a = 0; a < 3; a++)
		{
			auto t0 = fmin((axis(a).min - r.origin()[a]) / r.direction()[a],
				(axis(a).max - r.origin()[a]) / r.direction()[a]);
			auto t1 = fmax((axis(a).min - r.origin()[a]) / r.direction()[a],
				(axis(a).max - r.direction()[a]) / r.direction()[a]);

			ray_t.min = fmax(t0, ray_t.min);
			ray_t.max = fmin(t1, ray_t.max);
			if (ray_t.max <= ray_t.min)
				return false;
		}
		return true;
	}*/
	// Andrew Kensler at Pixar 
	bool hit(const ray& r, interval ray_t) const {
		for (int a = 0; a < 3; a++) {
			auto invD = 1 / r.direction()[a];
			auto orig = r.origin()[a];

			auto t0 = (axis(a).min - orig) * invD;
			auto t1 = (axis(a).max - orig) * invD;

			if (invD < 0)
				std::swap(t0, t1);

			if (t0 > ray_t.min) ray_t.min = t0;
			if (t1 < ray_t.max) ray_t.max = t1;

			if (ray_t.max <= ray_t.min)
				return false;
		}
		return true;
	}
};

aabb operator+(const aabb& bbox, const vec3& offset)
{
	return aabb(bbox.x + offset.x(), bbox.y + offset.y(), bbox.z + offset.z());
}

aabb operator+(const vec3& offset, const aabb& bbox) {
	return bbox + offset;
}

#endif