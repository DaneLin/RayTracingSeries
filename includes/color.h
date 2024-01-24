#pragma once
#ifndef COLOR_H
#define COLOR_H

#include "vec3.h"

#include <iostream>

using color = vec3;

inline double linear_to_gamma(double linear_component)
{
	return sqrt(linear_component);
}

void write_color(std::ostream& out, color pixel_color, int samplers_per_pixel) {
	auto r = pixel_color.x();
	auto g = pixel_color.y();
	auto b = pixel_color.z();

	// Replace NaN components with zero
	if (r != r) r = 0.0;
	if (g != g) g = 0.0;
	if (b != b) b = 0.0;

	// Divide the color by the numer of samples
	auto scale = 1.0f / samplers_per_pixel;
	r *= scale;
	g *= scale;
	b *= scale;

	// Apply the linear to gamma transform.
	r = linear_to_gamma(r);
	g = linear_to_gamma(g);
	b = linear_to_gamma(b);

	static const interval intensity(0.000, 0.999);

	// Write the translated [0,255] value of each color component.
	out << static_cast<int>(256 * intensity.clamp(r)) << ' '
		<< static_cast<int>(256 * intensity.clamp(g)) << ' '
		<< static_cast<int>(256 * intensity.clamp(b)) << '\n';
}

#endif