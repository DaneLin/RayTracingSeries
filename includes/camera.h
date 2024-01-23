#pragma once
#ifndef CAMERA_H
#define CAMERA_H

#include "rtweekend.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "color.h"
#include "hittable.h"
#include "material.h"

class camera
{
public:
	double aspect_ratio = 1.0f;
	int image_width = 100;
	int samplers_per_pixel = 10;// Count of random samples for each pixel
	int max_depth = 10;// Maximum number of ray bounces into scene
	color background; // Scene background color

	double vfov = 90;//Vertical view angle(field of view
	point3 lookfrom = point3(0, 0, -1);// Point camera is looking from
	point3 lookat = point3(0, 0, 0); // Point camera is looking at
	vec3 vup = vec3(0, 1, 0); // Camera-relative "up" direction

	double defocus_angle = 0; //Variation angle of rays through each pixel
	double focus_dist = 10; //distance from camera lookfrom point to plane of perfect focus

	void render(const hittable& world)
	{
		initialize();

		// Render

		std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";

		for (int j = 0; j < image_height; ++j)
		{
			std::clog << "\rScanlines remaining: " << (image_height - j) << " " << std::flush;
			for (int i = 0; i < image_width; ++i)
			{
				
				color pixel_color(0, 0, 0);
#pragma omp parallel
{
	#pragma omp for
	for (int sample = 0; sample < samplers_per_pixel; ++sample)
	{
					ray r = get_ray(i, j);
					pixel_color += ray_color(r, max_depth,world);
	}
}
				
				write_color(std::cout, pixel_color,samplers_per_pixel);
			}
			
		}

		std::clog << "\rDone.                  \n";
	}

private:
	int    image_height;   // Rendered image height
	point3 center;         // Camera center
	point3 pixel00_loc;    // Location of pixel 0, 0
	vec3   pixel_delta_u;  // Offset to pixel to the right
	vec3   pixel_delta_v;  // Offset to pixel below
	vec3   u, v, w; // Camera frame basis vectors
	vec3   defocus_disk_u;// Defocus disk horizontal radius
	vec3   defocus_disk_v;// Defocus disk vertical radius

	void initialize()
	{
		image_height = static_cast<int>(image_width / aspect_ratio);
		image_height = (image_height < 1) ? 1 : image_height;

		center = lookfrom;

		// Determine viewport dimensions.
		//auto focal_length = (lookfrom - lookat).length();
		auto theta = degrees_to_radians(vfov);
		auto h = tan(theta / 2);
		auto viewport_height = 2.0 * h * focus_dist;
		auto viewport_width = viewport_height * (static_cast<double>(image_width) / image_height);

		// Calculate the u,v,w unit basis vectors for the camera coordinate frame
		w = unit_vector(lookfrom - lookat);
		u = unit_vector(cross(vup, w));
		v = cross(w, u);

		// Calculate the vectors across the horizontal and down the vertical viewport edges.
		auto viewport_u = viewport_width * u;;
		auto viewport_v =viewport_height * -v;

		// Calculate the horizontal and vertical delta vectors from pixel to pixel.
		pixel_delta_u = viewport_u / image_width;
		pixel_delta_v = viewport_v / image_height;

		// Calculate the location of the upper left pixel.
		auto viewport_upper_left =
			center - (focus_dist * w) - viewport_u / 2 - viewport_v / 2;
		pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

		// Calculate the camera defocus disk basis vectors
		auto defocus_radius = focus_dist * tan(degrees_to_radians(defocus_angle / 2));
		defocus_disk_u = u * defocus_radius;
		defocus_disk_v = v * defocus_radius;
	}

	color ray_color(const ray& r, int depth, const hittable& world) const
	{
		hit_record rec;
		// if we've exceeded the ray bounce limit, no more light is gathered;
		if (depth <= 0)
			return color(0, 0, 0);

		if (world.hit(r, interval(0.001, infinity), rec)) {
			/*vec3 direction = rec.normal + random_in_unit_sphere();
			return 0.5 * ray_color(ray(rec.p, direction),depth -1, world);*/
			ray scattered;
			color attenuation;

			color color_from_emission = rec.mat->emitted(rec.u, rec.v, rec.p);

			if (!rec.mat->scatter(r, rec, attenuation, scattered))
				return color_from_emission;
			
			color color_from_scatter = attenuation * ray_color(scattered, depth - 1, world);
			return color_from_emission + color_from_scatter;
		}
		return background;
	}

	ray get_ray(int i, int j) const
	{
		// Get a randomly sampled camera ray for the pixel at location i,j, originating from
		// the camera defocus disk

		auto pixel_center = pixel00_loc + (i * pixel_delta_u) + (j * pixel_delta_v);
		auto pixel_sample = pixel_center + pixel_sample_square();

		auto ray_origin = (defocus_angle <= 0) ? center : defocus_disk_sample();
		auto ray_direction = pixel_sample - ray_origin;
		auto ray_time = random_double();

		return ray(ray_origin, ray_direction,ray_time);
	}

	vec3 pixel_sample_square() const
	{
		// Returns a random point in the square surrounding a pixel at the origin
		auto px = -0.5 + random_double();
		auto py = -0.5 + random_double();
		return (px * pixel_delta_u) + (py * pixel_delta_v);
	}

	point3 defocus_disk_sample() const
	{
		// returns a random point in the camera defocus disk
		auto p = random_in_unit_disk();
		return center + (p[0] * defocus_disk_u) + (p[1] * defocus_disk_v);
	}
};

#endif