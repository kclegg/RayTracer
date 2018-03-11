#include "CImg.h"
#include <iostream>

using namespace cimg_library;

struct vector{
	float x, y, z;
};

struct rgb{
	float r, g, b;
};

struct Sphere{
	vector coordinates;
	float red;
	float green;
	float blue;
	int radius;
};

struct Tetrahedron{
	vector vertices[4];
	float r;
	float g;
	float b;
};

vector triangleNormal(vector v1, vector v2, vector v3){
	vector normal;
	vector v;				// V = Vertex2 - Vertex1
	v.x = (v2.x - v1.x);
	v.y = (v2.y - v1.y);
	v.z = (v2.z - v1.z);
	vector w;				// W = Vertex3 - Vertex1
	w.x = (v3.x - v1.x);
	w.y = (v3.y - v1.y);
	w.z = (v3.z - v1.z);
	normal.x = (v.y * w.z) - (v.z * w.y);	// Nx = (Vy * Wz) - (Vz * Wy)
	normal.y = (v.z * w.x) - (v.x * w.z);	// Ny = (Vz * Wx) - (Vx * Wz)
	normal.z = (v.x * w.y) - (v.y * w.x);	// Nz = (Vx * Wy) - (Vy * Wx)
	return normal;
}

float triangleIntersection(vector vertex, vector rayOrigin, vector normal, vector direction){
	// t = ((a - p) . n) / (d . n)
	float t;
	vector temp;
	temp.x = vertex.x-rayOrigin.x;
	temp.y = vertex.y-rayOrigin.y;
	temp.z = vertex.z-rayOrigin.z;
	float var1 = (temp.x*normal.x) + (temp.y*normal.y) + (temp.z*normal.z);
	float var2 = (direction.x*normal.x)+(direction.y*normal.y)+(direction.z*normal.z);
	float var3 = var1/var2;
	if(t < 0)
		return -1;
	return var3;
}

vector triangleXValue(float t, vector rayOrigin, vector direction){
	// x = p + td
	vector x;
	x.x = (rayOrigin.x + (direction.x * t));
	x.y = (rayOrigin.y + (direction.y * t));
	x.z = (rayOrigin.z + (direction.z * t));
	return x;
}

bool insideTriangle(vector a, vector b, vector c, vector x, vector normal){
	// ((b - a) * (x - a)) . n > 0
	// ((c - b) * (x - b)) . n > 0
	// ((a - c) * (x - c)) . n > 0

	vector ba;
	ba.x = (b.x - a.x);
	ba.y = (b.y - a.y);
	ba.z = (b.z - a.z);
	vector xa;
	xa.x = (x.x - a.x);
	xa.y = (x.y - a.y);
	xa.z = (x.z - a.z);
	vector var1;
	var1.x = (ba.x * xa.x);
	var1.y = (ba.y * xa.y);
	var1.z = (ba.z * xa.z);
	float value1 = (var1.x * normal.x) + (var1.y * normal.y) + (var1.z * normal.z);
	///////////////////////////////////////////////////////////////////////////////
	vector cb;
	cb.x = (c.x - b.x);
	cb.y = (c.y - b.y);
	cb.z = (c.z - b.z);
	vector xb;
	xb.x = (x.x - b.x);
	xb.y = (x.y - b.y);
	xb.z = (x.z - b.z);
	vector var2;
	var2.x = (cb.x * xb.x);
	var2.y = (cb.y * xb.y);
	var2.z = (cb.z * xb.z);
	float value2 = (var2.x * normal.x) + (var2.y * normal.y) + (var2.z * normal.z);
	///////////////////////////////////////////////////////////////////////////////
	vector ac;
	ac.x = (a.x - c.x);
	ac.y = (a.y - c.y);
	ac.z = (a.z - c.z);
	vector xc;
	xc.x = (x.x - c.x);
	xc.y = (x.y - c.y);
	xc.z = (x.z - c.z);
	vector var3;
	var3.x = (ac.x * xc.x);
	var3.y = (ac.y * xc.y);
	var3.z = (ac.z * xc.z);
	float value3 = (var3.x * normal.x) + (var3.y * normal.y) + (var3.z * normal.z);

	if(value1 > 0 && value2 > 0 && value3 > 0)
		return true;
	return false;
}

float rayIntersection(vector coordinates, int radius, vector rayOrigin, vector direction){
	// t < 0 (There are 0 intersection points)
	// t = 0 (There are 1 intersection points)
	// t > 0 (There are 2 intersection points)

	// t = (-d.(e-c) +- sqrt(((d.(e-c))^2) - (d.d)((e-c).(e-c) - R^2))) / (d.d)

	vector d;
	d.x = direction.x;
	d.y = direction.y;
	d.z = direction.z;

	vector c;
	c.x = coordinates.x;
	c.y = coordinates.y;
	c.z = coordinates.z;

	vector e;
	e.x = rayOrigin.x;
	e.y = rayOrigin.y;
	e.z = rayOrigin.z;

	int R = radius;

	vector vectorEC;
	vectorEC.x = (e.x - c.x);
	vectorEC.y = (e.y - c.y);
	vectorEC.z = (e.z - c.z);

	float value1 = (pow(((d.x*vectorEC.x)+(d.y*vectorEC.y)+(d.z*vectorEC.z)),2)) - (((d.x*d.x)+(d.y*d.y)+(d.z*d.z))*(((vectorEC.x*vectorEC.x)+(vectorEC.y*vectorEC.y)+(vectorEC.z*vectorEC.z)) - (R*R)));

	if (value1 < 0){
		return -1;
	}
	else{
		value1 = pow(value1,0.5);

		float t1 = ((((-d.x*vectorEC.x)+(-d.y*vectorEC.y)+(-d.z*vectorEC.z))) + (value1)) / ((d.x*d.x)+(d.y*d.y)+(d.z*d.z));
		float t2 = ((((-d.x*vectorEC.x)+(-d.y*vectorEC.y)+(-d.z*vectorEC.z))) - (value1)) / ((d.x*d.x)+(d.y*d.y)+(d.z*d.z));

		if(t1 < 0 && t2 < 0){
			return -1;
		}
		else{
			if(t1 > 0 && t2 > 0){
				if(t1 > t2){
					return t2;
				}
				return t1;
			}
			else{
				if(t1 > 0){
					return t1;
				}
				return t2;
			}
			return -1;
		}
	}
}

vector computeNormal(vector coordinates, vector rayOrigin, vector direction, float t){
	// normal vector = (intersection point - sphere center) / radius
	// normal vector = 2*(intersection point - sphere center)
	// Intersection Point = rayOrigin + direction*t

	vector intersection;
	intersection.x = (rayOrigin.x + (direction.x * t));
	intersection.y = (rayOrigin.y + (direction.y * t));
	intersection.z = (rayOrigin.z + (direction.z * t));

	vector normal;
	normal.x = (2*(intersection.x - coordinates.x));
	normal.y = (2*(intersection.y - coordinates.y));
	normal.z = (2*(intersection.z - coordinates.z));
	return normal;
}


rgb diffuseShading(float r, float g, float b, vector light, vector normal){
	rgb values;
	float value = ((normal.x * light.x) + (normal.z * light.z) + (normal.z * light.z));
	if(value < 0) {value = 0;}
	value = (value)*(0.9);			// adjust this light intensity

	values.r = (r * value);
	values.g = (g * value);
	values.b = (b * value);
	return values;
}

rgb specularShading(float r, float g, float b, vector light, vector normal, vector direction){
	// (RGB vector) Ls = I * ks * max(0,n * Vh)^n
	// Vh = (Vl + Ve) / (||(Vl + Ve)||)
	// ||vector|| = sqrt(vx*vx + vy*vy + vz*vz)

	rgb values;
	float value;

	float halfwayVectorx = (light.x + direction.x);
	float halfwayVectory = (light.y + direction.y);
	float halfwayVectorz = (light.z + direction.z);

	float temp = pow(halfwayVectorx*halfwayVectorx + halfwayVectory*halfwayVectory + halfwayVectorz*halfwayVectorz, 0.5);

	if(temp > 0){
		halfwayVectorx = halfwayVectorx / temp;
		halfwayVectory = halfwayVectory / temp;
		halfwayVectorz = halfwayVectorz / temp;
	}

	value = (normal.x * halfwayVectorx) + (normal.y * halfwayVectory) + (normal.z * halfwayVectorz);
	if(value < 0) {value = 0;}
	value = pow(value,1);			// adjust n value
	value = (value)*(0.9);			// adjust light intensity

	values.r = (r * value);
	values.g = (g * value);
	values.b = (b * value);
	return values;
}

rgb ambientShading(float r, float g, float b){
	rgb values;
	values.r = (r)*(0.9);		// adjust this light intensity
	values.g = (g)*(0.9);
	values.b = (b)*(0.9);
	return values;
}

int main(){
	// Using planar projection
	// Varying origin, constant fixed direction

	CImg<float> image(256,256,1,3);

	Sphere sphere1;
	sphere1.coordinates.x = 40;
	sphere1.coordinates.y = 40;
	sphere1.coordinates.z = 100;
	sphere1.red = 250;
	sphere1.green = 70;
	sphere1.blue = 22;
	sphere1.radius = 30;

	Sphere sphere2;
	sphere2.coordinates.x = 170;
	sphere2.coordinates.y = 150;
	sphere2.coordinates.z = 100;
	sphere2.red = 0;
	sphere2.green = 33;
	sphere2.blue = 165;
	sphere2.radius = 65;

	vector LightDir;
	LightDir.x = -2;	// -2
	LightDir.y = 40;	// 40
	LightDir.z = -1;	// -1

	vector LookAt;	// camera origin
	LookAt.x = 0;
	LookAt.y = 0;
	LookAt.z = 1;

	vector rayOrigin;

	Tetrahedron shape;
	shape.r = 255;
	shape.g = 255;
	shape.b = 51;

	shape.vertices[0].x = 60;
	shape.vertices[0].y = 170;
	shape.vertices[0].z = 70;

	shape.vertices[1].x = 60;
	shape.vertices[1].y = 225;
	shape.vertices[1].z = 125;

	shape.vertices[2].x = 115;
	shape.vertices[2].y = 225;
	shape.vertices[2].z = 125;

	shape.vertices[3].x = 115;
	shape.vertices[3].y = 170;
	shape.vertices[3].z = 125;


	// create color per pixel

	for(int x = 0; x < 256; x++){
		for(int y = 0; y < 256; y++){
			// Test Image Below - Figure (b)
			//image(x,y,0) =  y % 256;		// red		(channel 0)
			//image(x,y,1) =  x % 256;		// green	(channel 1)
			//image(x,y,2) =  x*y % 256;	// blue		(channel 2)

			rayOrigin.x = x;
			rayOrigin.y = y;
			rayOrigin.z = 0;

			float tVar1 = rayIntersection(sphere1.coordinates,sphere1.radius,rayOrigin,LookAt);
			float tVar2 = rayIntersection(sphere2.coordinates,sphere2.radius,rayOrigin,LookAt);

			vector normal1 = triangleNormal(shape.vertices[0],shape.vertices[1],shape.vertices[2]);
			vector normal2 = triangleNormal(shape.vertices[0],shape.vertices[1],shape.vertices[3]);
			vector normal3 = triangleNormal(shape.vertices[1],shape.vertices[2],shape.vertices[3]);
			vector normal4 = triangleNormal(shape.vertices[2],shape.vertices[3],shape.vertices[0]);

			float face1 = triangleIntersection(shape.vertices[0],rayOrigin,normal1,LookAt);
			float face2 = triangleIntersection(shape.vertices[0],rayOrigin,normal2,LookAt);
			float face3 = triangleIntersection(shape.vertices[1],rayOrigin,normal3,LookAt);
			float face4 = triangleIntersection(shape.vertices[2],rayOrigin,normal4,LookAt);

			vector vect1;
			vect1.x = (triangleXValue(face1,rayOrigin,LookAt)).x;
			vect1.y = (triangleXValue(face1,rayOrigin,LookAt)).y;
			vect1.z = (triangleXValue(face1,rayOrigin,LookAt)).z;
			vector vect2;
			vect2.x = (triangleXValue(face2,rayOrigin,LookAt)).x;
			vect2.y = (triangleXValue(face2,rayOrigin,LookAt)).y;
			vect2.z = (triangleXValue(face2,rayOrigin,LookAt)).z;
			vector vect3;
			vect3.x = (triangleXValue(face3,rayOrigin,LookAt)).x;
			vect3.y = (triangleXValue(face3,rayOrigin,LookAt)).y;
			vect3.z = (triangleXValue(face3,rayOrigin,LookAt)).z;
			vector vect4;
			vect4.x = (triangleXValue(face4,rayOrigin,LookAt)).x;
			vect4.y = (triangleXValue(face4,rayOrigin,LookAt)).y;
			vect4.z = (triangleXValue(face4,rayOrigin,LookAt)).z;

			bool intersect1 = insideTriangle(shape.vertices[0],shape.vertices[1],shape.vertices[2],vect1,normal1);
			bool intersect2 = insideTriangle(shape.vertices[0],shape.vertices[1],shape.vertices[3],vect2,normal2);
			bool intersect3 = insideTriangle(shape.vertices[1],shape.vertices[2],shape.vertices[3],vect3,normal3);
			bool intersect4 = insideTriangle(shape.vertices[2],shape.vertices[3],shape.vertices[0],vect4,normal4);

			if(intersect1 || intersect2 || intersect3 || intersect4){
				// Figure (c)
				//image(x,y,0) = shape.r;
				//image(x,y,1) = shape.g;
				//image(x,y,2) = shape.b;

				vector normalDir;
				normalDir.x = computeNormal(shape.vertices[0],rayOrigin,LookAt,face1).x;
				normalDir.y = computeNormal(shape.vertices[0],rayOrigin,LookAt,face1).y;
				normalDir.z = computeNormal(shape.vertices[0],rayOrigin,LookAt,face1).z;

				rgb diffuse = diffuseShading(shape.r,shape.g,shape.b,LightDir,normalDir);
				rgb ambient = ambientShading(shape.r,shape.g,shape.b);
				rgb specular = specularShading(sphere1.red,shape.g,shape.b,LightDir,normalDir,rayOrigin);

				image(x,y,0) = (shape.r * ((diffuse.r + ambient.r) + (specular.r)));
				image(x,y,1) = (shape.g * ((diffuse.g + ambient.g) + (specular.g)));
				image(x,y,2) = (shape.b * ((diffuse.b + ambient.b) + (specular.b)));
			}

			if(tVar1 >= 0){
				// Figure (c)
				//image(x,y,0) = sphere1.red;
				//image(x,y,1) = sphere1.green;
				//image(x,y,2) = sphere1.blue;

				vector normalDir;
				normalDir.x = computeNormal(sphere1.coordinates,rayOrigin,LookAt,tVar1).x;
				normalDir.y = computeNormal(sphere1.coordinates,rayOrigin,LookAt,tVar1).y;
				normalDir.z = computeNormal(sphere1.coordinates,rayOrigin,LookAt,tVar1).z;

				rgb diffuse = diffuseShading(sphere1.red,sphere1.green,sphere1.blue,LightDir,normalDir);
				rgb ambient = ambientShading(sphere1.red,sphere1.green,sphere1.blue);
				rgb specular = specularShading(sphere1.red,sphere1.green,sphere1.blue,LightDir,normalDir,rayOrigin);	// lookat -> rayorigin

				image(x,y,0) = (sphere1.red * ((diffuse.r + ambient.r) + (specular.r)));
				image(x,y,1) = (sphere1.green * ((diffuse.g + ambient.g) + (specular.g)));
				image(x,y,2) = (sphere1.blue * ((diffuse.b + ambient.b) + (specular.b)));
			}
			else if(tVar2 >= 0){
				// Figure (c)
				//image(x,y,0) = sphere2.red;
				//image(x,y,1) = sphere2.green;
				//image(x,y,2) = sphere2.blue;

				vector normalDir;
				normalDir.x = computeNormal(sphere2.coordinates,rayOrigin,LookAt,tVar2).x;
				normalDir.y = computeNormal(sphere2.coordinates,rayOrigin,LookAt,tVar2).y;
				normalDir.z = computeNormal(sphere2.coordinates,rayOrigin,LookAt,tVar2).z;

				rgb diffuse = diffuseShading(sphere2.red,sphere2.green,sphere2.blue,LightDir,normalDir);
				rgb ambient = ambientShading(sphere2.red,sphere2.green,sphere2.blue);
				rgb specular = specularShading(sphere2.red,sphere2.green,sphere2.blue,LightDir,normalDir,rayOrigin);	// lookat -> rayorigin

				image(x,y,0) = (sphere2.red * ((diffuse.r + ambient.r) + (specular.r)));
				image(x,y,1) = (sphere2.green * ((diffuse.g + ambient.g) + (specular.g)));
				image(x,y,2) = (sphere2.blue * ((diffuse.b + ambient.b) + (specular.b)));
			}

		}
	}

	CImgDisplay main_disp(image,"Test Image");
	std::cin.ignore().get();
	//system("pause");
	return 0;
}

// L  = Ld + Ls + La
