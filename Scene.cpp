#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>

#include "Scene.h"
#include "Matrix4.h"
#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"

using namespace tinyxml2;
using namespace std;

/*
	Transformations, clipping, culling, rasterization are done here.
	You may define helper functions.
*/

Matrix4 camera_transformation(Camera *camera)
{
	// Function to move camera to the origin
	Matrix4 M_cam;
	M_cam.val[0][0] = camera->u.x;
	M_cam.val[0][1] = camera->u.y;
	M_cam.val[0][2] = camera->u.z;
	M_cam.val[0][3] = -1 * (camera->u.x * camera->pos.x + camera->u.y * camera->pos.y + camera->u.z * camera->pos.z);

	M_cam.val[1][0] = camera->v.x;
	M_cam.val[1][1] = camera->v.y;
	M_cam.val[1][2] = camera->v.z;
	M_cam.val[1][3] = -1 * (camera->v.x * camera->pos.x + camera->v.y * camera->pos.y + camera->v.z * camera->pos.z);

	M_cam.val[2][0] = camera->w.x;
	M_cam.val[2][1] = camera->w.y;
	M_cam.val[2][2] = camera->w.z;
	M_cam.val[2][3] = -1 * (camera->w.x * camera->pos.x + camera->w.y * camera->pos.y + camera->w.z * camera->pos.z);

	M_cam.val[3][0] = 0.0;
	M_cam.val[3][1] = 0.0;
	M_cam.val[3][2] = 0.0;
	M_cam.val[3][3] = 1.0;
	return M_cam;
}

Matrix4 orthographic_projection(Camera *camera)
{
	Matrix4 M_ortho;
	M_ortho.val[0][0] = (2.0) / (camera->right - camera->left);
	M_ortho.val[0][1] = 0.0;
	M_ortho.val[0][2] = 0.0;
	M_ortho.val[0][3] = (-1.0) * (camera->right + camera->left) / (camera->right - camera->left);

	M_ortho.val[1][0] = 0.0;
	M_ortho.val[1][1] = (2.0) / (camera->top - camera->bottom);
	M_ortho.val[1][2] = 0.0;
	M_ortho.val[1][3] = (-1.0) * (camera->top + camera->bottom) / (camera->top - camera->bottom);

	M_ortho.val[2][0] = 0.0;
	M_ortho.val[2][1] = 0.0;
	M_ortho.val[2][2] = (-2.0) / (camera->far - camera->near);
	M_ortho.val[2][3] = (-1.0) * (camera->far + camera->near) / (camera->far - camera->near);

	M_ortho.val[3][0] = 0.0;
	M_ortho.val[3][1] = 0.0;
	M_ortho.val[3][2] = 0.0;
	M_ortho.val[3][3] = 1.0;
	return M_ortho;
}

Matrix4 perspective_projection(Camera *camera)
{
	Matrix4 M_per;
	// Matrix4 M_ortho;
	// Matrix4 M_p2o;
	// M_p2o.val[0][0] = camera.near;
	// M_p2o.val[1][1] = camera.near;
	// M_p2o.val[2][2] = camera.far + camera.near;
	// M_p2o.val[2][3] = camera.far * camera.near;
	// M_p2o.val[3][2] = -1.0;

	// Matrix4 result = multiplyMatrixWithMatrix(M_ortho, M_p2o);
	// return result;
	M_per.val[0][0] = (2.0) * camera->near / (camera->right - camera->left);
	M_per.val[0][1] = 0.0;
	M_per.val[0][2] = (camera->right + camera->left) / (camera->right - camera->left);
	M_per.val[0][3] = 0.0;

	M_per.val[1][0] = 0.0;
	M_per.val[1][1] = (2.0) * camera->near / (camera->top - camera->bottom);
	M_per.val[1][2] = (camera->top + camera->bottom) / (camera->top - camera->bottom);
	M_per.val[1][3] = 0.0;

	M_per.val[2][0] = 0.0;
	M_per.val[2][1] = 0.0;
	M_per.val[2][2] = (-1.0) * (camera->far + camera->near) / (camera->far - camera->near);
	M_per.val[2][3] = (-2.0) * camera->far * camera->near / (camera->far - camera->near);

	M_per.val[3][0] = 0.0;
	M_per.val[3][1] = 0.0;
	M_per.val[3][2] = -1.0;
	M_per.val[3][3] = 0.0;
	return M_per;
}

void viewport_transformation(Camera *camera, double M_vp[3][4])
{
	M_vp[0][0] = camera->horRes / 2.0;
	M_vp[0][1] = 0.0;
	M_vp[0][2] = 0.0;
	M_vp[0][3] = (camera->horRes - 1.0) / 2.0;

	M_vp[1][0] = 0.0;
	M_vp[1][1] = camera->verRes / 2.0;
	M_vp[1][2] = 0.0;
	M_vp[1][3] = (camera->verRes - 1.0) / 2.0;

	M_vp[2][0] = 0.0;
	M_vp[2][1] = 0.0;
	M_vp[2][2] = 0.5;
	M_vp[2][3] = 0.5;
}

Matrix4 scaling_transformation(Scaling scaling)
{

	Matrix4 scaling_matrix;
	scaling_matrix.val[0][0] = scaling.sx;
	scaling_matrix.val[0][1] = 0.0;
	scaling_matrix.val[0][2] = 0.0;
	scaling_matrix.val[0][3] = 0.0;

	scaling_matrix.val[1][0] = 0.0;
	scaling_matrix.val[1][1] = scaling.sy;
	scaling_matrix.val[1][2] = 0.0;
	scaling_matrix.val[1][3] = 0.0;

	scaling_matrix.val[2][0] = 0.0;
	scaling_matrix.val[2][1] = 0.0;
	scaling_matrix.val[2][2] = scaling.sz;
	scaling_matrix.val[2][3] = 0.0;

	scaling_matrix.val[3][0] = 0.0;
	scaling_matrix.val[3][1] = 0.0;
	scaling_matrix.val[3][2] = 0.0;
	scaling_matrix.val[3][3] = 1.0;

	return scaling_matrix;
}

Matrix4 translation_transformation(Translation translation)
{
	Matrix4 translation_matrix;
	translation_matrix.val[0][0] = 1.0;
	translation_matrix.val[0][1] = 0.0;
	translation_matrix.val[0][2] = 0.0;
	translation_matrix.val[0][3] = translation.tx;

	translation_matrix.val[1][0] = 0.0;
	translation_matrix.val[1][1] = 1.0;
	translation_matrix.val[1][2] = 0.0;
	translation_matrix.val[1][3] = translation.ty;

	translation_matrix.val[2][0] = 0.0;
	translation_matrix.val[2][1] = 0.0;
	translation_matrix.val[2][2] = 1.0;
	translation_matrix.val[2][3] = translation.tz;

	translation_matrix.val[3][0] = 0.0;
	translation_matrix.val[3][1] = 0.0;
	translation_matrix.val[3][2] = 0.0;
	translation_matrix.val[3][3] = 1.0;

	return translation_matrix;
}

Matrix4 rotation_transformation(Rotation rotation)
{
	Vec3 u;
	u.x = rotation.ux;
	u.y = rotation.uy;
	u.z = rotation.uz;
	// assign edilmiyor olabilir.
	u = normalizeVec3(u);

	Vec3 v;
	double min_val = min(min(abs(u.x), abs(u.y)), abs(u.z));
	if (min_val == abs(u.x))
	{
		v.x = 0.0;
		v.y = (-1.0) * u.z;
		v.z = u.y;
	}
	else if (min_val == abs(u.y))
	{
		v.y = 0.0;
		v.x = (-1.0) * u.z;
		v.z = u.x;
	}
	else if (min_val == abs(u.z))
	{
		v.z = 0.0;
		v.x = (-1.0) * u.y;
		v.y = u.x;
	}
	// assign edilmiyor olabilir.
	Vec3 w;
	w = crossProductVec3(u, v);
	v = normalizeVec3(v);
	w = normalizeVec3(w);

	Matrix4 M_inverse;
	M_inverse.val[0][0] = u.x;
	M_inverse.val[0][1] = v.x;
	M_inverse.val[0][2] = w.x;
	M_inverse.val[1][0] = u.y;
	M_inverse.val[1][1] = v.y;
	M_inverse.val[1][2] = w.y;
	M_inverse.val[2][0] = u.z;
	M_inverse.val[2][1] = v.z;
	M_inverse.val[2][2] = w.z;
	M_inverse.val[3][3] = 1.0;

	Matrix4 M;
	M.val[0][0] = u.x;
	M.val[0][1] = u.y;
	M.val[0][2] = u.z;
	M.val[0][3] = 0.0;
	M.val[1][0] = v.x;
	M.val[1][1] = v.y;
	M.val[1][2] = v.z;
	M.val[1][3] = 0.0;
	M.val[2][0] = w.x;
	M.val[2][1] = w.y;
	M.val[2][2] = w.z;
	M.val[2][3] = 0.0;
	M.val[3][3] = 1.0;

	double theta = (rotation.angle * M_PI) / 180.0;
	Matrix4 M_Rx;
	M_Rx.val[0][0] = 1.0;
	M_Rx.val[1][1] = cos(theta);
	M_Rx.val[1][2] = (-1.0) * sin(theta);
	M_Rx.val[2][1] = sin(theta);
	M_Rx.val[2][2] = cos(theta);
	M_Rx.val[3][3] = 1.0;

	Matrix4 M_mid;
	M_mid = multiplyMatrixWithMatrix(M_Rx, M);
	Matrix4 rotation_matrix;
	rotation_matrix = multiplyMatrixWithMatrix(M_inverse, M_mid);
	return rotation_matrix;
}

Matrix4 Scene::modeling_transformation(Mesh *mesh)
{
	Matrix4 M_result = getIdentityMatrix();
	for (int i = 0; i < mesh->numberOfTransformations; i++)
	{
		int transformation_id = mesh->transformationIds[i];
		char transformation_type = mesh->transformationTypes[i];
		if (transformation_type == 'r')
		{
			Matrix4 M_r = rotation_transformation(*rotations[transformation_id - 1]);
			M_result = multiplyMatrixWithMatrix(M_r, M_result);
		}
		else if (transformation_type == 's')
		{
			Matrix4 M_s = scaling_transformation(*scalings[transformation_id - 1]);
			M_result = multiplyMatrixWithMatrix(M_s, M_result);
		}
		else if (transformation_type == 't')
		{
			Matrix4 M_t = translation_transformation(*translations[transformation_id - 1]);
			M_result = multiplyMatrixWithMatrix(M_t, M_result);
		}
	}
	return M_result;
}

void Scene::midpoint_algorithm(Vec3 v0, Vec3 v1)
{
	Color color;
	Color dc;
	int step;

	Vec3 v0_new;
	v0_new.x = v0.x;
	v0_new.y = v0.y;
	v0_new.colorId = v0.colorId;
	Vec3 v1_new;
	v1_new.x = v1.x;
	v1_new.y = v1.y;
	v1_new.colorId = v1.colorId;

	bool slope_control = abs(v1.y - v0.y) > abs(v1.x - v0.x);
	if (slope_control)
	{
		v0_new.x = v0.y;
		v0_new.y = v0.x;
		v1_new.x = v1.y;
		v1_new.y = v1.x;
	}
	if (v1_new.x < v0_new.x)
	{
		double temp_x = v0_new.x;
		double temp_y = v0_new.y;
		v0_new.x = v1_new.x;
		v0_new.y = v1_new.y;
		v1_new.x = temp_x;
		v1_new.y = temp_y;
		v1_new.colorId = v0.colorId;
		v0_new.colorId = v1.colorId;
	}
	if (v0_new.y < v1_new.y)
	{
		step = 1;
	}
	else
	{
		step = -1;
	}

	int x1_minus_x0 = v1_new.x - v0_new.x;
	int y1_minus_y0 = abs(v1_new.y - v0_new.y);
	int y0_minus_y1 = (-1) * y1_minus_y0;
	int y = v0_new.y;
	float d = y0_minus_y1 + 0.5 * (x1_minus_x0);

	// c=c0
	color.r = colorsOfVertices[v0_new.colorId - 1]->r;
	color.b = colorsOfVertices[v0_new.colorId - 1]->b;
	color.g = colorsOfVertices[v0_new.colorId - 1]->g;

	// dc = (c1-c0) / (x1-x0)
	dc.r = (colorsOfVertices[v1_new.colorId - 1]->r - color.r) / x1_minus_x0;
	dc.b = (colorsOfVertices[v1_new.colorId - 1]->b - color.b) / x1_minus_x0;
	dc.g = (colorsOfVertices[v1_new.colorId - 1]->g - color.g) / x1_minus_x0;

	for (int x = (int)v0_new.x; x <= (int)v1_new.x; x++)
	{
		// round(c)
		color.r = makeBetweenZeroAnd255(color.r);
		color.b = makeBetweenZeroAnd255(color.b);
		color.g = makeBetweenZeroAnd255(color.g);
		if (slope_control)
		{
			image[y][x] = color;
		}
		else
		{
			image[x][y] = color;
		}
		if (d < 0)
		{
			y = y + step;
			d += (y0_minus_y1 + x1_minus_x0);
		}
		else
		{
			d += y0_minus_y1;
		}
		color.r += dc.r;
		color.b += dc.b;
		color.g += dc.g;
	}
}

bool is_visible(double den, double num, double &te, double &tl)
{
	double t;
	if (den > 0)
	{
		t = num / den;
		if (t > tl)
			return false;
		if (t > te)
			te = t;
	}
	else if (den < 0)
	{
		t = num / den;
		if (t < te)
			return false;
		if (t < tl)
			tl = t;
	}
	else if (num > 0)
	{
		return false;
	}
	return true;
}

void Scene::line_rasterization(Vec3 v0, Vec3 v1, Camera *camera)
{

	double te = 0;
	double tl = 1;
	bool visible = false;
	int xmin = 0;
	int ymin = 0;
	int xmax = camera->horRes - 1;
	int ymax = camera->verRes - 1;
	v0.x = (int)v0.x;
	v1.x = (int)v1.x;
	v0.y = (int)v0.y;
	v1.y = (int)v1.y;

	int x0 = v0.x;
	int y0 = v0.y;
	int x1 = v1.x;
	int y1 = v1.y;

	double dx = x1 - x0;
	double dy = y1 - y0;

	if (is_visible(dx, (xmin - x0), te, tl))
	{
		if (is_visible(-1 * dx, (x0 - xmax), te, tl))
		{
			if (is_visible(dy, (ymin - y0), te, tl))
			{
				if (is_visible(-1 * dy, (y0 - ymax), te, tl))
				{
					visible = true;
					if (tl < 1)
					{
						v1.x = v0.x + dx * tl;
						v1.y = v0.y + dy * tl;
					}
					if (te > 0)
					{
						v0.x = v0.x + dx * te;
						v0.y = v0.y + dy * te;
					}
					midpoint_algorithm(v0, v1);
				}
			}
		}
	}
}

void Scene::triangle_rasterization(Vec3 v1, Vec3 v2, Vec3 v3, Camera *camera)
{
	int x_min = min_3(v1.x, v2.x, v3.x);
	int y_min = min_3(v1.y, v2.y, v3.y);
	int x_max = max_3(v1.x, v2.x, v3.x);
	int y_max = max_3(v1.y, v2.y, v3.y);

	int x0 = v1.x;
	int y0 = v1.y;
	int x1 = v2.x;
	int y1 = v2.y;
	int x2 = v3.x;
	int y2 = v3.y;
	Color lastColor;
	Color *color0 = colorsOfVertices[v1.colorId - 1];
	Color *color1 = colorsOfVertices[v2.colorId - 1];
	Color *color2 = colorsOfVertices[v3.colorId - 1];

	for (int y = y_min; y <= y_max; y++)
	{
		for (int x = x_min; x <= x_max; x++)
		{
			if (x < camera->horRes && x >= 0 && y < camera->verRes && y >= 0)
			{
				double alpha = line_equation(x2, y2, x1, y1, x, y) / line_equation(x2, y2, x1, y1, x0, y0);
				double beta = line_equation(x0, y0, x2, y2, x, y) / line_equation(x0, y0, x2, y2, x1, y1);
				double gamma = line_equation(x1, y1, x0, y0, x, y) / line_equation(x1, y1, x0, y0, x2, y2);

				if (alpha >= 0 && beta >= 0 && gamma >= 0)
				{
					lastColor.r = makeBetweenZeroAnd255(color0->r * alpha + color1->r * beta + color2->r * gamma);
					lastColor.b = makeBetweenZeroAnd255(color0->b * alpha + color1->b * beta + color2->b * gamma);
					lastColor.g = makeBetweenZeroAnd255(color0->g * alpha + color1->g * beta + color2->g * gamma);
					image[x][y] = lastColor;
				}
			}
		}
	}
}

void Scene::forwardRenderingPipeline(Camera *camera)
{
	// It is enough to calculate these matrix for every camera not every mesh:

	// camera transformation
	Matrix4 M_cam = camera_transformation(camera);

	// viewport transformation
	double M_vp[3][4];
	viewport_transformation(camera, M_vp);

	// projection transformation
	Matrix4 M_proj;
	// need to check for projection type
	// 1 for perspective, 0 for orthographic
	if (camera->projectionType)
	{
		M_proj = perspective_projection(camera);
	}
	else
	{
		M_proj = orthographic_projection(camera);
	}

	// create viewing transformation matrix without viewport transformation:
	Matrix4 M_proj_and_cam = multiplyMatrixWithMatrix(M_proj, M_cam);

	for (int i = 0; i < meshes.size(); i++)
	{
		Mesh *mesh = (meshes[i]);

		// calculate modeling transformation matrix for each mesh
		Matrix4 M_model = modeling_transformation(mesh);

		// create matrix for modeling, camera and projection transformation
		Matrix4 M_proj_cam_model = multiplyMatrixWithMatrix(M_proj_and_cam, M_model);

		for (int j = 0; j < mesh->numberOfTriangles; j++)
		{
			// each vertex of a triangle
			Vec4 v1;
			Vec4 v2;
			Vec4 v3;
			v1.x = vertices[mesh->triangles[j].vertexIds[0] - 1]->x;
			v1.y = vertices[mesh->triangles[j].vertexIds[0] - 1]->y;
			v1.z = vertices[mesh->triangles[j].vertexIds[0] - 1]->z;
			v1.t = 1.0;
			v1.colorId = vertices[mesh->triangles[j].vertexIds[0] - 1]->colorId;

			v2.x = vertices[mesh->triangles[j].vertexIds[1] - 1]->x;
			v2.y = vertices[mesh->triangles[j].vertexIds[1] - 1]->y;
			v2.z = vertices[mesh->triangles[j].vertexIds[1] - 1]->z;
			v2.t = 1.0;
			v2.colorId = vertices[mesh->triangles[j].vertexIds[1] - 1]->colorId;

			v3.x = vertices[mesh->triangles[j].vertexIds[2] - 1]->x;
			v3.y = vertices[mesh->triangles[j].vertexIds[2] - 1]->y;
			v3.z = vertices[mesh->triangles[j].vertexIds[2] - 1]->z;
			v3.t = 1.0;
			v3.colorId = vertices[mesh->triangles[j].vertexIds[2] - 1]->colorId;

			// apply modeling, camera and projection transformation to each vertex:
			Vec4 v1_transformed = multiplyMatrixWithVec4(M_proj_cam_model, v1);
			Vec4 v2_transformed = multiplyMatrixWithVec4(M_proj_cam_model, v2);
			Vec4 v3_transformed = multiplyMatrixWithVec4(M_proj_cam_model, v3);

			if (v1_transformed.t != 0)
			{
				v1_transformed.x = v1_transformed.x / v1_transformed.t;
				v1_transformed.y = v1_transformed.y / v1_transformed.t;
				v1_transformed.z = v1_transformed.z / v1_transformed.t;
				v1_transformed.t = 1;
			}
			if (v2_transformed.t != 0)
			{
				v2_transformed.x = v2_transformed.x / v2_transformed.t;
				v2_transformed.y = v2_transformed.y / v2_transformed.t;
				v2_transformed.z = v2_transformed.z / v2_transformed.t;
				v2_transformed.t = 1;
			}
			if (v3_transformed.t != 0)
			{
				v3_transformed.x = v3_transformed.x / v3_transformed.t;
				v3_transformed.y = v3_transformed.y / v3_transformed.t;
				v3_transformed.z = v3_transformed.z / v3_transformed.t;
				v3_transformed.t = 1;
			}

			Vec3 v1_transformed_vp;
			Vec3 v2_transformed_vp;
			Vec3 v3_transformed_vp;

			// apply viewport transformation to each vertex:
			multiply_3x4_MatrixWithVec4(v1_transformed_vp, M_vp, v1_transformed);
			multiply_3x4_MatrixWithVec4(v2_transformed_vp, M_vp, v2_transformed);
			multiply_3x4_MatrixWithVec4(v3_transformed_vp, M_vp, v3_transformed);

			if (cullingEnabled)
			{
				Vec3 edge1 = subtractVec3(v2_transformed_vp, v1_transformed_vp);
				Vec3 edge2 = subtractVec3(v3_transformed_vp, v1_transformed_vp);
				Vec3 normal = crossProductVec3(edge2, edge1);
				normal = normalizeVec3(normal);
				Vec3 viewing_vector = v1_transformed_vp;
				viewing_vector = normalizeVec3(viewing_vector);
				if (dotProductVec3(viewing_vector, normal) > 0)
				{
					continue;
				}
			}

			v1_transformed_vp.colorId = v1.colorId;
			v2_transformed_vp.colorId = v2.colorId;
			v3_transformed_vp.colorId = v3.colorId;

			// Solid
			if (mesh->type)
			{
				triangle_rasterization(v1_transformed_vp, v2_transformed_vp, v3_transformed_vp, camera);
			}
			// Wireframe
			else
			{
				line_rasterization(v1_transformed_vp, v2_transformed_vp, camera);
				line_rasterization(v2_transformed_vp, v3_transformed_vp, camera);
				line_rasterization(v3_transformed_vp, v1_transformed_vp, camera);
			}
		}
	}
}

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *pElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *pRoot = xmlDoc.FirstChild();

	// read background color
	pElement = pRoot->FirstChildElement("BackgroundColor");
	str = pElement->GetText();
	sscanf(str, "%lf %lf %lf", &(backgroundColor.r), &(backgroundColor.g), &(backgroundColor.b));

	// read culling
	pElement = pRoot->FirstChildElement("Culling");
	if (pElement != NULL)
	{
		str = pElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			cullingEnabled = true;
		}
		else
		{
			cullingEnabled = false;
		}
	}

	// read cameras
	pElement = pRoot->FirstChildElement("Cameras");
	XMLElement *pCamera = pElement->FirstChildElement("Camera");
	XMLElement *camElement;
	while (pCamera != NULL)
	{
		Camera *cam = new Camera();

		pCamera->QueryIntAttribute("id", &cam->cameraId);

		// read projection type
		str = pCamera->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			cam->projectionType = 0;
		}
		else
		{
			cam->projectionType = 1;
		}

		camElement = pCamera->FirstChildElement("Position");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &(cam->pos.x), &(cam->pos.y), &(cam->pos.z));

		camElement = pCamera->FirstChildElement("Gaze");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &(cam->gaze.x), &(cam->gaze.y), &(cam->gaze.z));

		camElement = pCamera->FirstChildElement("Up");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &(cam->v.x), &(cam->v.y), &(cam->v.z));

		cam->gaze = normalizeVec3(cam->gaze);
		cam->u = crossProductVec3(cam->gaze, cam->v);
		cam->u = normalizeVec3(cam->u);

		cam->w = inverseVec3(cam->gaze);
		cam->v = crossProductVec3(cam->u, cam->gaze);
		cam->v = normalizeVec3(cam->v);

		camElement = pCamera->FirstChildElement("ImagePlane");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &cam->left, &cam->right, &cam->bottom, &cam->top,
			   &cam->near, &cam->far, &cam->horRes, &cam->verRes);

		camElement = pCamera->FirstChildElement("OutputName");
		str = camElement->GetText();
		cam->outputFileName = string(str);

		cameras.push_back(cam);

		pCamera = pCamera->NextSiblingElement("Camera");
	}

	// read vertices
	pElement = pRoot->FirstChildElement("Vertices");
	XMLElement *pVertex = pElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (pVertex != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = pVertex->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = pVertex->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		vertices.push_back(vertex);
		colorsOfVertices.push_back(color);

		pVertex = pVertex->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	pElement = pRoot->FirstChildElement("Translations");
	XMLElement *pTranslation = pElement->FirstChildElement("Translation");
	while (pTranslation != NULL)
	{
		Translation *translation = new Translation();

		pTranslation->QueryIntAttribute("id", &translation->translationId);

		str = pTranslation->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		translations.push_back(translation);

		pTranslation = pTranslation->NextSiblingElement("Translation");
	}

	// read scalings
	pElement = pRoot->FirstChildElement("Scalings");
	XMLElement *pScaling = pElement->FirstChildElement("Scaling");
	while (pScaling != NULL)
	{
		Scaling *scaling = new Scaling();

		pScaling->QueryIntAttribute("id", &scaling->scalingId);
		str = pScaling->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		scalings.push_back(scaling);

		pScaling = pScaling->NextSiblingElement("Scaling");
	}

	// read rotations
	pElement = pRoot->FirstChildElement("Rotations");
	XMLElement *pRotation = pElement->FirstChildElement("Rotation");
	while (pRotation != NULL)
	{
		Rotation *rotation = new Rotation();

		pRotation->QueryIntAttribute("id", &rotation->rotationId);
		str = pRotation->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		rotations.push_back(rotation);

		pRotation = pRotation->NextSiblingElement("Rotation");
	}

	// read meshes
	pElement = pRoot->FirstChildElement("Meshes");

	XMLElement *pMesh = pElement->FirstChildElement("Mesh");
	XMLElement *meshElement;
	while (pMesh != NULL)
	{
		Mesh *mesh = new Mesh();

		pMesh->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = pMesh->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = 0;
		}
		else
		{
			mesh->type = 1;
		}

		// read mesh transformations
		XMLElement *pTransformations = pMesh->FirstChildElement("Transformations");
		XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

		while (pTransformation != NULL)
		{
			char transformationType;
			int transformationId;

			str = pTransformation->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			pTransformation = pTransformation->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *clone_str;
		int v1, v2, v3;
		XMLElement *pFaces = pMesh->FirstChildElement("Faces");
		str = pFaces->GetText();
		clone_str = strdup(str);

		row = strtok(clone_str, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		meshes.push_back(mesh);

		pMesh = pMesh->NextSiblingElement("Mesh");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFileName.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFileName << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}