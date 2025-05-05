// GasSeger.h : Include file for standard system include files,
// or project specific include files.

#pragma once
#include <pcl/filters/box_clipper3D.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/visualization/pcl_visualizer.h>
#include "include/on_nurbs/fitting_curve_pdm.h"
#include "include/on_nurbs/triangulation.h"
#include <omp.h>
#include <atomic>
#include <boost/filesystem.hpp>

#include <vtkBox.h>

#include <filesystem>
#define MPCLIP
using PointT = pcl::PointXYZ;
using PointCT = pcl::PointXYZRGB;
using namespace pcl;
using ViewPTR = visualization::PCLVisualizer::Ptr;

struct Node
{
	Eigen::Vector3f pos;
	Eigen::Vector3f tangent;
	Eigen::Vector3f normal;
	float width;
	float length;
	int generation;
	int total_generation;
	int id;
	int parent_id;
	int type; // 0: stem, 1: leaf
	std::vector<int> children_ids;


	Node()
	{
		pos = Eigen::Vector3f(0, 0, 0);
		tangent = Eigen::Vector3f(0, 0, 0);
		normal = Eigen::Vector3f(0, 0, 0);
		type = -1;
		id = -1;
		parent_id = -1;
		width = 0;
		length = 0;
		generation = 0;
		total_generation = 0;
	}
	Node(
		Eigen::Vector3f in_pos,
		Eigen::Vector3f in_tangent,
		Eigen::Vector3f in_normal,
		float in_width,
		float in_length,
		int in_generation,
		int in_total_generation,
		int in_id,
		int in_parent_id,
		int in_type)
	{
		pos = in_pos;
		tangent = in_tangent;
		normal = in_normal;
		width = in_width;
		length = in_length;
		generation = in_generation;
		total_generation = in_total_generation;
		id = in_id;
		parent_id = in_parent_id;
		type = in_type;
	}

};

class GaussianSeger
{
public:
	ViewPTR viewer;
	ON_NurbsCurve Stem_Curve;
	Eigen::Vector3f Stem_Center;
	bool bShowCurve;
	bool bShowCylinder;
	int LeafNum;
	int ViewportID;
	std::string CurrentCloudId;
	float SampleRadius;

	std::vector<Node> Nodes;
	std::map<std::string, std::vector<Node>> AllNodes;

	GaussianSeger();

	void DumpNodes();


	void SliceParts(
		PointCloud<PointT>::Ptr downsampled_cloud,
		ON_NurbsCurve CurrentCurve,
		Eigen::Vector3f UP,
		std::vector<PointCloud<PointT>::Ptr>& SelectedClouds,
		std::vector<Eigen::Vector3f>& SelectedTangent,
		std::vector<Eigen::Vector3f>& SelectedPos,
		std::vector<Eigen::Vector3f>& SelectedCurvature,
		int& StartIndex,
		float& curvatues,
		float& start_t,
		float& end_t,
		int& start_knot,
		int& AllPartNum);

	void StaticParts(
		std::vector<PointCloud<PointT>::Ptr> SelectedClouds,
		std::vector<Eigen::Vector3f> SelectedTangent,
		std::vector<Eigen::Vector3f> SelectedPos,
		int StartIndex,
		int AllPartNum,
		double& length,
		double& total_torsion,
		double& total_area,
		double& total_width,
		int& CulledNum);

	std::string CalculatePart(
		pcl::PointCloud<PointT>::Ptr cloud,
		std::string cloud_id,
		bool bStem);

	std::string CalculateStem(
		PointCloud<PointT>::Ptr cloud,
		std::string cloud_id,
		bool ShouldDownSample);
	PointCloud<PointT>::Ptr DownSample(
		PointCloud<PointT>::Ptr cloud,
		const float leaf_size = 0.01);
	ON_NurbsCurve SimpleFitCurve(
		PointCloud<PointT>::Ptr cloud,
		unsigned n_control_points,
		on_nurbs::vector_vec3d RefPoints);


	template <typename T>
	void UpdateSphere(
		T cur_point,
		double radius = 0.1,
		std::string name = "LeafIntersection");
	void LaunchViewer();

	void AddCloud(
		const std::string& Name,
		PointCloud<PointT>::Ptr cloud)
	{
		clouds.emplace(Name, cloud);
	}

	void InitCloudColor()
	{
		cloudColors.clear();
		for (std::pair<const std::string, shared_ptr<PointCloud<PointXYZ>>> pair : clouds)
		{
			float r = (rand() % 65536) / 140000.f + 0.2f;
			float g = (rand() % 65536) / 140000.f + 0.2f;
			float b = (rand() % 65536) / 140000.f + 0.2f;
			cloudColors.emplace(pair.first, Eigen::Vector3f({ r, g, b }));
		}
	}

	void ShowClouds();
	void InitializeViewer();

private:
	std::map<std::string, PointCloud<PointT>::Ptr> clouds;
	std::map<std::string, Eigen::Vector3f> cloudColors;
};

template <typename T>
void GaussianSeger::UpdateSphere(
	T cur_point,
	double radius,
	std::string name)
{
	if (!bShowCurve)
		return;
	if (!viewer->updateSphere(cur_point, radius, 1, 0, 0, name))
		viewer->addSphere(cur_point, radius, 1, 0, 0, name);
	//viewer->spinOnce(100);
}

class SphereClipper
{
public:
	double radius;
	Eigen::Vector3f pos;

	SphereClipper()
	{
	};

	SphereClipper(
		double in_radius,
		Eigen::Vector3f in_pos)
		: radius(in_radius),
		pos(in_pos)
	{
	};

	virtual bool clipPoint3D(
		const PointT point)
	{
		auto p = pos - point.getVector3fMap();
		return (p.squaredNorm()) < radius * radius;
	}

	void Clip3DPoints(
		PointCloud<PointT>::Ptr& cloud,
		Indices& indices)
	{
		indices.clear();
		PointCloud<PointT>& cloud_in = *cloud;
		if (indices.empty())
		{
			indices.reserve(cloud->size());
#ifdef MPCLIP
			std::atomic<int> idx(0);
			int num_t = cloud->size() > 1000 ? 8 : 2;
			indices.resize(cloud->size());
#pragma omp parallel for num_threads(num_t)
			for (std::size_t pIdx = 0; pIdx < cloud->size(); ++pIdx)
				if (clipPoint3D(cloud_in[pIdx]))
					indices[idx.fetch_add(1)] = pIdx;

			indices.resize(idx.load());
#elif
			for (std::size_t pIdx = 0; pIdx < cloud->size(); ++pIdx)
				if (clipPoint3D(cloud_in[pIdx]))
					indices.push_back(pIdx);
#endif
		}
	}
};

class CylinderClipper : public SphereClipper
{
public:
	double range;
	Eigen::Vector3f dir;

	CylinderClipper()
	{
	};

	CylinderClipper(
		double in_range,
		double in_radius,
		Eigen::Vector3f in_pos,
		Eigen::Vector3f in_dir)
		: SphereClipper(in_radius, in_pos),
		range(in_range),
		dir(in_dir)
	{
	};

	virtual bool clipPoint3D(
		const PointT point) override
	{
		auto p = pos - point.getVector3fMap();
		float local_z = dir.dot(p);

		return (p.squaredNorm() - local_z * local_z) < radius * radius && abs(local_z) < range;
	}
};


ON_NurbsCurve initNurbsCurvePCAOffset(
	int order,
	const on_nurbs::vector_vec3d& data,
	int ncps,
	double rf = 1.0);

void PointCloud2Vector2d(
	PointCloud<PointT>::Ptr cloud,
	on_nurbs::vector_vec3d& data);

void VisualizeCurve(
	ViewPTR viewer,
	std::string cloud_id,
	PointCloud<PointCT>::Ptr cloud,
	double r,
	double g,
	double b,
	bool show_cps);


void VisualTagent(
	ViewPTR viewer,
	std::string cloud_id,
	int id,
	ON_3dPoint& Point,
	ON_3dVector& Vector);

void IterSample(
	PointCloud<PointT>::Ptr ClonedCloud,
	on_nurbs::vector_vec3d& OutPoints,
	GaussianSeger* viewer,
	float Radius = 0.01);
