#include "GasSeger.h"

#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/io/obj_io.h>  
#include <pcl/io/ply_io.h>  
#include <pcl/PolygonMesh.h>
#include <pcl/conversions.h>  
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/features/don.h>
#include <string>

#include <pcl/features/normal_3d_omp.h>
#include <pcl/filters/conditional_removal.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/statistical_outlier_removal.h>

#include <pcl/gpu/segmentation/gpu_extract_clusters.h>
#include <pcl/gpu/octree/octree.hpp>
#include <pcl/segmentation/sac_segmentation.h>
// 可视化mesh模型和点云

#define USE_GPU 1
namespace fs = std::filesystem;

using namespace pcl;
#if USE_GPU
using PointT = pcl::PointXYZ;
#else
using PointT = pcl::PointXYZRGBNormal;
#endif
using PointNT = pcl::Normal;
using PointOutT = pcl::PointNormal;
using SearchPtr = pcl::search::Search<PointT>::Ptr;

template<typename ExtractPointType>
typename pcl::PointCloud<ExtractPointType>::Ptr Extract(typename pcl::PointCloud<ExtractPointType>::Ptr cloud, pcl::PointIndices::Ptr inliers_plane, bool Negtive = true)
{
	ExtractIndices<ExtractPointType > extract;
	extract.setInputCloud(cloud);
	extract.setIndices(inliers_plane);
	extract.setNegative(Negtive);
	typename pcl::PointCloud<ExtractPointType>::Ptr new_cloud(new pcl::PointCloud<ExtractPointType>());
	extract.filter(*new_cloud);
	return new_cloud;
}

void SegPart(pcl::PointCloud<PointT>::Ptr& cloud, pcl::PointCloud<PointNT>::Ptr& cloud_normals)
{
	cout << "Input cloud: " << cloud->points.size() << " data points." << endl;
	SACSegmentationFromNormals<PointT, PointNT> seg;
	pcl::ModelCoefficients::Ptr coefficients_plane(new pcl::ModelCoefficients);
	pcl::PointIndices::Ptr inliers_plane(new pcl::PointIndices);
	// ------------------------点云分割，提取平面上的点--------------------------
	seg.setOptimizeCoefficients(true);
	seg.setModelType(pcl::SACMODEL_PLANE);
	seg.setNormalDistanceWeight(0.2);
	seg.setMethodType(pcl::SAC_RANSAC);
	seg.setMaxIterations(100);
	seg.setDistanceThreshold(0.3);
	seg.setInputCloud(cloud);
	seg.setInputNormals(cloud_normals);
	seg.segment(*inliers_plane, *coefficients_plane);//获取平面模型系数和平面上的点
	cout << "Plane coefficients: " << *coefficients_plane << endl;
	//----------------------------------提取平面以外的---------------------------------
	cloud = Extract<PointT>(cloud, inliers_plane);
	cloud_normals = Extract<PointNT>(cloud_normals, inliers_plane);
	cout << "PointCloud representing the planar component: " << cloud->points.size() << " data points." << endl;
	return;

}
void meshPointCloudandViewer(pcl::PointCloud<PointT>::Ptr& cloud)
{
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("3D PointCloud Viewer"));
	int v2;
	viewer->createViewPort(0.0, 0.0, 1.0, 1.0, v2);  // 左侧窗口
	viewer->setBackgroundColor(0.0, 0.0, 0.0, v2);  // 黑色背景
	viewer->addText("Original PointCloud", 10, 10, "vp1_text", v2);  // 标题
	pcl::visualization::PointCloudColorHandlerRGBField<PointT> cloud_color_handler(cloud);  // 绿色
	viewer->addPointCloud(cloud, "original_cloud", v2);
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.5, 0.8, .5, "original_cloud");

	// 添加坐标系
   /* viewer->addCoordinateSystem(0.1);
	viewer->initCameraParameters();*/

	// 可视化循环
	while (!viewer->wasStopped())
	{
		viewer->spinOnce(100);
	}

}

template<typename PointType>
void sor_filter(typename pcl::PointCloud<PointType>::Ptr& src, int num)
{
	//滤波离群点
	pcl::StatisticalOutlierRemoval<PointType> sor;
	sor.setInputCloud(src);
	sor.setMeanK(num);
	sor.setStddevMulThresh(1.0);
	sor.filter(*src);
}

void SaveCluster(std::vector<pcl::PointIndices> cluster_indices, pcl::PointCloud<PointT>::Ptr cloud, std::string outfile, std::string prefix = "cluster")
{
	std::cout << "No of clusters formed are " << cluster_indices.size() << std::endl;
	pcl::PLYWriter writer;

	// Saving the clusters in separate pcd files
	int j = 0;
	for (const auto& cluster : cluster_indices)
	{
		pcl::PointCloud<PointT>::Ptr cloud_cluster(new pcl::PointCloud<PointT>);
		for (const auto& index : cluster.indices) {
			cloud_cluster->push_back((*cloud)[index]);
		}
		cloud_cluster->width = cloud_cluster->size();
		cloud_cluster->height = 1;
		cloud_cluster->is_dense = true;

		std::cout << "PointCloud representing the Cluster using xyzn: " << cloud_cluster->size() << " data " << std::endl;
		std::stringstream ss;
		ss << "./" << outfile << "/" << prefix << "_" << j << ".ply";
		writer.write<PointT>(ss.str(), *cloud_cluster, false);
		++j;
	}
}

void Euclidean(pcl::PointCloud<PointT>::Ptr cloud, float tolerance, unsigned int min_cluster_size, unsigned int max_cluster_size, std::vector<pcl::PointIndices>& cluster_indices)
{
	cout << "Euclidean cluster extraction :: " << cloud->points.size() << " data points." << endl;
	pcl::gpu::Octree::PointCloud cloud_device;
	cloud_device.upload(cloud->points);

	pcl::gpu::Octree::Ptr octree_device(new pcl::gpu::Octree);
	octree_device->setCloud(cloud_device);
	octree_device->build();

	pcl::gpu::EuclideanClusterExtraction<pcl::PointXYZ> gec;
	gec.setClusterTolerance(tolerance); // 2cm
	gec.setMinClusterSize(min_cluster_size);
	gec.setMaxClusterSize(max_cluster_size);
	gec.setSearchMethod(octree_device);
	gec.setHostCloud(cloud);
	gec.extract(cluster_indices);
}

void DisplayCloud(bool ShouldDisplay, pcl::PointCloud<PointT>::Ptr cloud)
{
	if (ShouldDisplay)
		meshPointCloudandViewer(cloud);
}

int main(int argc, char* argv[])
{
	if (argc < 2) {
		std::cerr << "Expected 2 arguments: inputfile outputfile" << std::endl;
	}

	///The file to read from.
	std::string infile = argv[1];

	///The file to output to.
	std::string outfile = argv[2];
	// Extracting Euclidean clusters using cloud and its normals
	std::vector<pcl::PointIndices> cluster_indices;
	float tolerance = 0.01f;
	constexpr double eps_angle = 5 * (M_PI / 180.0);
	constexpr unsigned int min_cluster_size = 50;
	unsigned int max_cluster_size = 250000;
	unsigned int max_large_cluster_size = 500000;
	bool ShouldDisplay = false;
	if (argc >= 4)
		tolerance = atof(argv[3]);
	if (argc >= 5)
		max_cluster_size = atoi(argv[4]);
	if (argc >= 6)
		max_large_cluster_size = atoi(argv[5]);
	if (argc >= 7)
		ShouldDisplay = atoi(argv[6]) == 1;


	// Load cloud in blob format
	pcl::PCLPointCloud2 blob;
	pcl::io::loadPLYFile(infile, blob);

	pcl::PointCloud<PointT>::Ptr cloud(new pcl::PointCloud<PointT>);
	pcl::PointCloud<PointNT>::Ptr normal(new pcl::PointCloud<PointNT>);
	std::cout << "Loading point cloud...";
	pcl::fromPCLPointCloud2(blob, *cloud);
	pcl::fromPCLPointCloud2(blob, *normal, { {16,0,4},{20,4,4},{24,8,4} });
	std::cout << "done." << std::endl;

	SegPart(cloud, normal);
	DisplayCloud(ShouldDisplay, cloud);
	SegPart(cloud, normal);
	DisplayCloud(ShouldDisplay, cloud);

	sor_filter<PointT>(cloud, 15);



#if USE_GPU
	Euclidean(cloud, tolerance, min_cluster_size, max_cluster_size, cluster_indices);
#else

	pcl::KdTree<PointT>::Ptr tree_ec(new pcl::KdTreeFLANN<PointT>());
	tree_ec->setInputCloud(cloud);
	pcl::extractEuclideanClusters(*cloud, *cloud, tolerance, tree_ec, cluster_indices, eps_angle, min_cluster_size);
#endif
	SaveCluster(cluster_indices, cloud, outfile);

	PointIndices::Ptr all_ind_ptr(new PointIndices());
	for (auto& indices : cluster_indices)
	{
		all_ind_ptr->indices.insert(all_ind_ptr->indices.end(), indices.indices.begin(), indices.indices.end());
	}

	cloud = Extract<PointT>(cloud, all_ind_ptr);
	DisplayCloud(ShouldDisplay, cloud);
	cluster_indices.clear();

	sor_filter<PointT>(cloud, 15);
	Euclidean(cloud, tolerance, min_cluster_size, max_large_cluster_size, cluster_indices);
	SaveCluster(cluster_indices, cloud, outfile, "large_cluster");
}


/*
int main(int argc, char** argv)
{

	// 读取OBJ文件
	pcl::PolygonMesh mesh;
	pcl::PCLPointCloud2 PCLcloud;
	if (pcl::io::loadPLYFile(argv[1], PCLcloud) == -1)  // 文件路径
	{
		PCL_ERROR("Couldn't read OBJ file\n");
		return -1;
	}

	// 网格顶点转换为点云
	pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
	pcl::fromPCLPointCloud2(PCLcloud, *cloud);


	pcl::PointCloud<PointNT>::Ptr doncloud (new pcl::PointCloud<PointNT>);
	copyPointCloud<pcl::PointXYZRGBNormal, PointNT>(*cloud, *doncloud);

	pcl::DifferenceOfNormalsEstimation<PointT, PointNT, PointNT> don;
	don.setInputCloud (cloud);
	don.setNormalScaleLarge (normals_large_scale);
	don.setNormalScaleSmall (normals_small_scale); don.initCompute ();
	// Compute DoN
	don.computeFeature (*doncloud);



	meshPointCloudandViewer(mesh, cloud);

	return 0;
}
*/