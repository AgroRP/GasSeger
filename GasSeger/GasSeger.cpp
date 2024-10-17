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
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/filters/voxel_grid.h>

#include <pcl/gpu/segmentation/gpu_extract_clusters.h>
#include <pcl/gpu/segmentation/gpu_extract_labeled_clusters.h>
#include <pcl/gpu/octree/octree.hpp>
// 可视化mesh模型和点云

#define USE_GPU 1

using namespace pcl;
#if USE_GPU
using PointT = pcl::PointXYZ;
#else
using PointT = pcl::PointXYZRGBNormal;
#endif
using PointNT = pcl::PointNormal;
using PointOutT = pcl::PointNormal;
using SearchPtr = pcl::search::Search<PointT>::Ptr;

void meshPointCloudandViewer(pcl::PointCloud<PointT>::Ptr& cloud)
{
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("3D PointCloud Viewer"));
	int v2;
	viewer->createViewPort(0.0, 0.0, 1.0, 1.0, v2);  // 左侧窗口
	viewer->setBackgroundColor(0.0, 0.0, 0.0, v2);  // 黑色背景
	viewer->addText("Original PointCloud", 10, 10, "vp1_text", v2);  // 标题
	pcl::visualization::PointCloudColorHandlerRGBField<PointT> cloud_color_handler(cloud);  // 绿色
	viewer->addPointCloud(cloud, cloud_color_handler, "original_cloud", v2);

	// 添加坐标系
   /* viewer->addCoordinateSystem(0.1);
	viewer->initCameraParameters();*/

	// 可视化循环
	while (!viewer->wasStopped())
	{
		viewer->spinOnce(100);
	}

}

int main(int argc, char* argv[])
{
	///The smallest scale to use in the DoN filter.
	constexpr double scale1 = 0.2;

	///The largest scale to use in the DoN filter.
	constexpr double scale2 = 2.0;

	///The minimum DoN magnitude to threshold by
	constexpr double threshold = 0.25;

	///segment scene into clusters with given distance tolerance using euclidean clustering
	double segradius = 0.2;

	//voxelization factor of pointcloud to use in approximation of normals
	bool approx = false;
	constexpr double decimation = 100;

	if (argc < 2) {
		std::cerr << "Expected 2 arguments: inputfile outputfile" << std::endl;
	}

	///The file to read from.
	std::string infile = "sampled_point_cloud.ply";

	///The file to output to.
	std::string outfile = "seg_sampled_point_cloud.ply";

	// Load cloud in blob format
	pcl::PCLPointCloud2 blob;
	pcl::io::loadPLYFile(infile, blob);

	pcl::PointCloud<PointT>::Ptr cloud(new pcl::PointCloud<PointT>);
	std::cout << "Loading point cloud...";
	pcl::fromPCLPointCloud2(blob, *cloud);
	std::cout << "done." << std::endl;

	//meshPointCloudandViewer(cloud);


	// Extracting Euclidean clusters using cloud and its normals
	std::vector<pcl::PointIndices> cluster_indices;
	constexpr float tolerance = 0.5f; // 50cm tolerance in (x, y, z) coordinate system
	constexpr double eps_angle = 5 * (M_PI / 180.0); // 5degree tolerance in normals
	constexpr unsigned int min_cluster_size = 50;



#if USE_GPU
	constexpr unsigned int max_cluster_size = 5000;
	pcl::gpu::Octree::PointCloud cloud_device;
	cloud_device.upload(cloud->points);

	pcl::gpu::Octree::Ptr octree_device(new pcl::gpu::Octree);
	octree_device->setCloud(cloud_device);
	octree_device->build();

	pcl::gpu::EuclideanClusterExtraction<pcl::PointXYZ> gec;
	gec.setClusterTolerance(0.02); // 2cm
	gec.setMinClusterSize(100);
	gec.setMaxClusterSize(25000);
	gec.setSearchMethod(octree_device);
	gec.setHostCloud(cloud);
	gec.extract(cluster_indices);
#else

	pcl::KdTree<PointT>::Ptr tree_ec(new pcl::KdTreeFLANN<PointT>());
	tree_ec->setInputCloud(cloud);
	pcl::extractEuclideanClusters(*cloud, *cloud, tolerance, tree_ec, cluster_indices, eps_angle, min_cluster_size);
#endif


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
		ss << "./Clusters/cloud_cluster_" << j << ".ply";
		writer.write<PointT>(ss.str(), *cloud_cluster, false);
		++j;
	}
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