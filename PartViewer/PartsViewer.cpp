#include "PartsViewer.h"

#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/io/obj_io.h>  
#include <pcl/io/ply_io.h>  
#include <pcl/PolygonMesh.h>
#include <pcl/conversions.h>  
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/features/don.h>
#include <string>
#include <vtkPicker.h>
#include <vtkPropPicker.h>
#include <vtkBox.h>
#include <vtkVector.h>
#include <iostream>
#include <filesystem>
#include <vtkRenderer.h>

namespace fs = std::filesystem;
char IntersectBox(const double bounds[6], const double origin[3], const double dir[3],
	double coord[3], double& t, double tolerance = 0.0)
{
	bool inside = true;
	char quadrant[3];
	int i, whichPlane = 0;
	double maxT[3], candidatePlane[3];

	// Make sure tolerance is non-zero
	double tol = (tolerance <= 0 ? FLT_EPSILON : tolerance);

	// Make sure bounds are not degenerate (i.e., have positive volume); pad by
	// tolerance if it is.
	double bds[6];
	for (i = 0; i < 3; ++i)
	{
		if ((bounds[2 * i + 1] - bounds[2 * i]) > 0)
		{
			bds[2 * i] = bounds[2 * i];
			bds[2 * i + 1] = bounds[2 * i + 1];
		}
		else
		{
			bds[2 * i] = bounds[2 * i] - tol;
			bds[2 * i + 1] = bounds[2 * i + 1] + tol;
		}
	}

	//  First find closest planes
	//
	for (i = 0; i < 3; i++)
	{
		if (origin[i] < bds[2 * i])
		{
			quadrant[i] = 0;
			candidatePlane[i] = bds[2 * i];
			inside = false;
		}
		else if (origin[i] > bds[2 * i + 1])
		{
			quadrant[i] = 1;
			candidatePlane[i] = bds[2 * i + 1];
			inside = false;
		}
		else
		{
			quadrant[i] = 2;
		}
	}

	//  Check whether origin of ray is inside bbox
	//
	if (inside)
	{
		coord[0] = origin[0];
		coord[1] = origin[1];
		coord[2] = origin[2];
		t = 0;
		return 1;
	}

	//  Calculate parametric distances to plane
	//
	for (i = 0; i < 3; i++)
	{
		if (quadrant[i] != 2 && dir[i] != 0.0)
		{
			maxT[i] = (candidatePlane[i] - origin[i]) / dir[i];
		}
		else
		{
			maxT[i] = -1.0;
		}
	}

	//  Find the largest parametric value of intersection
	//
	for (i = 0; i < 3; i++)
	{
		if (maxT[whichPlane] < maxT[i])
		{
			whichPlane = i;
		}
	}

	//  Check for valid intersection along line
	//
	if (maxT[whichPlane] > 1.0 || maxT[whichPlane] < 0.0)
	{
		return 0;
	}
	else
	{
		t = maxT[whichPlane];
	}

	//  Intersection point along line is okay.  Check bbox.
	//
	for (i = 0; i < 3; i++)
	{
		if (whichPlane != i)
		{
			coord[i] = origin[i] + maxT[whichPlane] * dir[i];
			if ((coord[i] < bds[2 * i] - tol) || (coord[i] > bds[2 * i + 1] + tol))
			{
				return 0;
			}
		}
		else
		{
			coord[i] = candidatePlane[i];
		}
	}

	return 1;
}

using namespace pcl;
using PointT = pcl::PointXYZ;

void ClusterViewer(std::vector<pcl::PointCloud<PointT>::Ptr> clouds)
{
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("3D PointCloud Viewer"));
	int v2;
	viewer->createViewPort(0.0, 0.0, 1.0, 1.0, v2);
	viewer->setBackgroundColor(0.0, 0.0, 0.0, v2);
	viewer->addText("Original PointCloud", 10, 10, "Clusters Viewer", v2);
	int ID = 0;
	double r, g, b;
	vtkActor* LastActor = nullptr;
	std::string LastActor_name = "Null";
	double SelectionDistance = 100;

	viewer->registerMouseCallback([&viewer, &v2, &LastActor, &r, &b, &g, &LastActor_name, &SelectionDistance](const pcl::visualization::MouseEvent& event)
		{
			if (event.getButton() == visualization::MouseEvent::LeftButton && event.getType() == visualization::MouseEvent::MouseButtonPress)
			{
				std::vector<visualization::Camera> Cams;
				viewer->getCameras(Cams);

				auto Cam = Cams[0];
				double* WS = Cam.window_size;
				double offset[2] = { tan(Cam.fovy * 0.5) * (event.getY() / WS[1] - 0.5) * 2.0,tan(Cam.fovy * 0.5) * ((event.getX() - 0.5 * WS[0]) / WS[1]) * 2.0 };
				double dir[3] = { (Cam.focal[0] - Cam.pos[0]),(Cam.focal[1] - Cam.pos[1]),(Cam.focal[2] - Cam.pos[2]) };
				double len = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
				double forward_dir[3] = { dir[0] / len,dir[1] / len,dir[2] / len };
				double up[3] = { Cam.view[0],Cam.view[1],Cam.view[2] };
				double right[3] = { forward_dir[1] * up[2] - forward_dir[2] * up[1],forward_dir[2] * up[0] - forward_dir[0] * up[2],forward_dir[0] * up[1] - forward_dir[1] * up[0] };
				double qdir[3] = { forward_dir[0] + offset[0] * up[0] + offset[1] * right[0],forward_dir[1] + offset[0] * up[1] + offset[1] * right[1],forward_dir[2] + offset[0] * up[2] + offset[1] * right[2] };
				qdir[0] *= SelectionDistance;
				qdir[1] *= SelectionDistance;
				qdir[2] *= SelectionDistance;
				vtkActor* actor = nullptr;
				std::string actor_name = "Null";
				if (event.getKeyboardModifiers() == visualization::KeyboardEvent::Alt || event.getKeyboardModifiers() == visualization::KeyboardEvent::Ctrl)
					for (auto ele : *viewer->getCloudActorMap())
					{
						double HitPos[3];
						double t;
						double* bound = ele.second.actor->GetBounds();
						double bounds[6] = { bound[0],bound[1],bound[2],bound[3],bound[4],bound[5] };
						;
						if (IntersectBox(ele.second.actor->GetBounds(), Cam.pos, qdir, HitPos, t))
						{
							actor_name = ele.first;
							actor = ele.second.actor;
							break;
						}
					}
				if (event.getKeyboardModifiers() == visualization::KeyboardEvent::Ctrl)
					if (LastActor)
					{
						viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, r, g, b, LastActor_name);
						LastActor = nullptr;
						LastActor_name = "Null";
					}
				if (actor)
				{
					if (event.getKeyboardModifiers() == visualization::KeyboardEvent::Alt)
					{
						viewer->removePointCloud(actor_name);
						return;
					}
					LastActor = actor;
					LastActor_name = actor_name;
					std::cout << actor_name << std::endl;
					viewer->getPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, r, g, b, actor_name);
					viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, std::min(r + 0.1, 1.), std::min(g + 0.1, 1.), std::min(b + 0.1, 1.), actor_name);
				}
			}
		});
	for (auto& cloud : clouds)
	{
		float r = (rand() % 65536) / 140000.f + 0.2f;
		float g = (rand() % 65536) / 140000.f + 0.2f;
		float b = (rand() % 65536) / 140000.f + 0.2f;
		char sID[10];
		sprintf(sID, "%d", ID++);
		std::string actor_name = std::string("cluster_").append(sID);
		viewer->addPointCloud(cloud, actor_name, v2);
		viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, r, g, b, actor_name);
	}

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


	fs::path folderPath = argv[1];
	std::vector<pcl::PointCloud<PointT>::Ptr> PCs;
	for (const auto& entry : fs::directory_iterator(folderPath)) {
		if (!fs::is_directory(entry)) {
			// Load cloud in blob format
			pcl::PCLPointCloud2 blob;
			pcl::io::loadPLYFile(entry.path().string(), blob);
			pcl::PointCloud<PointT>::Ptr cloud(new pcl::PointCloud<PointT>);
			std::cout << "Loading " << entry.path().filename().string() << std::endl;
			pcl::fromPCLPointCloud2(blob, *cloud);
			std::cout << "done.\n";
			PCs.push_back(cloud);
		}
	}

	ClusterViewer(PCs);
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