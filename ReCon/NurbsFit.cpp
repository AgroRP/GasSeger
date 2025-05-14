#include <filesystem>
#include <vtkBox.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/box_clipper3D.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/surface/on_nurbs/fitting_curve_pdm.h>
#include <pcl/surface/on_nurbs/triangulation.h>
#include <pcl/visualization/pcl_visualizer.h>

namespace fs = std::filesystem;
using namespace pcl;

pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("Curve Fitting 3D"));
using PointT = pcl::PointXYZ;
using PointCT = pcl::PointXYZRGB;

void PointCloud2Vector2d(
	pcl::PointCloud<PointT>::Ptr cloud,
	pcl::on_nurbs::vector_vec3d& data)
{
	for (const auto& p : *cloud)
	{
		if (!std::isnan(p.x) && !std::isnan(p.y) && !std::isnan(p.z))
			data.emplace_back(p.x, p.y, p.z);
	}
}


void VisualizeCurve(
	pcl::PointCloud<PointCT>::Ptr cloud,
	double r,
	double g,
	double b,
	bool show_cps)
{
	for (std::size_t i = 0; i < cloud->size() - 1; i++)
	{
		PointCT& p1 = cloud->at(i);
		PointCT& p2 = cloud->at(i + 1);
		std::ostringstream os;
		os << "line_" << r << "_" << g << "_" << b << "_" << i;
		viewer->addLine<PointCT>(p1, p2, r, g, b, os.str());
	}
}

void VisualDir(
	double x,
	double y,
	double z,
	double x2,
	double y2,
	double z2,
	double r,
	double g,
	double b,
	std::string id)
{
}

void VisualTagent(
	int id,
	ON_3dPoint& Point,
	ON_3dVector& Vector)
{
	auto VectorEnd = Point + 0.1 * Vector;
	std::ostringstream os;
	os << "Tangent_" << id;
	viewer->addLine<PointCT>(PointCT(Point.x, Point.y, Point.z, 255, 255, 255),
	                         PointCT(VectorEnd.x, VectorEnd.y, VectorEnd.z, 255, 255, 255),
	                         1.0,
	                         1.0,
	                         0,
	                         os.str());
}

class CylinderClipper
{
public:
	double range;
	double radius;
	Eigen::Vector3f pos;
	Eigen::Vector3f dir;

	CylinderClipper()
	{
	};

	bool clipPoint3D(
		const PointT point)
	{
		auto p = pos - point.getVector3fMap();
		return p.norm() < radius && abs(dir.dot(p)) < range;
	}

	void Clip3DPoints(
		pcl::PointCloud<PointT>::Ptr& cloud,
		pcl::Indices& indices)
	{
		indices.clear();
		pcl::PointCloud<PointT>& cloud_in = *cloud;
		if (indices.empty())
		{
			indices.reserve(cloud->size());
			for (std::size_t pIdx = 0; pIdx < cloud->size(); ++pIdx)
				if (clipPoint3D(cloud_in[pIdx]))
					indices.push_back(pIdx);
		}
	}
};

bool bShowCylinder = false;
bool bShowCurve = false;

std::vector<pcl::PointCloud<PointT>::Ptr> CachedClouds;
std::vector<pcl::PointCloud<PointCT>::Ptr> CachedColorClouds;

void CachePtr(
	pcl::PointCloud<PointT>::Ptr cloud)
{
	CachedClouds.push_back(cloud);
}

void CachePtr(
	pcl::PointCloud<PointCT>::Ptr cloud)
{
	CachedColorClouds.push_back(cloud);
}


int CalculatePart(
	pcl::PointCloud<PointT>::Ptr cloud)
{
	pcl::VoxelGrid<PointT> voxel_grid;
	const float leaf_size = 0.01f;
	voxel_grid.setLeafSize(leaf_size, leaf_size, leaf_size);

	voxel_grid.setInputCloud(cloud);
	pcl::PointCloud<PointT>::Ptr downsampled_cloud(new pcl::PointCloud<PointT>);
	CachePtr(downsampled_cloud);
	voxel_grid.filter(*downsampled_cloud);

	viewer->addPointCloud<PointT>(downsampled_cloud, "down cloud");


	// convert to NURBS data structure
	pcl::on_nurbs::NurbsDataCurve data;
	PointCloud2Vector2d(downsampled_cloud, data.interior);

	// #################### CURVE PARAMETERS #########################
	unsigned order(3);
	unsigned n_control_points(20);

	pcl::on_nurbs::FittingCurve::Parameter curve_params;
	curve_params.smoothness = 0.000001;

	// #################### CURVE FITTING #########################
	ON_NurbsCurve curve = pcl::on_nurbs::FittingCurve::initNurbsCurvePCA(
		order,
		data.interior,
		n_control_points);

	pcl::on_nurbs::FittingCurve fit(&data, curve);
	fit.assemble(curve_params);
	fit.solve();
	if (bShowCurve)
	{
		pcl::PointCloud<PointCT>::Ptr CurveCloud(new pcl::PointCloud<PointCT>);
		pcl::on_nurbs::Triangulation::convertCurve2PointCloud(fit.m_nurbs, CurveCloud, 8);
		VisualizeCurve(CurveCloud, 1.0, 0.0, 0.0, false);
		CachePtr(CurveCloud);
	}


	int id = 0;
	auto UP = Eigen::Vector3f(0, 0, 1);
	double r, g, b;
	r = 0.7;
	g = 0.2;
	b = 0.0;
	std::vector<pcl::PointCloud<PointT>::Ptr> SelectedClouds;
	std::vector<Eigen::Vector3f> SelectedTangent;
	std::vector<Eigen::Vector3f> SelectedPos;
	std::vector<Eigen::Vector3f> SelectedCurvature;
	int cp_red = fit.m_nurbs.Order() - 2;
	int resolution = 8;
	bool bLastValid;
	int StartIndex;
	float curvatues = 0;
	for (int i = 1; i < fit.m_nurbs.KnotCount() - 1 - cp_red; i++)
	{
		double dr = 1.0 / (resolution - 1);
		double xi0 = fit.m_nurbs.m_knot[i];
		double xid = (fit.m_nurbs.m_knot[i + 1] - xi0);

		for (unsigned j = 0; j < resolution; j++, id++)
		{
			double t = (xi0 + dr * xid * j);
			ON_3dPoint Point;
			ON_3dVector Tangent;
			fit.m_nurbs.EvTangent(t, Point, Tangent);


			pcl::ExtractIndices<PointT> extract_indices;
			pcl::Indices indices;
			Eigen::Translation3f Translation(Point.x, Point.y, Point.z);
			Eigen::Vector3f Scale = Eigen::Vector3f({1, 1, 0.01});
			auto vec_T = Eigen::Vector3f(Tangent.x, Tangent.y, Tangent.z);
			vec_T.normalize();
			auto axis = UP.cross(vec_T);
			auto right = axis.cross(UP);
			double inter = UP.dot(vec_T);
			double dir_dot = right.dot(vec_T);
			auto rotation = Eigen::AngleAxisf(dir_dot > 0 ? acos(inter) : -acos(inter), axis);
			auto Quaternion = Eigen::Quaternionf(rotation);


			auto partclipper = CylinderClipper();
			partclipper.radius = 1;
			partclipper.range = 0.005;
			partclipper.dir = vec_T;
			partclipper.pos = Translation.vector();
			partclipper.Clip3DPoints(downsampled_cloud, indices);
			if (indices.size() > 0)
			{
				VisualTagent(id, Point, Tangent);

				r = (rand() % 65536) / 140000.f + 0.2f;
				g = (rand() % 65536) / 140000.f + 0.2f;
				b = (rand() % 65536) / 140000.f + 0.2f;

				pcl::PointCloud<PointT>::Ptr cloud_out(new pcl::PointCloud<PointT>);
				extract_indices.setInputCloud(downsampled_cloud);
				extract_indices.setIndices(pcl::make_shared<pcl::Indices>(indices));
				extract_indices.filter(*cloud_out);
				std::ostringstream os;
				if (bShowCylinder)
				{
					os << "Cube_" << id;
					auto cube_name = os.str();
					viewer->addCube(Translation.vector(), Quaternion, 1, 1, 0.05, cube_name);
					viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, r, g, b, cube_name);
					viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, 0.5, cube_name);
					os.clear();
				}
				os << "Cloud_" << id;
				auto cloud_name = os.str();
				std::cout << cloud_name << "size ::" << cloud_out->size() << std::endl;
				viewer->addPointCloud(cloud_out, cloud_name);
				viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, r, g, b, cloud_name);
				SelectedClouds.push_back(cloud_out);
				SelectedTangent.push_back(vec_T);
				SelectedPos.push_back(Translation.vector());
				auto curvature = fit.m_nurbs.CurvatureAt(t);
				curvatues += curvature.Length();
				SelectedCurvature.push_back(Eigen::Vector3f(curvature.x, curvature.y, curvature.z));
				if (!bLastValid)
				{
					bLastValid = true;
					StartIndex = SelectedClouds.size() - 1;
				}
			}
			else
			{
				bLastValid = false;
			}
			viewer->spinOnce(1);
		}
	}

	int AllPartNum = SelectedClouds.size();
	double length = 0;
	Eigen::Vector3f LastPos;
	Eigen::Vector3f LastTangent;
	Eigen::Vector3f LastStart;
	Eigen::Vector3f LastEnd;
	Eigen::Vector3d LastNormal;
	double total_torsion = 0.;
	double total_area = 0.;
	double total_width = 0.;
	for (int i = 0; i < AllPartNum; i++)
	{
		int ind = (StartIndex + i) % AllPartNum;
		pcl::on_nurbs::NurbsDataCurve CurData;
		auto CurCloud = SelectedClouds[ind];
		PointCloud2Vector2d(CurCloud, CurData.interior);
		Eigen::Vector3d mean;
		Eigen::Matrix3d eigenvectors;
		Eigen::Vector3d eigenvalues;
		pcl::on_nurbs::NurbsTools::pca(CurData.interior, mean, eigenvectors, eigenvalues);
		Eigen::Vector3d CurNormal;
		double min_v = 10000;
		double max_v = -10000;
		if (eigenvalues[0] > eigenvalues[1] && eigenvalues[0] > eigenvalues[2])
		{
			CurNormal = eigenvectors.col(0);
		}
		else if (eigenvalues[1] > eigenvalues[0] && eigenvalues[1] > eigenvalues[2])
		{
			CurNormal = eigenvectors.col(1);
		}
		else
		{
			CurNormal = eigenvectors.col(2);
		}
		CurNormal.normalize();
		Eigen::Vector3f CurNormalf(CurNormal.x(), CurNormal.y(), CurNormal.z());
		auto CurPos = SelectedPos[ind];

		for (std::size_t pIdx = 0; pIdx < CurCloud->size(); ++pIdx)
		{
			double v = ((*CurCloud)[pIdx].getVector3fMap() - CurPos).dot(CurNormalf);
			if (v > max_v)
				max_v = v;
			if (v < min_v)
				min_v = v;
		}

		std::ostringstream os;
		os << "Normal_" << i;
		total_width += max_v - min_v;
		auto cube_name = os.str();
		Eigen::Vector3f VectorEnd = CurPos + max_v * CurNormalf;
		Eigen::Vector3f VectorStart = CurPos + min_v * CurNormalf;
		viewer->addLine<PointCT>(PointCT(VectorStart.x(), VectorStart.y(), VectorStart.z(), 255, 255, 255),
		                         PointCT(VectorEnd.x(), VectorEnd.y(), VectorEnd.z(), 255, 255, 255),
		                         1.0,
		                         0.0,
		                         1.0,
		                         os.str());


		auto CurTangent = SelectedTangent[ind];
		if (i > 0)
		{
			auto CurLength = (LastPos - CurPos).norm();
			length += CurLength;
			total_torsion += 1 - abs(LastNormal.dot(CurNormal));
			total_area += ((VectorEnd - VectorStart).cross(LastStart - VectorStart)).norm();
			total_area += ((LastEnd - LastStart).cross(VectorEnd - LastStart)).norm();
		}
		LastPos = CurPos;
		LastNormal = CurNormal;
		LastTangent = CurTangent;
		LastStart = VectorStart;
		LastEnd = VectorEnd;
	}


	std::ostringstream static_info;
	static_info << "average curvature	::  " << curvatues / AllPartNum << "\n";
	static_info << "total length		::  " << length << "\n";
	static_info << "total torsion		::  " << total_torsion << "\n";
	static_info << "average torsion		::  " << total_torsion / AllPartNum << "\n";
	static_info << "total area			::  " << total_area / 2.0 << "\n";
	static_info << "average width		::  " << total_width / AllPartNum << "\n";
	viewer->removePointCloud("down cloud");
	viewer->updateText(static_info.str(), 0, 300, "static_info");

	for (auto& ptr : SelectedClouds)
	{
		CachePtr(ptr);
	}
	return 0;
}


void LaunchViewer(
	std::map<std::string, pcl::PointCloud<PointT>::Ptr> clouds)
{
	int v2;
	viewer->createViewPort(0.0, 0.0, 1.0, 1.0, v2);
	viewer->setBackgroundColor(0.0, 0.0, 0.0, v2);
	viewer->addText("Original PointCloud", 10, 10, "Clusters Viewer", v2);
	int ID = 0;
	double r, g, b;
	vtkActor* LastActor = nullptr;
	std::string LastActor_name = "Null";
	double SelectionDistance = 100;

	viewer->registerMouseCallback([&v2, &LastActor, &r, &b, &g, &LastActor_name, &SelectionDistance](
		const pcl::visualization::MouseEvent& event)
		{
			if (event.getButton() == visualization::MouseEvent::LeftButton && event.getType() ==
				visualization::MouseEvent::MouseButtonPress)
			{
				std::vector<visualization::Camera> Cams;
				viewer->getCameras(Cams);

				auto Cam = Cams[0];
				double* WS = Cam.window_size;
				double offset[2] = {
					tan(Cam.fovy * 0.5) * (event.getY() / WS[1] - 0.5) * 2.0,
					tan(Cam.fovy * 0.5) * ((event.getX() - 0.5 * WS[0]) / WS[1]) * 2.0
				};
				double dir[3] = {(Cam.focal[0] - Cam.pos[0]), (Cam.focal[1] - Cam.pos[1]), (Cam.focal[2] - Cam.pos[2])};
				double len = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
				double forward_dir[3] = {dir[0] / len, dir[1] / len, dir[2] / len};
				double up[3] = {Cam.view[0], Cam.view[1], Cam.view[2]};
				double right[3] = {
					forward_dir[1] * up[2] - forward_dir[2] * up[1], forward_dir[2] * up[0] - forward_dir[0] * up[2],
					forward_dir[0] * up[1] - forward_dir[1] * up[0]
				};
				double qdir[3] = {
					forward_dir[0] + offset[0] * up[0] + offset[1] * right[0],
					forward_dir[1] + offset[0] * up[1] + offset[1] * right[1],
					forward_dir[2] + offset[0] * up[2] + offset[1] * right[2]
				};
				qdir[0] *= SelectionDistance;
				qdir[1] *= SelectionDistance;
				qdir[2] *= SelectionDistance;
				vtkActor* actor = nullptr;
				std::string actor_name = "Null";
				if (event.getKeyboardModifiers() == visualization::KeyboardEvent::Alt || event.getKeyboardModifiers() ==
					visualization::KeyboardEvent::Ctrl)
					for (auto ele : *viewer->getCloudActorMap())
					{
						double HitPos[3];
						double t;
						double* bound = ele.second.actor->GetBounds();
						double bounds[6] = {bound[0], bound[1], bound[2], bound[3], bound[4], bound[5]};;
						if (vtkBox::IntersectBox(ele.second.actor->GetBounds(), Cam.pos, qdir, HitPos, t))
						{
							actor_name = ele.first;
							actor = ele.second.actor;
							break;
						}
					}
				if (event.getKeyboardModifiers() == visualization::KeyboardEvent::Ctrl)
					if (LastActor)
					{
						viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR,
						                                         r,
						                                         g,
						                                         b,
						                                         LastActor_name);
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
					viewer->getPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR,
					                                         r,
					                                         g,
					                                         b,
					                                         actor_name);
					viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR,
					                                         std::min(r + 0.1, 1.),
					                                         std::min(g + 0.1, 1.),
					                                         std::min(b + 0.1, 1.),
					                                         actor_name);
				}
			}
		});
	bool bCal = false;
	viewer->registerKeyboardCallback([&bCal, &LastActor, &clouds, &LastActor_name](
		const pcl::visualization::KeyboardEvent event)
		{
			if (event.keyDown() && event.getKeyCode() == 'n')
			{
				if (LastActor && !bCal)
				{
					bCal = true;
					auto center = LastActor->GetCenter();
					viewer->setCameraPosition(center[0],
					                          center[1],
					                          center[2] + 5,
					                          center[0],
					                          center[1],
					                          center[2],
					                          0,
					                          1,
					                          0);
					std::vector<pcl::visualization::Camera> Cams;
					viewer->getCameras(Cams);
					auto Cam = Cams[0];
					double up[3] = {Cam.view[0], Cam.view[1], Cam.view[2]};
					double* cam_pos = LastActor->GetCenter();
					double text_offset = 10;
					double textpos[3] = {
						up[0] * text_offset + cam_pos[0], up[1] * text_offset + cam_pos[1],
						up[2] * text_offset + cam_pos[2]
					};
					viewer->addText("static_info", 0, 300, 12, 1, 0, 1, "static_info");
					double* boundary = LastActor->GetBounds();
					viewer->addCube(boundary[0], boundary[1], boundary[2], boundary[3], boundary[4], boundary[5]);
					auto cloud = clouds[LastActor_name];
					CalculatePart(cloud);
					bCal = false;
				}
			}
		});
	bool bWire = false;
	viewer->registerKeyboardCallback([&bWire](
		const pcl::visualization::KeyboardEvent event)
		{
			if (event.keyDown() && event.getKeyCode() == 'a')
			{
				for (auto& actor : *viewer->getShapeActorMap())
				{
					if (!bWire)
						viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_REPRESENTATION,
						                                    pcl::visualization::PCL_VISUALIZER_REPRESENTATION_WIREFRAME,
						                                    actor.first);
					else
						viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_REPRESENTATION,
						                                    pcl::visualization::PCL_VISUALIZER_REPRESENTATION_SURFACE,
						                                    actor.first);
					bWire = !bWire;
				}
			}
			if (event.keyDown() && (event.getKeyCode() == 'b' || event.getKeyCode() == 'm'))
			{
				for (auto& actor : *viewer->getCloudActorMap())
				{
					double psize;
					viewer->getPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE,
					                                         psize,
					                                         actor.first);
					viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE,
					                                         psize + (event.getKeyCode() == 'b' ? 1 : -1),
					                                         actor.first);
				}
			}
		});
	for (auto& cloud : clouds)
	{
		float r = (rand() % 65536) / 140000.f + 0.2f;
		float g = (rand() % 65536) / 140000.f + 0.2f;
		float b = (rand() % 65536) / 140000.f + 0.2f;
		viewer->addPointCloud(cloud.second, cloud.first, v2);
		viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, r, g, b, cloud.first);
	}

	while (!viewer->wasStopped())
	{
		viewer->spinOnce(100);
	}
}

int main(
	int argc,
	char* argv[])
{
	fs::path folderPath = argv[1];
	std::map<std::string, pcl::PointCloud<PointT>::Ptr> PCs;
	for (const auto& entry : fs::directory_iterator(folderPath))
	{
		if (!fs::is_directory(entry))
		{
			// Load cloud in blob format
			pcl::PCLPointCloud2 blob;
			pcl::io::loadPCDFile(entry.path().string(), blob);
			pcl::PointCloud<PointT>::Ptr cloud(new pcl::PointCloud<PointT>);
			std::cout << "Loading " << entry.path().filename().string() << std::endl;
			pcl::fromPCLPointCloud2(blob, *cloud);
			std::cout << "done.\n";
			PCs.emplace(entry.path().filename().string(), cloud);
			break;
		}
	}

	LaunchViewer(PCs);
}
