#include "GasSeger.h"

namespace fs = boost::filesystem;
using namespace pcl;


void VisualizeCurve(
	ViewPTR viewer,
	std::string cloud_id,
	PointCloud<PointCT>::Ptr cloud,
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
		os << cloud_id << "cur_line_" << r << "_" << g << "_" << b << "_" << i;
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
	ViewPTR viewer,
	std::string cloud_id,
	int id,
	ON_3dPoint& Point,
	ON_3dVector& Vector)
{
	auto VectorEnd = Point + 0.1 * Vector;
	std::ostringstream os;
	os << cloud_id << "_Tangent_" << id;
	viewer->addLine<PointCT>(
		PointCT(Point.x, Point.y, Point.z, 255, 255, 255),
		PointCT(VectorEnd.x, VectorEnd.y, VectorEnd.z, 255, 255, 255),
		1.0,
		1.0,
		0,
		os.str());
}

void GaussianSeger::ShowClouds()
{
	InitCloudColor();
	for (auto& cloud : clouds)
	{
		auto color = cloudColors[cloud.first];
		viewer->addPointCloud(cloud.second, cloud.first, ViewportID);
		viewer->setPointCloudRenderingProperties(
			visualization::PCL_VISUALIZER_COLOR,
			color[0],
			color[1],
			color[2],
			cloud.first);
	}
}


//bool bShowCylinder = false;
//bool bShowCurve = true;


//std::vector<PointCloud<PointT>::Ptr> CachedClouds;
//std::vector<PointCloud<PointCT>::Ptr> CachedColorClouds;
//ON_NurbsCurve Stem_Curve;
//Eigen::Vector3f Stem_Center;


void GaussianSeger::InitializeViewer()
{
	viewer->createViewPort(0.0, 0.0, 1.0, 1.0, ViewportID);
	viewer->setBackgroundColor(0.0, 0.0, 0.0, ViewportID);
	viewer->addText("Original PointCloud", 10, 10, "Clusters Viewer", ViewportID);

	viewer->addText("statistical info", 0, 300, 18, 1, 0, 1, "static_info");
}

void GaussianSeger::LaunchViewer()
{
	InitializeViewer();

	std::map<std::string, std::string> static_infos;
	vtkActor* LastActor = nullptr;
	std::string LastActor_name = "Null";
	double SelectionDistance = 100;


	// select cloud
	viewer->registerMouseCallback(
		[&LastActor,
		&LastActor_name,
		&SelectionDistance,
		&static_infos,
		this](
			const visualization::MouseEvent& event)
		{
			if (event.getButton() == visualization::MouseEvent::LeftButton &&
				event.getType() == visualization::MouseEvent::MouseButtonPress)
			{
				std::vector<visualization::Camera> Cams;
				viewer->getCameras(Cams);

				auto Cam = Cams[0];
				double* WS = Cam.window_size;
				double offset[2] = {
					tan(Cam.fovy * 0.5) * (event.getY() / WS[1] - 0.5) * 2.0,
					tan(Cam.fovy * 0.5) * ((event.getX() - 0.5 * WS[0]) / WS[1]) *
					2.0
				};
				double dir[3] = {
					(Cam.focal[0] - Cam.pos[0]),
					(Cam.focal[1] - Cam.pos[1]),
					(Cam.focal[2] - Cam.pos[2])
				};
				double len = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
				double forward_dir[3] = { dir[0] / len, dir[1] / len, dir[2] / len };
				double up[3] = { Cam.view[0], Cam.view[1], Cam.view[2] };
				double right[3] = {
					forward_dir[1] * up[2] - forward_dir[2] * up[1],
					forward_dir[2] * up[0] - forward_dir[0] * up[2],
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
				if (event.getKeyboardModifiers() == visualization::KeyboardEvent::Alt ||
					event.getKeyboardModifiers() == visualization::KeyboardEvent::Ctrl)
					for (auto ele : *viewer->getCloudActorMap())
					{
						if (ele.first[0] == '_' || ele.first[0] == 'g')
							continue;
						double HitPos[3];
						double t;
						double* bound = ele.second.actor->GetBounds();
						double bounds[6] = {
							bound[0], bound[1], bound[2], bound[3], bound[4], bound[5]
						};
						if (vtkBox::IntersectBox(
							ele.second.actor->GetBounds(),
							Cam.pos,
							qdir,
							HitPos,
							t))
						{
							actor_name = ele.first;
							actor = ele.second.actor;
							break;
						}
					}
				if (event.getKeyboardModifiers() == visualization::KeyboardEvent::Ctrl)
					if (LastActor)
					{
						auto last_cloud = viewer->getCloudActorMap()->find(LastActor_name);
						if (last_cloud != viewer->getCloudActorMap()->end())
						{
							auto color = cloudColors[LastActor_name];
							viewer->setPointCloudRenderingProperties(
								visualization::PCL_VISUALIZER_COLOR,
								color[0],
								color[1],
								color[2],
								LastActor_name);
							LastActor = nullptr;
							LastActor_name = "Null";
						}
						else
						{
							auto color = cloudColors[LastActor_name];
							viewer->addPointCloud(clouds[LastActor_name], LastActor_name, ViewportID);
							viewer->setPointCloudRenderingProperties(
								visualization::PCL_VISUALIZER_COLOR,
								color[0],
								color[1],
								color[2],
								LastActor_name);
						}
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
					auto info = static_infos.find(LastActor_name);
					if (info != static_infos.end())
					{
						viewer->updateText(info->second, 0, 300, "static_info");
					}
					else
					{
						viewer->updateText("no valid statistical info", 0, 300, "static_info");
					}
					std::cout << actor_name << std::endl;
					auto color = cloudColors[actor_name];

					viewer->setPointCloudRenderingProperties(
						visualization::PCL_VISUALIZER_COLOR,
						std::min(color[0] + 0.1, 1.),
						std::min(color[1] + 0.1, 1.),
						std::min(color[2] + 0.1, 1.),
						actor_name);
				}
			}
		});

	// do calculation 
	bool bCal = false;
	viewer->registerKeyboardCallback(
		[&bCal, &LastActor, this, &LastActor_name, &static_infos](
			const visualization::KeyboardEvent event)
		{
			if (event.keyDown() && event.getKeyCode() == 'n')
			{
				if (LastActor && !bCal)
				{
					std::string actor_name = LastActor_name;
					if(actor_name.starts_with("Stem"))
					{
						static_infos["Stem.pcd"] = CalculateStem(clouds["Stem.pcd"], "Stem", true);
						return;
					}
					bCal = true;
					auto center = LastActor->GetCenter();

					std::vector<visualization::Camera> Cams;
					viewer->getCameras(Cams);
					auto Cam = Cams[0];
					double up[3] = { Cam.view[0], Cam.view[1], Cam.view[2] };
					double* cam_pos = LastActor->GetCenter();
					double text_offset = 10;
					double textpos[3] = {
						up[0] * text_offset + cam_pos[0],
						up[1] * text_offset + cam_pos[1],
						up[2] * text_offset + cam_pos[2]
					};
					double* boundary = LastActor->GetBounds();
					std::ostringstream os;
					os << "_Cube_" << actor_name;

					auto color = cloudColors[LastActor_name];
					viewer->addCube(boundary[0],
						boundary[1],
						boundary[2],
						boundary[3],
						boundary[4],
						boundary[5],
						color[0],
						color[1],
						color[2],
						os.str());
					auto cloud = clouds[actor_name];
					viewer->setShapeRenderingProperties(
						visualization::PCL_VISUALIZER_REPRESENTATION,
						visualization::PCL_VISUALIZER_REPRESENTATION_WIREFRAME,
						os.str());
					viewer->removePointCloud(actor_name);
					auto info = CalculatePart(cloud, actor_name, false);

					static_infos[actor_name] = info;
					bCal = false;
				}
			}
			else if (event.keyDown() && event.getKeyCode() == 'd')
			{
				DumpNodes();
			}
		});

	// switch visual mode
	bool bWire = false;
	viewer->registerKeyboardCallback(
		[&bWire, this](
			const visualization::KeyboardEvent event)
		{
			if (event.keyDown() && event.getKeyCode() == 'a' && !event.isAltPressed())
			{
				for (auto& actor : *viewer->getShapeActorMap())
				{
					if (!bWire)
						viewer->setShapeRenderingProperties(
							visualization::PCL_VISUALIZER_REPRESENTATION,
							visualization::PCL_VISUALIZER_REPRESENTATION_WIREFRAME,
							actor.first);
					else
						viewer->setShapeRenderingProperties(
							visualization::PCL_VISUALIZER_REPRESENTATION,
							visualization::PCL_VISUALIZER_REPRESENTATION_SURFACE,
							actor.first);
					bWire = !bWire;
				}
			}
			else if (event.keyDown() &&
				(event.getKeyCode() == 'b' || event.getKeyCode() == 'm'))
			{
				for (auto& actor : *viewer->getCloudActorMap())
				{
					double psize;
					viewer->getPointCloudRenderingProperties(
						visualization::PCL_VISUALIZER_POINT_SIZE,
						psize,
						actor.first);
					viewer->setPointCloudRenderingProperties(
						visualization::PCL_VISUALIZER_POINT_SIZE,
						psize + (event.getKeyCode() == 'b' ? 1 : -1),
						actor.first);
				}
			}
			else if (event.keyDown() && event.getKeyCode() == 'a' && event.isAltPressed())
			{
				viewer->removeAllShapes();
				viewer->addText("statistical info", 0, 300, 18, 1, 0, 1, "static_info");
			}
		});
	ShowClouds();

	//static_infos["Stem.pcd"] = CalculateStem(clouds["Stem.pcd"], "Stem", true);

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
	GaussianSeger Seger = GaussianSeger();
	for (const auto& entry : fs::directory_iterator(folderPath))
	{
		if (!is_directory(entry))
		{
			// Load cloud in blob format
			auto file_name = entry.path().filename().string();
			if (file_name.back() != 'd')
				continue;
			if (file_name[0] == 'l'||file_name[0] == 'L')
				//continue;
				Seger.LeafNum += 1;

			PCLPointCloud2 blob;
			io::loadPCDFile(entry.path().string(), blob);
			PointCloud<PointT>::Ptr cloud(new PointCloud<PointT>);
			std::cout << "Loading " << entry.path().filename().string() << std::endl;
			fromPCLPointCloud2(blob, *cloud);
			std::cout << "done.\n";
			Seger.AddCloud(entry.path().filename().string(), cloud);
		}
	}
	time_t cur_time;
	time(&cur_time);
	srand(cur_time);
	Seger.LaunchViewer();
}
