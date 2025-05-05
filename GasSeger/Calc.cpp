#include <valarray>

#include "GasSeger.h"
#include <pcl/search/kdtree.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/extract_clusters.h>

class FixedCurve : public pcl::on_nurbs::FittingCurve
{
	void updateCurve(
		double damp) //override
	{
		int cp_red = m_nurbs.m_order - 2;
		int ncp = m_nurbs.m_cv_count - 2 * cp_red;

		for (int j = 2; j < ncp - 1; j++)
		{
			ON_3dPoint cp;
			cp.x = Solver().x(j, 0);
			cp.y = Solver().x(j, 1);
			cp.z = Solver().x(j, 2);

			m_nurbs.SetCV(j + cp_red, cp);
		}

		//for (int j = 0; j < cp_red; j++)
		//{
		//	ON_3dPoint cp;
		//	m_nurbs.GetCV(m_nurbs.m_cv_count - 1 - cp_red + j, cp);
		//	m_nurbs.SetCV(j, cp);

		//	m_nurbs.GetCV(cp_red - j, cp);
		//	m_nurbs.SetCV(m_nurbs.m_cv_count - 1 - j, cp);
		//}
	}
};

ON_NurbsCurve initNurbsCurvePCAOffset(
	int order,
	const on_nurbs::vector_vec3d& data,
	int ncps,
	double rf)
{
	if (data.empty())
		printf("[FittingCurve::initNurbsCurvePCA] Warning, no boundary parameters "
			"available\n");

	Eigen::Vector3d mean;
	Eigen::Matrix3d eigenvectors;
	Eigen::Vector3d eigenvalues;

	unsigned s = unsigned(data.size());

	on_nurbs::NurbsTools::pca(data, mean, eigenvectors, eigenvalues);

	eigenvalues /=
		s; // seems that the eigenvalues are dependent on the number of points (???)

	double r = rf * eigenvalues(0);
	double aux_r1 = rf * eigenvalues(1);
	double aux_r2 = rf * eigenvalues(2);
	bool bShouldOffset = false;
	if (aux_r1 < 0.1 * r && aux_r2 < 0.1 * r)
	{
		bShouldOffset = true;
	}

	if (ncps < 2 * order)
		ncps = 2 * order;

	ON_NurbsCurve nurbs = ON_NurbsCurve(3, false, order, ncps);
	nurbs.MakePeriodicUniformKnotVector(1.0 / (ncps - order + 1));

	double dcv = (2.0 * M_PI) / (ncps - order + 1);
	Eigen::Vector3d cv, cv_t;
	if (bShouldOffset)
	{
		dcv = (M_PI / 3);
		for (int j = 0; j < 3; j++)
		{
			cv(0) = 4 * r * sin(dcv * j - M_PI / 6);
			cv(1) = 4 * r * (std::cos(dcv * j - M_PI / 6) + sqrt(3.0) / 2.0);
			cv(2) = 0.0;
			cv_t = eigenvectors * cv + mean;
			nurbs.SetCV(j, ON_3dPoint(cv_t(0), cv_t(1), cv_t(2)));
		}
		dcv = (M_PI / 3) / (ncps - 6);
		for (int j = 0; j < ncps - 6; j++)
		{
			cv(0) = 4 * r * sin(dcv * j + 5 * M_PI / 6);
			cv(1) = 4 * r * (std::cos(dcv * j + 5 * M_PI / 6) + sqrt(3.0) / 2.0);
			cv(2) = 0.0;
			cv_t = eigenvectors * cv + mean;
			nurbs.SetCV(j + 3, ON_3dPoint(cv_t(0), cv_t(1), cv_t(2)));
		}
		dcv = (M_PI / 3);
		for (int j = 0; j < 3; j++)
		{
			cv(0) = 4 * r * sin(dcv * j + 7 * M_PI / 6);
			cv(1) = 4 * r * (std::cos(dcv * j + 7 * M_PI / 6) + sqrt(3.0) / 2.0);
			cv(2) = 0.0;
			cv_t = eigenvectors * cv + mean;
			nurbs.SetCV(ncps - order - 2 + j, ON_3dPoint(cv_t(0), cv_t(1), cv_t(2)));
		}
	}
	else
	{
		for (int j = 0; j < ncps; j++)
		{
			cv(0) = r * sin(dcv * j);
			cv(1) = r * (std::cos(dcv * j));
			cv(2) = 0.0;
			cv_t = eigenvectors * cv + mean;
			nurbs.SetCV(j, ON_3dPoint(cv_t(0), cv_t(1), cv_t(2)));
		}
	}

	// for (int j = 0; j < ncps; j++) {
	//   cv(0) = r * sin(dcv * j);
	//   cv(1) = r * (std::cos(dcv * j) + (bShouldOffset ? 2.0 : 0.0));
	//   cv(2) = 0.0;
	//   cv_t = eigenvectors * cv + mean;
	//   nurbs.SetCV(j, ON_3dPoint(cv_t(0), cv_t(1), cv_t(2)));
	// }

	return nurbs;
}


ON_NurbsCurve initNurbsCurveRef(
	int order,
	const on_nurbs::vector_vec3d& data,
	int ncps,
	on_nurbs::vector_vec3d ref)
{
	if (data.empty())
		printf("[FittingCurve::initNurbsCurvePCA] Warning, no boundary parameters "
			"available\n");

	Eigen::Vector3d mean;
	Eigen::Matrix3d eigenvectors;
	Eigen::Vector3d eigenvalues;

	unsigned s = unsigned(ref.size());

	on_nurbs::NurbsTools::pca(ref, mean, eigenvectors, eigenvalues);

	eigenvalues /=
		s; // seems that the eigenvalues are dependent on the number of points (???)

	auto offset = eigenvalues[0] * (eigenvectors.col(1).normalized());

	ncps = ref.size() + 4;
	ON_NurbsCurve nurbs = ON_NurbsCurve(3, false, order, ncps);
	nurbs.MakePeriodicUniformKnotVector(1.0 / (ncps - order + 1));

	for (int i = 0; i < ref.size(); ++i)
	{
		auto& ref_point = ref[i];
		nurbs.SetCV(i + 2, ON_3dPoint(ref_point(0), ref_point(1), ref_point(2)));
	}
	nurbs.SetCV(ncps - 2,
		ON_3dPoint(ref[ncps - 5](0) * 2.f - ref[ncps - 6](0),
			ref[ncps - 5](1) * 2.f - ref[ncps - 6](1),
			ref[ncps - 5](2) * 2.f - ref[ncps - 6](2)));
	nurbs.SetCV(ncps - 1,
		ON_3dPoint(ref[ncps - 5](0) * 5.f - ref[ncps - 6](0)*4.f,
			ref[ncps - 5](1) * 5.f - ref[ncps - 6](1)*4.f,
			ref[ncps - 5](2) * 5.f - ref[ncps - 6](2)*4.f));
	//nurbs.SetCV(ncps - 1,
		//ON_3dPoint(offset(0) + ref[ncps - 5](0), offset(1) + ref[ncps - 5](1), offset(2) + ref[ncps - 5](2)));
	nurbs.SetCV(1, ON_3dPoint(ref[0](0) * 2.f - ref[1](0), ref[0](1) * 2.f - ref[1](1), ref[0](2) * 2.f - ref[1](2)));
	nurbs.SetCV(0, ON_3dPoint(ref[0](0) * 5.f - 4.f*ref[1](0), ref[0](1) * 5.f - 4.f*ref[1](1), ref[0](2) * 5.f - 4.f*ref[1](2)));
	//nurbs.SetCV(0, ON_3dPoint(offset(0) + ref[0](0), offset(1) + ref[0](1), offset(2) + ref[0](2)));

	// for (int j = 0; j < ncps; j++) {
	//   cv(0) = r * sin(dcv * j);
	//   cv(1) = r * (std::cos(dcv * j) + (bShouldOffset ? 2.0 : 0.0));
	//   cv(2) = 0.0;
	//   cv_t = eigenvectors * cv + mean;
	//   nurbs.SetCV(j, ON_3dPoint(cv_t(0), cv_t(1), cv_t(2)));
	// }

	return nurbs;
}

void PointCloud2Vector2d(
	PointCloud<PointT>::Ptr cloud,
	on_nurbs::vector_vec3d& data)
{
	for (const auto& p : *cloud)
	{
		if (!std::isnan(p.x) && !std::isnan(p.y) && !std::isnan(p.z))
			data.emplace_back(p.x, p.y, p.z);
	}
}


GaussianSeger::GaussianSeger()
	: viewer(new visualization::PCLVisualizer("Curve Fitting 3D")),
	bShowCurve(true),
	bShowCylinder(true),
	LeafNum(0)
{
	Nodes.clear();
}

void GaussianSeger::DumpNodes()
{
	Nodes.clear();
	for (auto& node : AllNodes["Stem.pcd"])
	{
		Nodes.push_back(node);
	}

	for(auto &PartNodes:AllNodes)
	{
		if (PartNodes.first == "Stem.pcd")
		{
			continue;
		}
		int i = 0;
		for (auto& node : PartNodes.second)
		{
			if (i!=0)
			{
				node.parent_id = Nodes.size() - 1;
			}
			i++;
			node.id = Nodes.size();
			node.total_generation = PartNodes.second.size();
			Nodes.push_back(node);
		}
	}

	std::cout << "{\n\"pos\":[" << std::endl;
	bool first = true;
	for (auto& node : Nodes)
	{
		if (first)
		{
			first = false;
		}
		else
		{
			std::cout << "," << std::endl;
		}
		std::cout << "{\"x\":" << node.pos.x() << ",\"y\":" << node.pos.y() << ",\"z\":" << node.pos.z()<<"}"<< std::endl;
	}
	std::cout << "],\"dir\":[" << std::endl;
	first = true;
	for (auto& node : Nodes)
	{
		if (first)
		{
			first = false;
		}
		else
		{
			std::cout << "," << std::endl;
		}
		std::cout << "{\"x\":" << node.tangent.x() << ",\"y\":" << node.tangent.y() << ",\"z\":" << node.tangent.z() << "}" << std::endl;
	}
	std::cout << "],\"nor\":[" << std::endl;
	first = true;
	for (auto& node : Nodes)
	{
		if (first)
		{
			first = false;
		}
		else
		{
			std::cout << "," << std::endl;
		}
		std::cout << "{\"x\":" << node.normal.x() << ",\"y\":" << node.normal.y() << ",\"z\":" << node.normal.z() << "}" << std::endl;
	}
	std::cout << "],\"pro\":[" << std::endl;
	first = true;
	for (auto& node : Nodes)
	{
		if (first)
		{
			first = false;
		}
		else
		{
			std::cout << "," << std::endl;
		}
		std::cout << "{\"x\":" << node.length << ",\"y\":" << node.width << ",\"z\":" << node.type + 2 << "}" << std::endl;
	}
	std::cout << "],\"gen\":[" << std::endl;
	first = true;
	for (auto& node : Nodes)
	{
		if (first)
		{
			first = false;
		}
		else
		{
			std::cout << "," << std::endl;
		}
		std::cout << "{\"x\":" << node.total_generation << ",\"y\":" << node.generation << ",\"z\":" << node.type + 2 << "}" << std::endl;
	}
	std::cout << "],\"ids\":[" << std::endl;
	first = true;
	for (auto& node : Nodes)
	{
		if (first)
		{
			first = false;
		}
		else
		{
			std::cout << "," << std::endl;
		}
		std::cout << "{\"x\":" << node.id << ",\"y\":" << node.parent_id << ",\"z\":" << node.type + 2 << "}" << std::endl;
	}


	std::cout << "\n]\n} "<< std::endl;
}

void GaussianSeger::StaticParts(
	std::vector<PointCloud<PointT>::Ptr> SelectedClouds,
	std::vector<Eigen::Vector3f> SelectedTangent,
	std::vector<Eigen::Vector3f> SelectedPos,
	int StartIndex,
	int AllPartNum,
	double& length,
	double& total_torsion,
	double& total_area,
	double& total_width,
	int& CulledNum)
{
	if (AllPartNum == 0)
	{
		return;
	}

	Eigen::Vector3f LastPos;
	Eigen::Vector3f LastTangent;
	Eigen::Vector3f LastStart = { 0, 0, 0 };
	Eigen::Vector3f LastEnd = { 0, 0, 0 };
	Eigen::Vector3d LastNormal = { 1, 0, 0 };
	total_torsion = 0.;
	total_area = 0.;
	total_width = 0.;
	bool LastValid = false;
	CulledNum = 0;
	auto CurNodes = AllNodes.find(CurrentCloudId);
	if (CurNodes != AllNodes.end())
	{
		CurNodes->second.clear();
	}
	else
	{
		AllNodes[CurrentCloudId] = std::vector<Node>();
		CurNodes = AllNodes.find(CurrentCloudId);
	}
	CurNodes->second.reserve(std::min(AllPartNum,20));
	float node_step =std::min(AllPartNum,20)/ (float)AllPartNum;
	float node_t = 0;

	auto startPos = SelectedPos[StartIndex];
	int endIndex = (StartIndex - 1 + AllPartNum) % AllPartNum;
	auto endPos = SelectedPos[endIndex];
	if ((endPos - Stem_Center).norm() < (startPos - Stem_Center).norm())
	{
		startPos = endPos;
	}
	float min_dist = FLT_MAX;
	int start_id = 0;
	for(auto node: AllNodes["Stem.pcd"])
	{
		auto cur_dist =
			(startPos - node.pos).norm();
		if (cur_dist < min_dist)
		{
			min_dist = cur_dist;
			start_id = node.id;
		}
	}


	for (int i = 0; i < AllPartNum; i++)
	{
		int ind = (StartIndex + i) % AllPartNum;
		on_nurbs::NurbsDataCurve CurData;
		auto CurCloud = SelectedClouds[ind];
		auto CurPos = SelectedPos[ind];
		Eigen::Vector3d CurNormal;
		PointCloud2Vector2d(CurCloud, CurData.interior);
		Eigen::Vector3f CurNormalf;
		double min_v = 10000;
		double max_v = -10000;
		if (CurCloud->size() > 5)
		{
			Eigen::Vector3d mean;
			Eigen::Matrix3d eigenvectors;
			Eigen::Vector3d eigenvalues;
			on_nurbs::NurbsTools::pca(CurData.interior, mean, eigenvectors, eigenvalues);
			CurNormal = eigenvectors.col(0);
			CurNormal.normalize();
			CurNormalf = SelectedTangent[ind].cross(SelectedTangent[ind].cross(CurNormal.cast<float>())).normalized();
			CurNormalf = CurNormal.dot(LastNormal) <0 ? -CurNormalf : CurNormalf;
		}
		else
		{
			CulledNum++;
			continue;
			CurNormal = LastNormal;
		}

		for (std::size_t pIdx = 0; pIdx < CurCloud->size(); ++pIdx)
		{
			double v = ((*CurCloud)[pIdx].getVector3fMap() - CurPos).dot(CurNormalf);
			if (v > max_v)
				max_v = v;
			if (v < min_v)
				min_v = v;
		}
		std::ostringstream os;
		os << "width" << CurrentCloudId << "__" << i;
		total_width += max_v - min_v;
		auto cube_name = os.str();
		Eigen::Vector3f VectorEnd = CurPos + max_v * CurNormalf;
		Eigen::Vector3f VectorStart = CurPos + min_v * CurNormalf;
		viewer->addLine<PointCT>(
			PointCT(VectorStart.x(), VectorStart.y(), VectorStart.z(), 255, 255, 255),
			PointCT(VectorEnd.x(), VectorEnd.y(), VectorEnd.z(), 255, 255, 255),
			1.0,
			0.0,
			1.0,
			os.str());

		auto CurTangent = SelectedTangent[ind];
		if (i > 0 && LastValid)
		{
			auto CurLength = (LastPos - CurPos).norm();
			length += CurLength;
			total_torsion += 1 - abs(LastNormal.dot(CurNormal));
			double area = ((VectorEnd - VectorStart).cross(LastStart - VectorStart)).norm();
			if (!isnan(area))
				total_area += area;
			area = ((LastEnd - LastStart).cross(VectorEnd - LastStart)).norm();
			if (!isnan(area))
				total_area += area;
		}

		if(node_t + node_step > std::ceil(node_t))
		{
			Node node;
			node.pos = CurPos;
			node.tangent = CurTangent;
			node.normal =  CurNormalf;
			node.length = CurNodes->second.size()>0? (CurPos -CurNodes->second.at(CurNodes->second.size()-1).pos).norm() : 0;
			node.width = max_v - min_v;
			node.type = 1;
			node.generation = std::floor(node_t);
			node.parent_id = -1;
			CurNodes->second.push_back(node);
		}
		node_t += node_step;
		LastPos = CurPos;
		LastNormal = CurNormalf.cast<double>();
		LastTangent = CurTangent;
		LastStart = VectorStart;
		LastEnd = VectorEnd;
		LastValid = true;
	}
	if(CurNodes->second.size()>0)
	{
		CurNodes->second[0].parent_id = start_id;
		CurNodes->second[0].length = AllNodes["Stem.pcd"][start_id].length;
	}
}

void GaussianSeger::SliceParts(
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
	int& AllPartNum)
{
	int id = 0;
	double r, g, b;
	r = 0.7;
	g = 0.2;
	b = 0.0;
	int cp_red = CurrentCurve.Order() - 2;
	int resolution = 4;
	bool bLastValid = false;
	StartIndex = 0;
	curvatues = 0;
	start_t = 0.0f;
	end_t = 0.0f;
	start_knot = 0;
	int end_knot = 0;
	for (int i = 1; i < CurrentCurve.KnotCount() - 1 - cp_red; i++)
	{
		double dr = 1.0 / (resolution - 1);
		double xi0 = CurrentCurve.m_knot[i];
		double xid = (CurrentCurve.m_knot[i + 1] - xi0);
		double step = std::max(0.008, xid * dr);
		for (unsigned j = 0; j < resolution && (j * step) <= xid; j++, id++)
		{
			double t = xi0 + j * step;
			ON_3dPoint Point;
			ON_3dVector Tangent;
			CurrentCurve.EvTangent(t, Point, Tangent);

			ExtractIndices<PointT> extract_indices;
			Indices indices;
			Eigen::Translation3f Translation(Point.x, Point.y, Point.z);
			Eigen::Vector3f Scale = Eigen::Vector3f({ 1, 1, 0.01 });
			auto vec_T = Eigen::Vector3f(Tangent.x, Tangent.y, Tangent.z);
			vec_T.normalize();
			auto axis = UP.cross(vec_T);
			auto right = axis.cross(UP);
			double inter = UP.dot(vec_T);
			double dir_dot = right.dot(vec_T);
			auto rotation = Eigen::AngleAxisf(dir_dot > 0 ? acos(inter) : -acos(inter), axis);
			auto Quaternion = Eigen::Quaternionf(rotation);

			auto partclipper = CylinderClipper(SampleRadius * 0.2, 0.1, Translation.vector(), vec_T);
			partclipper.Clip3DPoints(downsampled_cloud, indices);
			if (indices.size() > 0)
			{
				// VisualTagent(cloud_id, id, Point, Tangent);

				r = (rand() % 65536) / 140000.f + 0.2f;
				g = (rand() % 65536) / 140000.f + 0.2f;
				b = (rand() % 65536) / 140000.f + 0.2f;

				PointCloud<PointT>::Ptr cloud_out(new PointCloud<PointT>);
				extract_indices.setInputCloud(downsampled_cloud);
				extract_indices.setIndices(std::make_shared<Indices>(indices));
				extract_indices.filter(*cloud_out);
				std::ostringstream os;
				if (bShowCylinder && false)
				{
					os << CurrentCloudId << "_Cube_" << id;
					auto cube_name = os.str();
					viewer->addCube(Translation.vector(), Quaternion, 1, 1, SampleRadius * 0.2, cube_name);
					viewer->setShapeRenderingProperties(
						visualization::PCL_VISUALIZER_COLOR,
						r,
						g,
						b,
						cube_name);
					viewer->setShapeRenderingProperties(
						visualization::PCL_VISUALIZER_OPACITY,
						0.5,
						cube_name);
					viewer->setShapeRenderingProperties(
						visualization::PCL_VISUALIZER_REPRESENTATION,
						visualization::PCL_VISUALIZER_REPRESENTATION_WIREFRAME,
						cube_name);
					os.clear();
				}
				os << "_Cloud_" << id;
				//auto cloud_name = os.str();
				//std::cout << cloud_name << "  size ::" << cloud_out->size() << std::endl;
				//viewer->addPointCloud(cloud_out, cloud_name);
				//viewer->setPointCloudRenderingProperties(
				//	visualization::PCL_VISUALIZER_COLOR,
				//	r,
				//	g,
				//	b,
				//	cloud_name);
				SelectedClouds.push_back(cloud_out);
				SelectedTangent.push_back(vec_T);
				SelectedPos.push_back(Translation.vector());
				auto curvature = CurrentCurve.CurvatureAt(t);
				curvatues += curvature.Length();
				SelectedCurvature.push_back(
					Eigen::Vector3f(curvature.x, curvature.y, curvature.z));
				if (!bLastValid)
				{
					bLastValid = true;
					StartIndex = SelectedClouds.size() - 1;
					start_t = t;
					start_knot = i;
				}
			}
			else
			{
				if (bLastValid)
				{
					end_knot = i;
					end_t = t;
				}
				bLastValid = false;
			}
		}
		viewer->spinOnce(1);
	}

	AllPartNum = SelectedClouds.size();
}

PointCloud<PointT>::Ptr GaussianSeger::DownSample(
	PointCloud<PointT>::Ptr cloud,
	const float leaf_size)
{
	VoxelGrid<PointT> voxel_grid;
	voxel_grid.setLeafSize(leaf_size, leaf_size, leaf_size);

	voxel_grid.setInputCloud(cloud);
	PointCloud<PointT>::Ptr downsampled_cloud(cloud);
	voxel_grid.filter(*downsampled_cloud);

	viewer->addPointCloud<PointT>(downsampled_cloud, "down cloud");
	return downsampled_cloud;
}

ON_NurbsCurve GaussianSeger::SimpleFitCurve(
	PointCloud<PointT>::Ptr cloud,
	unsigned n_control_points,
	on_nurbs::vector_vec3d RefPoints)

{
	// #################### CURVE PARAMETERS #########################
	unsigned order(3);
	on_nurbs::FittingCurve::Parameter curve_params;
	curve_params.smoothness = 0.000001;

	on_nurbs::NurbsDataCurve data;
	PointCloud2Vector2d(cloud, data.interior);

	// #################### CURVE FITTING #########################
	ON_NurbsCurve curve = RefPoints.size() < 3
		? initNurbsCurvePCAOffset(order, data.interior, n_control_points)
		: initNurbsCurveRef(order, data.interior, n_control_points, RefPoints);

	on_nurbs::FittingCurve fit(&data, curve);
	fit.assemble(curve_params);
	fit.solve();

	return fit.m_nurbs;
}

std::string GaussianSeger::CalculatePart(
	PointCloud<PointT>::Ptr cloud,
	std::string cloud_id,
	bool ShouldDownSample)
{
	CurrentCloudId = cloud_id;
	PointCloud<PointT>::Ptr InputCloud = cloud;
	if (ShouldDownSample)
		InputCloud = DownSample(cloud);
	on_nurbs::vector_vec3d posss;

	// statistic ellapsed time
	 auto StartTime = std::chrono::high_resolution_clock::now();

	IterSample(cloud, posss, this);
	if (posss.size() <= 6)
	{
		IterSample(cloud, posss, this, SampleRadius * 0.5);
	}
	auto EndTime = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(EndTime - StartTime);
	std::cout << "ex depart Time: " << duration.count() << " ms" << std::endl;
	StartTime = EndTime;
	auto CurrentCurve = SimpleFitCurve(InputCloud, 20, posss);
	EndTime = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(EndTime - StartTime);
	std::cout << "Curve Fitting Time: " << duration.count() << " ms" << std::endl;


	if (bShowCurve)
	{
		PointCloud<PointCT>::Ptr CurveCloud(new PointCloud<PointCT>);
		on_nurbs::Triangulation::convertCurve2PointCloud(CurrentCurve, CurveCloud, 8);
		// io::savePCDFile("d:/GS/ExStem/stem_curve.pcd", *CurveCloud, true);
		VisualizeCurve(viewer, cloud_id, CurveCloud, 1.0, 0.0, 0.0, false);
	}
	auto UP = Eigen::Vector3f(0, 0, 1);

	std::vector<PointCloud<PointT>::Ptr> SelectedClouds;
	std::vector<Eigen::Vector3f> SelectedTangent;
	std::vector<Eigen::Vector3f> SelectedPos;
	std::vector<Eigen::Vector3f> SelectedCurvature;
	int StartIndex;
	float curvatues;
	float start_t;
	float end_t;
	int start_knot;
	int AllPartNum;
	SliceParts(InputCloud,
		CurrentCurve,
		UP,
		SelectedClouds,
		SelectedTangent,
		SelectedPos,
		SelectedCurvature,
		StartIndex,
		curvatues,
		start_t,
		end_t,
		start_knot,
		AllPartNum);
	double length = 0;
	double LiqAngle;

	bool bReverse = false;

	{
		auto startTangent = SelectedTangent[StartIndex];
		auto startPos = SelectedPos[StartIndex];
		int endIndex = (StartIndex - 1 + AllPartNum) % AllPartNum;
		auto endPos = SelectedPos[endIndex];
		if ((endPos - Stem_Center).norm() < (startPos - Stem_Center).norm())
		{
			startPos = endPos;
			startTangent = SelectedTangent[endIndex];
			bReverse = true;
		}

		// Stem_Curve.
		double firstKnot = Stem_Curve.m_knot[1];
		double lastKnot = Stem_Curve.m_knot[Stem_Curve.KnotCount() - 2];
		int closestKnotId = 0;
		double min_dist = 1000;
		double KnotStep = (lastKnot - firstKnot);
		for (int i = 0; i < 20; i++)
		{
			ON_3dPoint curPoint;
			Stem_Curve.EvPoint(i * 0.05 * KnotStep + firstKnot, curPoint);
			auto cur_dist =
				(startPos - Eigen::Vector3f(curPoint.x, curPoint.y, curPoint.z)).norm();
			if (cur_dist < min_dist)
			{
				min_dist = cur_dist;
				closestKnotId = i;
				UpdateSphere(curPoint);
			}
		}
		firstKnot = Stem_Curve.m_knot[closestKnotId  >= 2 ? closestKnotId - 1 : 1];
		lastKnot = Stem_Curve.m_knot[closestKnotId  <= Stem_Curve.KnotCount() -3
			? closestKnotId + 1
			: Stem_Curve.KnotCount() - 2];
		KnotStep = (lastKnot - firstKnot);
		min_dist = 1000;
		double closest_t = 0;
		for (int i = 0; i < 20; i++)
		{
			ON_3dPoint curPoint;
			double cur_t = i * 0.05 * KnotStep + firstKnot;
			Stem_Curve.EvPoint(cur_t, curPoint);
			auto cur_dist =
				(startPos - Eigen::Vector3f(curPoint.x, curPoint.y, curPoint.z)).norm();
			if (cur_dist < min_dist)
			{
				min_dist = cur_dist;
				closest_t = cur_t;
				UpdateSphere(curPoint);
			}
		}
		ON_3dPoint Point;
		ON_3dVector ClosetTangent_on_Stem;
		Stem_Curve.EvTangent(closest_t, Point, ClosetTangent_on_Stem);
		Eigen::Vector3f closestTangent(
			ClosetTangent_on_Stem.x,
			ClosetTangent_on_Stem.y,
			ClosetTangent_on_Stem.z);
		UpdateSphere(Point);
		LiqAngle = 180.f * acos(abs(closestTangent.dot(startTangent))) / M_PI;
	}

	double total_torsion;
	double total_area;
	double total_width;
	int CulledNum;
	StaticParts(
		SelectedClouds,
		SelectedTangent,
		SelectedPos,
		StartIndex,
		AllPartNum,
		length,
		total_torsion,
		total_area,
		total_width,
		CulledNum);

	std::ostringstream static_info;
	static_info << "total length		  ::  " << length << "\n";
	AllPartNum -= CulledNum;
	{
		static_info << "average curvature	::  " << curvatues / AllPartNum << "\n";
		static_info << "total torsion		  ::  " << total_torsion << "\n";
		static_info << "average torsion		::  " << total_torsion / AllPartNum << "\n";
		static_info << "total area			  ::  " << total_area / 2.0 << "\n";
		static_info << "average width		  ::  " << total_width / AllPartNum << "\n";
		static_info << "Leaf inclination  ::  " << LiqAngle << "\n";
	}

	viewer->removePointCloud("down cloud");
	viewer->updateText(static_info.str(), 0, 300, "static_info");
	std::cout << static_info.str() << std::endl;


	return static_info.str();
}

Eigen::Vector3d GetMean(
	PointCloud<PointT>::Ptr InCloud)
{
	on_nurbs::NurbsDataCurve Data;
	PointCloud2Vector2d(InCloud, Data.interior);
	return on_nurbs::NurbsTools::computeMean(Data.interior);
}

// 拆分采样结果

int SplitClusters(
	PointCloud<PointT>::Ptr InputCloud,
	float Tolerance,
	int MinSize,
	int Maxsize,
	Eigen::Vector3d InputPos,
	Eigen::Vector3d& OutPos_1,
	Eigen::Vector3d& OutPos_2
)
{
	if (MinSize > InputCloud->size())
		return -1;
	search::KdTree<PointT>::Ptr tree(new search::KdTree<PointT>());
	tree->setInputCloud(InputCloud);

	std::vector<pcl::PointIndices> cluster_indices;
	pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
	ec.setClusterTolerance(Tolerance);
	ec.setMinClusterSize(MinSize);
	ec.setMaxClusterSize(Maxsize);
	ec.setSearchMethod(tree);
	ec.setInputCloud(InputCloud);
	ec.extract(cluster_indices);

	if (cluster_indices.size() <= 0)
	{
		std::cout << "Failed to split cluster" << std::endl;
		SplitClusters(InputCloud, Tolerance * 2, MinSize, Maxsize * 2, InputPos, OutPos_1, OutPos_2);
		return 0;
	}
	else if (cluster_indices.size() >= 4)
	{
		return SplitClusters(InputCloud, Tolerance * 2, MinSize, Maxsize * 2, InputPos, OutPos_1, OutPos_2);
	}

	std::vector<PointCloud<PointT>::Ptr> SplitClouds;
	for (const auto& cluster : cluster_indices)
	{
		PointCloud<PointT>::Ptr SplitCloud(new PointCloud<PointT>);
		pcl::copyPointCloud(*InputCloud, cluster, *SplitCloud);
		SplitClouds.push_back(SplitCloud);
	}
	OutPos_1 = Eigen::Vector3d(0, 0, 0);
	OutPos_2 = Eigen::Vector3d(0, 0, 0);
	if (cluster_indices.size() == 1)
	{
		OutPos_1 = GetMean(SplitClouds[0]);
		return 1;
	}
	if (cluster_indices.size() == 2)
	{
		OutPos_1 = GetMean(SplitClouds[0]);
		OutPos_2 = GetMean(SplitClouds[1]);
	}
	if (cluster_indices.size() == 3)
	{
		if (SplitClouds[0]->size() < SplitClouds[1]->size() / 4 && SplitClouds[0]->size() < SplitClouds[2]->size() / 4)
		{
			OutPos_1 = GetMean(SplitClouds[1]);
			OutPos_2 = GetMean(SplitClouds[2]);
		}
		else if (SplitClouds[1]->size() < SplitClouds[0]->size() / 4 && SplitClouds[1]->size() < SplitClouds[2]->size()
			/ 4)
		{
			OutPos_1 = GetMean(SplitClouds[0]);
			OutPos_2 = GetMean(SplitClouds[2]);
		}
		else if (SplitClouds[2]->size() < SplitClouds[0]->size() / 4 && SplitClouds[2]->size() < SplitClouds[1]->size()
			/ 4)
		{
			OutPos_1 = GetMean(SplitClouds[0]);
			OutPos_2 = GetMean(SplitClouds[1]);
		}
		else
		{
			Eigen::Vector3d dir_1 = GetMean(SplitClouds[0]) - InputPos;
			Eigen::Vector3d dir_2 = GetMean(SplitClouds[1]) - InputPos;
			Eigen::Vector3d dir_3 = GetMean(SplitClouds[2]) - InputPos;
			dir_1.normalize();
			dir_2.normalize();
			dir_3.normalize();
			auto angle_1 = dir_1.dot(dir_2);
			auto angle_2 = dir_2.dot(dir_3);
			auto angle_3 = dir_3.dot(dir_1);

			if (angle_1 >= angle_2 && angle_1 >= angle_3)
			{
				OutPos_1 = GetMean(SplitClouds[1]);
				OutPos_2 = GetMean(SplitClouds[2]);
			}
			else if (angle_2 >= angle_3 && angle_2 >= angle_1)
			{
				OutPos_1 = GetMean(SplitClouds[0]);
				OutPos_2 = GetMean(SplitClouds[2]);
			}
			else if (angle_3 >= angle_1 && angle_3 >= angle_2)
			{
				OutPos_1 = GetMean(SplitClouds[0]);
				OutPos_2 = GetMean(SplitClouds[1]);
			}
		}
	}
	return 2;
}


void ExtractBoth(
	PointCloud<PointT>::Ptr Input,
	Indices InIndices,
	PointCloud<PointT>::Ptr InRange,
	PointCloud<PointT>::Ptr OutRange)
{
	pcl::copyPointCloud(*Input, InIndices, *InRange);
	if (OutRange.get())
	{
		ExtractIndices<PointT> extract_indices;
		extract_indices.setInputCloud(Input);
		extract_indices.setIndices(std::make_shared<Indices>(InIndices));
		extract_indices.setNegative(true);
		extract_indices.filter(*OutRange);
	}
}

void IterSampleFromPoint(
	PointCloud<PointT>::Ptr InCloud,
	Eigen::Vector3d init_pos,
	Eigen::Vector3d init_dir,
	float Radius,
	float Tolerance,
	on_nurbs::vector_vec3d& OutPoints,
	GaussianSeger* viewer
)
{
	Eigen::Vector3d CurPos = init_pos;
	Eigen::Vector3d CurDir = init_dir;
	Eigen::Vector3d Pos_1, Pos_2;
	bool valid = false;
	do
	{
		SphereClipper Clipper(Radius, CurPos.cast<float>());
		Indices IndicesInRange;
		Clipper.Clip3DPoints(InCloud, IndicesInRange);
		PointCloud<PointT>::Ptr CloudInRange(new PointCloud<PointT>);

		ExtractBoth(InCloud, IndicesInRange, CloudInRange, InCloud);
		int SplitNum = SplitClusters(CloudInRange, Tolerance, 10, INT_MAX, init_pos, Pos_1, Pos_2);
		if (SplitNum < 0)
		{
			Clipper.radius = Radius * 2;
			Clipper.Clip3DPoints(InCloud, IndicesInRange);
			if (IndicesInRange.size() < 10)
				return;
			else
			{
				Radius *= 2;
				ExtractBoth(InCloud, IndicesInRange, CloudInRange, InCloud);
				SplitNum = SplitClusters(CloudInRange, Tolerance, 10, INT_MAX, init_pos, Pos_1, Pos_2);
				if (SplitNum < 0)
					return;
			}
		}
		if (SplitNum <= 1)
		{
			if (SplitNum == 0)
				Tolerance *= 2;
			CurPos = Pos_1;
			CurDir = (CurPos - init_pos).normalized();

			if (CurDir.dot(init_dir) > 0.73)
			{
				OutPoints.push_back(CurPos);
				init_dir = CurDir;
				init_pos = CurPos;
			}
			else
			{
				return;
			}
		}
		if (SplitNum == 2)
		{
			if ((Pos_1 - init_pos).normalized().dot((Pos_2 - init_pos).normalized()) > 0.5)
			{
				CurPos = (Pos_1 + Pos_2) * 0.5;
				CurDir = (CurPos - init_pos).normalized();
				OutPoints.push_back(CurPos);
				init_dir = CurDir;
				init_pos = CurPos;
			}
			else
			{
				return;
			}
		}
		ON_3dPoint curPoint(CurPos[0], CurPos[1], CurPos[2]);
		viewer->UpdateSphere(curPoint, Radius, "IterRange");
	} while (true);
}

void IterSample(
	PointCloud<PointT>::Ptr InCloud,
	on_nurbs::vector_vec3d& OutPoints,
	GaussianSeger* viewer,
	float Radius
)
{
	PointCloud<PointT>::Ptr ClonedCloud(new PointCloud<PointT>);
	ClonedCloud->resize(InCloud->size());
	pcl::detail::copyPointCloudMemcpy(*InCloud, *ClonedCloud);

	int rand_idx = rand() % ClonedCloud->size();
	Eigen::Vector3d init_pos = (*ClonedCloud)[rand_idx].getVector3fMap().cast<double>();

	//initialize sample radius
	float tolerance = 0.05;
	Eigen::Vector3d Pos_1, Pos_2;
	bool valid = false;
	bool RigonValid = false;
	PointCloud<PointT>::Ptr CloudExtracted(new PointCloud<PointT>);

	do
	{
		Radius *= 2;
		if (Radius > 1000)
			return;

		SphereClipper Clipper(Radius, init_pos.cast<float>());
		Indices IndicesInRange;
		Clipper.Clip3DPoints(ClonedCloud, IndicesInRange);
		PointCloud<PointT>::Ptr CloudInRange(new PointCloud<PointT>);

		ExtractBoth(ClonedCloud, IndicesInRange, CloudInRange, ClonedCloud);
		int SplitNum = SplitClusters(CloudInRange, tolerance, 10, INT_MAX, init_pos, Pos_1, Pos_2);
		if (SplitNum == 0)
			tolerance *= 2;
		valid = SplitNum >= 2;
		if (SplitNum == 1)
		{
			pcl::concatenate(*CloudExtracted, *CloudInRange, *CloudExtracted);
			on_nurbs::vector_vec3d CurData;
			PointCloud2Vector2d(CloudExtracted, CurData);
			Eigen::Vector3d mean;
			Eigen::Matrix3d eigenvectors;
			Eigen::Vector3d eigenvalues;
			on_nurbs::NurbsTools::pca(
				CurData,
				mean,
				eigenvectors,
				eigenvalues);
			Eigen::Vector3d nor_dir = eigenvectors.col(0);
			valid = (Pos_1 - init_pos).normalized().dot(nor_dir) > std::cos(3.14159 / 6.0);
		}
	} while (!valid);
	viewer->SampleRadius = Radius;
	on_nurbs::vector_vec3d SampledPoints = { init_pos };

	if (Pos_1.norm() > 0.000001)
	{
		IterSampleFromPoint(ClonedCloud,
			Pos_1,
			(Pos_1 - init_pos).normalized(),
			Radius,
			tolerance,
			SampledPoints,
			viewer);
	}
	if (Pos_2.norm() > 0.000001)
	{
		on_nurbs::vector_vec3d ReversedSampledPoints;
		IterSampleFromPoint(ClonedCloud,
			Pos_2,
			(Pos_2 - init_pos).normalized(),
			Radius,
			tolerance,
			ReversedSampledPoints,
			viewer);

		for (int i = SampledPoints.size() - 1; i >= 0; --i)
		{
			OutPoints.push_back(SampledPoints[i]);
		}

		OutPoints.insert(OutPoints.end(), ReversedSampledPoints.begin(), ReversedSampledPoints.end());
	}
	else
	{
		OutPoints = SampledPoints;
	}
	if (viewer->bShowCylinder&&false)
		for (int i = 0; i < OutPoints.size(); i++)
		{
			ON_3dPoint curPoint(OutPoints[i][0], OutPoints[i][1], OutPoints[i][2]);
			char name[20];
			sprintf(name, "IterRange%d", i);
			std::cout << name << std::endl;
			viewer->UpdateSphere(curPoint, 0.05, name);
		}
	viewer->viewer->removeShape("IterRange");
}

void ExtractLocalProperty(
	Eigen::Vector3d last_mean,
	float& total_offset,
	Indices indices,
	PointCloud<PointT>::Ptr InCloud,
	on_nurbs::vector_vec3d CurData,
	Eigen::Vector3d& mean,
	Eigen::Vector3d tangent,
	Eigen::Matrix3d& eigenvectors,
	Eigen::Vector3d& eigenvalues,
	double& max_v,
	double& min_v,
	double& max_v2,
	double& min_v2)
{
	auto curve_mean = mean.cast<float>();
	on_nurbs::NurbsTools::pca(
		CurData,
		mean,
		eigenvectors,
		eigenvalues);

	max_v = FLT_MIN;
	min_v = FLT_MAX;
	max_v2 = FLT_MIN;
	min_v2 = FLT_MAX;
	Eigen::Vector3d nor_dir = eigenvectors.col(0).normalized();
	nor_dir = tangent.cross(tangent.cross(nor_dir).normalized()).normalized();
	Eigen::Vector3f nor_dirf = nor_dir.cast<float>();
	for (std::size_t pIdx = 0; pIdx < InCloud->size(); ++pIdx)
	{
		double v = ((*InCloud)[pIdx].getVector3fMap()-curve_mean).dot(nor_dirf);
		if (v > max_v)
			max_v = v;
		if (v < min_v)
			min_v = v;
	}

	total_offset += (mean - last_mean).norm();
	eigenvalues /= indices.size();
}

std::string GaussianSeger::CalculateStem(
	PointCloud<PointT>::Ptr cloud,
	std::string cloud_id,
	bool ShouldDownSample)
{
	PointCloud<PointT>::Ptr InputCloud = cloud;
	if (ShouldDownSample)
		InputCloud = DownSample(cloud);
	std::cout<<"Stem cloud size"<<InputCloud->size()<<std::endl;
	auto StartTime  = std::chrono::high_resolution_clock::now();

	on_nurbs::vector_vec3d posss;
	IterSample(cloud, posss, this);
	auto EndTime = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(EndTime - StartTime);
	std::cout << "ex depart Time: " << duration.count() << " ms" << std::endl;
	StartTime = EndTime;
	auto CurrentCurve = SimpleFitCurve(InputCloud, 80, posss);
	EndTime = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(EndTime - StartTime);
	std::cout << "Curve Fitting Time: " << duration.count() << " ms" << std::endl;
	Stem_Curve = CurrentCurve;
	on_nurbs::vector_vec3d data;
	PointCloud2Vector2d(InputCloud, data);
	auto mean = on_nurbs::NurbsTools::computeMean(data);
	Stem_Center = Eigen::Vector3f(mean.x(), mean.y(), mean.z());

	if (bShowCurve)
	{
		PointCloud<PointCT>::Ptr CurveCloud(new PointCloud<PointCT>);
		on_nurbs::Triangulation::convertCurve2PointCloud(CurrentCurve, CurveCloud, 8);
		// io::savePCDFile("d:/GS/ExStem/stem_curve.pcd", *CurveCloud, true);
		VisualizeCurve(viewer, cloud_id, CurveCloud, 1.0, 0.0, 0.0, false);
	}
	auto UP = Eigen::Vector3f(0, 0, 1);

	std::vector<PointCloud<PointT>::Ptr> SelectedClouds;
	std::vector<Eigen::Vector3f> SelectedTangent;
	std::vector<Eigen::Vector3f> SelectedPos;
	std::vector<Eigen::Vector3f> SelectedCurvature;
	int StartIndex;
	float curvatues;
	float start_t;
	float end_t;
	int start_knot;
	int AllPartNum;
	SliceParts(InputCloud,
		CurrentCurve,
		UP,
		SelectedClouds,
		SelectedTangent,
		SelectedPos,
		SelectedCurvature,
		StartIndex,
		curvatues,
		start_t,
		end_t,
		start_knot,
		AllPartNum);
	double length = 0;
	double LiqAngle;

	bool bReverse = false;
	AllNodes.clear();
	AllNodes["Stem.pcd"] = std::vector<Node>();
	auto CurNodes = AllNodes.find("Stem.pcd");
	{
		int cur_not = start_knot;
		float total_t = CurrentCurve.m_knot[CurrentCurve.KnotCount() - 1];
		if (end_t < start_t)
			end_t += total_t;
		float step = (end_t - start_t) / 35.f;
		Eigen::Vector3d last_mean = Eigen::Vector3d(SelectedPos[StartIndex].x(),
			SelectedPos[StartIndex].y(),
			SelectedPos[StartIndex].z());
		Eigen::Vector3f last_nor = Eigen::Vector3f(0,1,0);
		float total_offset = 0;
		int cut_idx = 0;
		std::vector<float> R_15;
		std::vector<float> N_15;
		R_15.reserve(15);
		N_15.reserve(15);
		int generation = 0;
		for (float t = start_t; t <= total_t; t += step,++generation)
		{
			ON_3dPoint Point;
			ON_3dVector Tangent;
			float cur_t = fmod(t, total_t);
			CurrentCurve.EvTangent(cur_t, Point, Tangent);

			Indices indices;
			Eigen::Translation3f Translation(Point.x, Point.y, Point.z);
			auto vec_T = Eigen::Vector3f(Tangent.x, Tangent.y, Tangent.z);
			vec_T.normalize();

			auto partclipper = CylinderClipper(SampleRadius * 0.15, 0.5, Translation.vector(),vec_T);
			partclipper.Clip3DPoints(InputCloud, indices);

			if (indices.size() > 5)
			{
				ExtractIndices<PointT> extract_indices;
				PointCloud<PointT>::Ptr cloud_out(new PointCloud<PointT>);
				extract_indices.setInputCloud(InputCloud);
				extract_indices.setIndices(std::make_shared<Indices>(indices));
				extract_indices.filter(*cloud_out);
				on_nurbs::vector_vec3d CurData;
				PointCloud2Vector2d(cloud_out, CurData);
				Eigen::Vector3d mean {Translation.x(),Translation.y(),Translation.z()};
				Eigen::Matrix3d eigenvectors;
				Eigen::Vector3d eigenvalues;
				double max_v;
				double min_v;
				double max_v2;
				double min_v2;
				ExtractLocalProperty(last_mean,
					total_offset,
					indices,
					cloud_out,
					CurData,
					mean,
					vec_T.cast<double>(),
					eigenvectors,
					eigenvalues,
					max_v,
					min_v,
					max_v2,
					min_v2);

				Eigen::Vector3d VectorStart =  Eigen::Vector3d{Translation.x(),Translation.y(),Translation.z()}+ eigenvectors.col(0)*max_v;
				Eigen::Vector3d VectorEnd =  Eigen::Vector3d{Translation.x(),Translation.y(),Translation.z()}+ eigenvectors.col(0)*min_v;
				std::ostringstream os;
				//os << "width" << "Stem"<< "__" << t;
				//viewer->addLine<PointCT>(
				//PointCT(VectorStart.x(), VectorStart.y(), VectorStart.z(), 255, 255, 255),
				//PointCT(VectorEnd.x(), VectorEnd.y(), VectorEnd.z(), 255, 255, 255),
				//1.0,
				//0.0,
				//1.0,
				//os.str()
				//);
				//VectorStart =  Eigen::Vector3d{Translation.x(),Translation.y(),Translation.z()}+ eigenvectors.col(0);
				//VectorEnd =  Eigen::Vector3d{Translation.x(),Translation.y(),Translation.z()};
				os << "widths" << "Stem"<< "__" << t;
				viewer->addLine<PointCT>(
				PointCT(VectorStart.x(), VectorStart.y(), VectorStart.z(), 255, 255, 255),
				PointCT(VectorEnd.x(), VectorEnd.y(), VectorEnd.z(), 255, 255, 255),
				0.0,
				1.0,
				1.0,
				os.str()
				);
				/*std::cout << total_offset << "," << sqrt(eigenvalues[0]) << ","
					<< sqrt(eigenvalues[1]) << ",";
				std::cout << eigenvectors.col(2)[0] << "," << eigenvectors.col(2)[1] << ","
					<< eigenvectors.col(2)[2] << "," << indices.size() << ","
					<< CurrentCurve.CurvatureAt(cur_t).Length() << "," << max_v - min_v
					<< "," << max_v2 - min_v2 << std::endl;*/
				std::cout<<max_v-min_v<<std::endl;
				Eigen::Vector3f nor = eigenvectors.col(0).normalized().cast<float>();
				CurNodes->second.push_back(
						Node(
						Translation.cast<float>().vector(),
						vec_T,
						nor.dot(last_nor)<0? -nor:nor,
						(max_v-min_v)/2,
						(last_mean-mean).norm(),
						generation,
						35,
						CurNodes->second.size(),
						CurNodes->second.size()-1,
						0)
					);
				last_mean = mean;
				last_nor = nor.dot(last_nor)<0? -nor:nor;
			}
		}
				length = total_offset;
	}

	//double total_torsion;
	//double total_area;
	//double total_width;
	//int CulledNum;
	//StaticParts(
	//	SelectedClouds,
	//	SelectedTangent,
	//	SelectedPos,
	//	StartIndex,
	//	AllPartNum,
	//	length,
	//	total_torsion,
	//	total_area,
	//	total_width,
	//	CulledNum);

	std::ostringstream static_info;
	static_info << "total length		  ::  " << length << "\n";
	{
		static_info << "Leaf Num          ::  " << LeafNum << "\n";
	}
	viewer->removePointCloud("down cloud");
	viewer->updateText(static_info.str(), 0, 300, "static_info");
	std::cout << static_info.str() << std::endl;


	return static_info.str();
}
