#include <list>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <igraph/igraph.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkBitArray.h>
#include <vtkKdTreePointLocator.h>
#include <vtkSplineFilter.h>
#include <vtkParametricSpline.h>
#include <vtkKochanekSpline.h>
#include <vtkParametricFunctionSource.h>

#define _min(a,b) ((a<b)?a:b)
#define _max(a,b) ((a>b)?a:b)
#define _eps 1E-6

	class _points {
	private:
		double Bounds[6];
		vtkSmartPointer<vtkPoints> Points;
		vtkSmartPointer<vtkPolyData> PolyPoints;
		vtkSmartPointer<vtkKdTreePointLocator> Tree;

	public:
		std::vector<double> D1;
		std::vector<double> D2;
		std::vector<double> D12;
		std::vector<double> Area;
		_points(vtkPolyData *Skeleton);
		void CalculateRelativeDistances();
		double GetX(int i);
		double GetY(int i);
		double GetZ(int i);
		double GetBoundXmin() {return Bounds[0];}
		double GetBoundXmax() {return Bounds[1];}
		double GetBoundYmin() {return Bounds[2];}
		double GetBoundYmax() {return Bounds[3];}
		double GetBoundZmin() {return Bounds[4];}
		double GetBoundZmax() {return Bounds[5];}
	};

	class _curve {
	private:
		vtkIdType n;
		int dsl, dsr;
		double _spacing, R2;
		vtkSmartPointer<vtkPolyData> Spline;
		vtkSmartPointer<vtkFloatArray> Tx;
		vtkSmartPointer<vtkFloatArray> Ty;
		vtkSmartPointer<vtkFloatArray> Tz;
		vtkSmartPointer<vtkFloatArray> Kr;
		void SetTangent(vtkIdType id, double t[3]);
		void SetTangent(vtkIdType id, double tx, double ty, double tz);
		void GetTangent(vtkIdType id, double t[3]);
	public:
		std::vector<double> TCorr;
		_curve(vtkPolyData *Skeleton, vtkIdType cell_id, double spacing);
		~_curve() {
			TCorr.clear();
		}
		vtkIdType GetNumberOfPoints() {return n;}
		void CalculateTangentVersors(int order);
		void CalculateAverageTangentCorrelation();
		void CalculateEndFluctuations();
		void CalculateCurvature();
		double GetSpacing() {return _spacing;}
		double GetLength() {return n * _spacing;}
		double GetR2() {return R2;}
		void GetSplinePolyData(vtkPolyData *Line);
		void ExportCurvatures(FILE *f);

	};

	void GeometricalAnalysis(const char FilePrefix[]);
	void ExportResults(std::vector<_curve> Curves, _points Nodes, const char FilePrefix[]);
	double DotProduct(double to[3], double tf[3]);
	double GetLength(vtkPolyData *Skeleton, vtkIdType cell_id, bool plot);

/* ====================================================================
	AUXILIAR FUNCTIONS
======================================================================*/

double GetLength(vtkPolyData *Skeleton, vtkIdType cell_id, bool plot) {
	double ds = 0.0, d2 = 0.0;
	double r[3], u[3], d, length = 0.0;
	vtkCell *Line = Skeleton -> GetCell(cell_id);
	for (int i = 1; i < Line -> GetNumberOfPoints(); i++) {
		Skeleton -> GetPoints() -> GetPoint(Line->GetPointId(i-1),u);
		Skeleton -> GetPoints() -> GetPoint(Line->GetPointId(i-0),r);
		d = sqrt(pow(r[0]-u[0],2)+pow(r[1]-u[1],2)+pow(r[2]-u[2],2));
		length += d;
		ds += d;
		d2 += d*d;
	}
	if (plot) {
		ds /= Line -> GetNumberOfPoints()-1;
		d2 /= Line -> GetNumberOfPoints()-1;
		d2 = sqrt(d2-ds*ds);
		printf("\t%f\t+/-%f\n",ds,d2);
	}
	return length;
}

double DotProduct(double to[3], double tf[3]) {
	double dotp = to[0]*tf[0] + to[1]*tf[1] + to[2]*tf[2];
	dotp = (dotp> 1.00) ? 1.00 : dotp;
	dotp = (dotp<-1.00) ?-1.00 : dotp;
	return dotp;
}

void GeometricalAnalysis(const char FilePrefix[]) {
	
	/* PARAMETRIC CURVES */

	char _FullPath[256];
	sprintf(_FullPath,"%s_skeleton.vtk",FilePrefix);

	vtkSmartPointer<vtkPolyDataReader> PolyReader = vtkSmartPointer<vtkPolyDataReader>::New();
	PolyReader -> SetFileName(_FullPath);
	PolyReader -> Update();

	vtkPolyData *Skell = PolyReader -> GetOutput();

	int NLines = Skell -> GetNumberOfCells();

	std::vector<_curve> Curves;

	int line;
	double edge_length;
	for (line = 0; line < NLines; line++) {
		edge_length = GetLength(Skell,line,false);
		if ( edge_length > 1.5 ) {
			#ifdef DEBUG
				printf("Edge: %d (%1.4fum)\n",line,edge_length);
			#endif
			_curve curve(Skell,line,0.1);
			curve.CalculateTangentVersors(1);
			curve.CalculateAverageTangentCorrelation();
			curve.CalculateCurvature();
			curve.CalculateEndFluctuations();
			Curves.push_back(curve);
		}
	}

	/* POINTS IN SPACE */
	
	_points Nodes(Skell);
	Nodes.CalculateRelativeDistances();

	/* EXPORTING RESULTS */

	if (Curves.size()) {
		ExportResults(Curves,Nodes,FilePrefix);
	} else {
		printf("No edges long enough found.\n");
	}

	Curves.clear();
}

void ExportResults(std::vector<_curve> Curves, _points Nodes, const char FilePrefix[]) {

	/* PARAMETRIC CURVES */

	#ifdef DEBUG
		printf("Calculating tangent correlation decay...\n");
	#endif

	int p, c;
	double s = Curves[0].GetSpacing();

	#ifdef DEBUG
		printf("\t#Curves: %d\n",(int)Curves.size());
	#endif

	int max_length = 0;
	for (c = 0; c < Curves.size(); c++) {
		max_length = (Curves[c].TCorr.size()>max_length) ? Curves[c].TCorr.size() : max_length;
	}

	int *Count = new int[max_length];
	double *TCorr = new double[max_length];

	for (p = 0; p < max_length; p++) {
		Count[p] = 0;
		TCorr[p] = 0.0;
	}

	for (c = 0; c < Curves.size(); c++) {
		for (p = 0; p < Curves[c].TCorr.size(); p++) {
			Count[p] ++;
			TCorr[p] += Curves[c].TCorr[p];
		}
	}

	char _FullPath[256];
	sprintf(_FullPath,"%s.tcorr",FilePrefix);
	FILE *f = fopen(_FullPath,"w");
	fprintf(f,"Vectors Spacing (um)\tTangent Correlation\n");
	fprintf(f,"%1.4f\t%1.6f\n",0.0,1.0);
	for (p = 0; p < max_length; p++) {
		fprintf(f,"%1.4f\t%1.6f\n",p*s,TCorr[p]/Count[p]);
	}
	fclose(f);
	delete[] TCorr;
	delete[] Count;

	vtkSmartPointer<vtkAppendPolyData> Append = vtkSmartPointer<vtkAppendPolyData>::New();

	for (c = 0; c < Curves.size(); c++) {
		vtkSmartPointer<vtkPolyData> Line = vtkSmartPointer<vtkPolyData>::New();
		Curves[c].GetSplinePolyData(Line);
		Append -> AddInputData(Line);
	}
	Append -> Update();

	sprintf(_FullPath,"%s_lines.vtk",FilePrefix);
	vtkSmartPointer<vtkPolyDataWriter> Writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	Writer -> SetInputData(Append->GetOutput());
	Writer -> SetFileName(_FullPath);
	Writer -> Write();

	sprintf(_FullPath,"%s.curv",FilePrefix);
	FILE *fc = fopen(_FullPath,"w");
	fprintf(fc,"Curvature\n");
	for (c = 0; c < Curves.size(); c++) {
		Curves[c].ExportCurvatures(fc);
	}
	fclose(fc);

	sprintf(_FullPath,"%s.r2",FilePrefix);
	FILE *fr = fopen(_FullPath,"w");
	fprintf(fr,"Length (um) End-To-Edge Distance (um) End-To-Edge Distance^2 (um^2)\n");
	for (c = 0; c < Curves.size(); c++) {
		fprintf(fr,"%1.6f\t%1.6f\t%1.6f\n",Curves[c].GetLength(),sqrt(Curves[c].GetR2()),Curves[c].GetR2());
	}
	fclose(fr);

	/* POINTS IN SPACE */

	sprintf(_FullPath,"%s.pts",FilePrefix);
	FILE *fg = fopen(_FullPath,"w");
	fprintf(fg,"Header here\n");
	double avg_d1 = 0.0, avg_d2 = 0.0, avg_d12 = 0.0, avg_area = 0.0;
	double max_d1 = 0.0, max_d2 = 0.0, max_d12 = 0.0, max_area = 0.0;
	for (p = 0; p < Nodes.D1.size(); p++) {
		avg_d1 += Nodes.D1[p];
		avg_d2 += Nodes.D2[p];
		avg_d12 += Nodes.D12[p];
		avg_area += Nodes.Area[p];
		max_d1 = (Nodes.D1[p]>max_d1) ? Nodes.D1[p] : max_d1;
		max_d2 = (Nodes.D2[p]>max_d2) ? Nodes.D2[p] : max_d2;
		max_d12 = (Nodes.D12[p]>max_d12) ? Nodes.D12[p] : max_d12;
		max_area = (Nodes.Area[p]>max_area) ? Nodes.Area[p] : max_area;
	}
	fprintf(fg,"%1.6f\t%1.6f\t%1.6f\t%1.6f\t%1.6f\t%1.6f\t%1.6f\t%1.6f\n",avg_d1/Nodes.D1.size(),avg_d2/Nodes.D1.size(),avg_d12/Nodes.D1.size(),avg_area/Nodes.D1.size(),max_d1,max_d2,max_d12,max_area);
	fclose(fg);

	double x, y, z, r[3];
	sprintf(_FullPath,"%s.normcoo",FilePrefix);
	FILE *fn = fopen(_FullPath,"w");
	for (p = 0; p < Nodes.D1.size(); p++) {
		x = Nodes.GetX(p);
		y = Nodes.GetY(p);
		z = Nodes.GetZ(p);
		x = (x-Nodes.GetBoundXmin()) / (Nodes.GetBoundXmax()-Nodes.GetBoundXmin());
		y = (y-Nodes.GetBoundYmin()) / (Nodes.GetBoundYmax()-Nodes.GetBoundYmin());
		z = (z-Nodes.GetBoundZmin()) / (Nodes.GetBoundZmax()-Nodes.GetBoundZmin());
		fprintf(fn,"%1.4f\t%1.4f\t%1.4f\n",x,y,z);
	}
	fclose(fn);

	#ifdef DEBUG
		printf("\tDone!\n");
	#endif

}

/* ====================================================================
	_POINTS METHODS
======================================================================*/


_points::_points(vtkPolyData *Skeleton) {
	Points = vtkSmartPointer<vtkPoints>::New();
	PolyPoints = vtkSmartPointer<vtkPolyData>::New();
	Tree = vtkSmartPointer<vtkKdTreePointLocator>::New();
	double *B = Skeleton -> GetPoints() -> GetBounds();

	int i, j, line, N = Skeleton->GetPoints()->GetNumberOfPoints();
	bool *Used = new bool[N];
	vtkCell *Line;
	double ri[3], rj[3];
	for (i = N; i--;) Used[i] = 0;
	for (i = 6; i--;) Bounds[i] = B[i];

	#ifdef DEBUG
		printf("\tSkeleton Bounds: %1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n",Bounds[0],Bounds[1],Bounds[2],Bounds[3],Bounds[4],Bounds[5]);
	#endif

	for (line = 0; line < Skeleton->GetNumberOfCells(); line++) {
		Line = Skeleton -> GetCell(line);
		i = Line -> GetPointId(0);
		j = Line -> GetPointId(Line->GetNumberOfPoints()-1);
		if (!Used[i]) {
			Skeleton -> GetPoints() -> GetPoint(i,ri);
			Points -> InsertNextPoint(ri);
			Used[i] = 1;
		}
		if (!Used[j]) {
			Skeleton -> GetPoints() -> GetPoint(j,rj);
			Points -> InsertNextPoint(rj);
			Used[j] = 1;
		}
	}
	delete[] Used;

	PolyPoints -> SetPoints(Points);
	Tree -> SetDataSet(PolyPoints);
	Tree -> BuildLocator();
}

void _points::CalculateRelativeDistances() {

	int i, j, k;
	double dij, dik, djk, a, ri[3], rj[3], rk[3], rij[3], rik[3];
	vtkSmartPointer<vtkIdList> List = vtkSmartPointer<vtkIdList>::New();

	for (int i = 0; i < Points->GetNumberOfPoints(); i++) {
		Points -> GetPoint(i,ri);
		Tree -> FindClosestNPoints(3,ri,List);
		j = List -> GetId(1);
		Points -> GetPoint(j,rj);
		k = List -> GetId(2);
		Points -> GetPoint(k,rk);
	
		dij = sqrt(pow(rj[0]-ri[0],2)+pow(rj[1]-ri[1],2)+pow(rj[2]-ri[2],2));
		D1.push_back(dij);
		dik = sqrt(pow(rk[0]-ri[0],2)+pow(rk[1]-ri[1],2)+pow(rk[2]-ri[2],2));
		D2.push_back(dik);
		djk = sqrt(pow(rk[0]-rj[0],2)+pow(rk[1]-rj[1],2)+pow(rk[2]-rj[2],2));
		D12.push_back(djk);

		rij[0] = (rj[0]-ri[0]) / dij;
		rij[1] = (rj[1]-ri[1]) / dij;
		rij[2] = (rj[2]-ri[2]) / dij;

		rik[0] = (rk[0]-ri[0]) / dik;
		rik[1] = (rk[1]-ri[1]) / dik;
		rik[2] = (rk[2]-ri[2]) / dik;

		a = 0.5 * dij * dik * sin( acos( DotProduct(rij,rik) ) );

		Area.push_back(a);

	}
}

double _points::GetX(int i) {
	double r[3];
	Points -> GetPoint(i,r);
	return r[0];
}

double _points::GetY(int i) {
	double r[3];
	Points -> GetPoint(i,r);
	return r[1];
}

double _points::GetZ(int i) {
	double r[3];
	Points -> GetPoint(i,r);
	return r[2];
}


/* ====================================================================
	_CURVE METHODS
======================================================================*/

_curve::_curve(vtkPolyData *Skeleton, vtkIdType cell_id, double spacing) {

	#ifdef DEBUG
		printf("\tGenerating new spline representation...\n");
	#endif

	vtkIdType i;
	double r[3];
	_spacing = spacing;

	vtkPoints *Points = Skeleton -> GetPoints();
	vtkCell *Line = Skeleton -> GetCell(cell_id);
	vtkSmartPointer<vtkPoints> LinePoints = vtkSmartPointer<vtkPoints>::New();
	for (i = 0; i < Line->GetNumberOfPoints(); i++) {
		Points -> GetPoint(Line->GetPointId(i),r);
		LinePoints -> InsertNextPoint(r);
	}
	vtkSmartPointer<vtkCellArray> LineCellArray = vtkSmartPointer<vtkCellArray>::New();
	LineCellArray -> InsertNextCell(Line->GetNumberOfPoints());
	for (i = 0; i < Line->GetNumberOfPoints(); i++) {
		LineCellArray -> InsertCellPoint(i);
	}
	vtkSmartPointer<vtkPolyData> PolyDataLine = vtkSmartPointer<vtkPolyData>::New();
	PolyDataLine -> SetPoints(LinePoints);
	PolyDataLine -> SetLines(LineCellArray);

	vtkSmartPointer<vtkKochanekSpline> SplineModel = vtkSmartPointer<vtkKochanekSpline>::New();

	vtkSmartPointer<vtkSplineFilter> SplineFilter = vtkSmartPointer<vtkSplineFilter>::New();
	SplineFilter->SetSpline(SplineModel);
	SplineFilter->SetInputData(PolyDataLine);
	SplineFilter->SetSubdivideToLength();
	SplineFilter->SetLength(_spacing);
	SplineFilter->Update();

	Spline = SplineFilter -> GetOutput();
	n = Spline -> GetNumberOfPoints();

	#ifdef DEBUG
		printf("\tDone!\n");
	#endif
}

void _curve::SetTangent(vtkIdType id, double t[3]) {
	double tnorm = sqrt(pow(t[0],2)+pow(t[1],2)+pow(t[2],2));
	Tx -> SetTuple1(id,t[0]/tnorm);
	Ty -> SetTuple1(id,t[1]/tnorm);
	Tz -> SetTuple1(id,t[2]/tnorm);
}

void _curve::SetTangent(vtkIdType id, double tx, double ty, double tz) {
	double tnorm = sqrt(pow(tx,2)+pow(ty,2)+pow(tz,2));
	Tx -> SetTuple1(id,tx/tnorm);
	Ty -> SetTuple1(id,ty/tnorm);
	Tz -> SetTuple1(id,tz/tnorm);
	//printf("\t\t%f\t%f\n",tx/tnorm,ty/tnorm);
}

void _curve::GetTangent(vtkIdType id, double t[3]) {
	t[0] = Tx -> GetTuple1(id);
	t[1] = Ty -> GetTuple1(id);
	t[2] = Tz -> GetTuple1(id);
}

void _curve::CalculateTangentVersors(int order) {
	
	dsl = order;
	dsr = (order) ? order : 1;
	
	Tx = vtkSmartPointer<vtkFloatArray>::New();
	Ty = vtkSmartPointer<vtkFloatArray>::New();
	Tz = vtkSmartPointer<vtkFloatArray>::New();
	Tx -> SetNumberOfComponents(1);
	Ty -> SetNumberOfComponents(1);
	Tz -> SetNumberOfComponents(1);
	Tx -> SetNumberOfTuples(n);
	Ty -> SetNumberOfTuples(n);
	Tz -> SetNumberOfTuples(n);

	vtkIdType id;
	double r[3], u[3];
	for (id = dsl; id < n - dsr; id++) {
		Spline -> GetPoints() -> GetPoint(id+dsr,r);
		Spline -> GetPoints() -> GetPoint(id-dsl,u);
		SetTangent(id,r[0]-u[0],r[1]-u[1],r[2]-u[2]);
	}

}

void _curve::CalculateAverageTangentCorrelation() {
	int s;
	vtkIdType id;
	double to[3], tf[3], tcorr;
	int ndist = (int)(0.9*GetLength()/_spacing)-dsr-dsl;
	#ifdef DEBUG
		printf("\tTangent correlation...\n");
		printf("\t>NDist = %d\n",ndist);
	#endif

	for (s = 1; s < ndist; s++) {
		tcorr = 0.0;
		for (id = dsl; id < n - dsr - s; id++) {
			GetTangent(id+0,to);
			GetTangent(id+s,tf);
			tcorr += DotProduct(to,tf);
		}
		TCorr.push_back(tcorr / ((n-dsr-s) - dsl));
	}
}

void _curve::CalculateCurvature() {
	#ifdef DEBUG
		printf("\tCurvature...\n");
	#endif

	vtkIdType id;
	double to[3], tf[3], k;
	Kr = vtkSmartPointer<vtkFloatArray>::New();
	Kr -> SetNumberOfComponents(1);
	Kr -> SetNumberOfTuples(n);
	Kr -> FillComponent(0,0);
	for (id = dsl; id < n - dsr - 1; id++) {
		GetTangent(id+0,to);
		GetTangent(id+1,tf);
		//k = ||dT/ds||
		k = sqrt(pow(tf[0]-to[0],2)+pow(tf[1]-to[1],2)+pow(tf[2]-to[2],2)) / _spacing;
		Kr -> SetTuple1(id,k);
	}
}
void _curve::GetSplinePolyData(vtkPolyData *Line) {
	double r[3];
	vtkSmartPointer<vtkPoints> LinePoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> LineCellArray = vtkSmartPointer<vtkCellArray>::New();
	LineCellArray -> InsertNextCell(n-dsr-dsl);
	vtkSmartPointer<vtkFloatArray> Curvature = vtkSmartPointer<vtkFloatArray>::New();
	Curvature -> SetNumberOfComponents(1);
	Curvature -> SetNumberOfTuples(n-dsr-dsl);
	Curvature -> FillComponent(0,0);
	Curvature -> SetName("Curvature");

	for (vtkIdType id = dsl; id < n - dsr; id++) {
		Spline -> GetPoint(id,r);
		LinePoints -> InsertNextPoint(r);
		LineCellArray -> InsertCellPoint(id-dsl);
		Curvature -> SetTuple1(id-dsl,Kr->GetTuple1(id));
	}

	Line -> SetPoints(LinePoints);
	Line -> SetLines(LineCellArray);
	Line -> GetPointData() -> SetScalars(Curvature);
}

void _curve::ExportCurvatures(FILE *f) {
	for (vtkIdType id = dsl; id < n - dsr - 1; id++) {
		fprintf(f,"%1.6f\n",Kr->GetTuple1(id));
	}
}

void _curve::CalculateEndFluctuations() {
	double ro[3], rf[3];
	Spline -> GetPoints() -> GetPoint(0,ro);
	Spline -> GetPoints() -> GetPoint(Spline->GetNumberOfPoints()-1,rf);
	R2 = pow(rf[0]-ro[0],2) + pow(rf[1]-ro[1],2) + pow(rf[2]-ro[2],2);
}

/* =================================================================
   MAIN
   =================================================================*/


int main(int argc, char *argv[]) {     

	srand(getpid());
	char _RootFolder[256] = {""};

	for (int i = 0; i < argc; i++) {
		if (!strcmp(argv[i],"-path")) {
			sprintf(_RootFolder,"%s//",argv[i+1]);
		}
	}

	// Generating list of files to run
	char _cmd[256];
	sprintf(_cmd,"ls %s*_skeleton.vtk | sed -e 's/_skeleton.vtk//' > %smitographgeo.files",_RootFolder,_RootFolder);
	system(_cmd);

	char _PrefixFile[256];
	char _PrefixList[256];
	
	// List of files to run

	sprintf(_PrefixList,"%smitographgeo.files",_RootFolder);	
	FILE *f = fopen(_PrefixList,"r");

	// Main loop

	while (fgets(_PrefixFile,256, f) != NULL) {
		_PrefixFile[strcspn(_PrefixFile, "\n" )] = '\0';
		
		printf("%s\n",_PrefixFile);
		GeometricalAnalysis(_PrefixFile);
	
	}
	fclose(f);

	return 0;
}