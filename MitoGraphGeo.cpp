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
	void ExportResults(std::vector<_curve> Curves, const char FilePrefix[]);
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
	
	if (Curves.size()) {
		ExportResults(Curves,FilePrefix);
	} else {
		printf("No edges long enough found.\n");
	}

	Curves.clear();
}

void ExportResults(std::vector<_curve> Curves, const char FilePrefix[]) {

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

	#ifdef DEBUG
		printf("\tDone!\n");
	#endif

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