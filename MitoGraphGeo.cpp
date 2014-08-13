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

double GetLength(vtkCell *Line, vtkPoints *Points, bool plot) {
	double ds = 0.0, d2 = 0.0;
	double r[3], u[3], d, length = 0.0;
	for (int i = 1; i < Line -> GetNumberOfPoints(); i++) {
		Points -> GetPoint(Line->GetPointId(i-1),u);
		Points -> GetPoint(Line->GetPointId(i-0),r);
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
/*
struct _point{
	double x, y, z, s;
	
	double tx, ty, tz;
	double k, tcorr;

	_point (double _x, double _y, double _z, double _s) {
		x = _x; y = _y; z = _z; s = _s;
	}
	_point (double r[3], double _s) {
		x = r[0]; y = r[1]; z = r[2]; s = _s;
	}
};

class _curve {
private:
	int dsl, dsr;
public:
	double length;
	std::vector<_point> Points;
	void SetOrder(int order) {
		dsl = order;
		dsr = (order) ? order : 1;
	}
	_curve(vtkCell *Line, vtkPoints *VtkPts);
	void CalculateTangentVectors();
	void CalculateTangentCorrelation();
	void PrintProperties();
	void CalculateCurvature();
	void ExportPolyData(const char FileName[]);
};

_curve::_curve(vtkCell *Line, vtkPoints *VtkPts) {
	vtkSmartPointer<vtkPoints> SplinePoints = vtkSmartPointer<vtkPoints>::New();

	int i;
	length = 0.0;		
	dsl = 0; dsr = 1;
	double r[3], u[3];

	for (i = 0; i < Line->GetNumberOfPoints(); i++) {
		VtkPts -> GetPoint(Line->GetPointId(i-0),r);
		SplinePoints -> InsertNextPoint(r);
	}

	vtkSmartPointer<vtkCellArray> Line = vtkSmartPointer<vtkCellArray>::New();
	Line -> InsertNextCell(c);
	for (i = 0; i < c; i++) {
		Line -> InsertCellPoint(i);
		Scalar -> SetTuple1(i,Points[i+dsl].k);
	}
	vtkSmartPointer<vtkPolyData> PolyLine = vtkSmartPointer<vtkPolyData>::New();
	PolyLine -> SetPoints(vtkPts);
	PolyLine -> SetLines(Line);
	PolyLine -> GetPointData() -> SetScalars(Scalar);


	vtkSmartPointer<vtkParametricSpline> Spline = vtkSmartPointer<vtkParametricSpline>::New();
	Spline -> SetPoints(SplinePoints);
	Spline -> ClosedOff();
	Spline -> ParameterizeByLengthOn();

	vtkSmartPointer<vtkSplineFilter> SplineFilter = vtkSmartPointer<vtkSplineFilter>::New();
	SplineFilter->SetSpline(Spline);
	SplineFilter->SetInputConnection(SplinePoints);
	SplineFilter->SetSubdivideToLength();
	SplineFilter->SetLength( 0.18648 );
	SplineFilter->Update( );

	u[0] = 0.0;
	do {
		Spline -> Evaluate(u,r,NULL);
		_point point(r,0.0);
		Points.push_back(point);
		u[0] += 0.01;
	} while (fabs(1.0-u[0])>1E-5);

	double dij = 0.0;
	double dijs = 0.0;
	double dij2 = 0.0;
	for (i = 1; i < Points.size(); i++) {
		dij = sqrt(pow(Points[i].x-Points[i-1].x,2)+pow(Points[i].y-Points[i-1].y,2)+pow(Points[i].z-Points[i-1].z,2));
		dijs += dij;
		dij2 += dij*dij;
	}

	dijs /= Points.size() - 1;
	dij2 /= Points.size() - 1;

	printf("%f\t%f\n",dijs,sqrt(dij2-dij*dij));

	//VtkPts -> GetPoint(Line->GetPointId(0),r);
	//_point point(r,s);
	//Points.push_back(point);

	//double v[3];

}

void _curve::CalculateTangentVectors() {
	int i;
	double Tn;
	for (i = dsl; i < Points.size()-dsr; i++) {
		Points[i].tx = Points[i+dsr].x - Points[i-dsl].x;
		Points[i].ty = Points[i+dsr].y - Points[i-dsl].y;
		Points[i].tz = Points[i+dsr].z - Points[i-dsl].z;
		Tn = sqrt(pow(Points[i].tx,2)+pow(Points[i].ty,2)+pow(Points[i].tz,2));
		Points[i].tx /= Tn;
		Points[i].ty /= Tn;
		Points[i].tz /= Tn;
		printf("Tn = %f\n",Tn);
	}
}

void _curve::CalculateTangentCorrelation() {
	int i;
	for (i = dsl; i < Points.size()-dsr; i++) {
		Points[i].tcorr = Points[dsl].tx * Points[i].tx + Points[dsl].ty * Points[i].ty + Points[dsl].tz * Points[i].tz;
		Points[i].tcorr = (Points[i].tcorr > 1.00) ? 1.00 : Points[i].tcorr;
	}
}

void _curve::CalculateCurvature() {
	int i;
	for (i = dsl; i < Points.size()-dsr; i++) {
		Points[i].k = fabs( acos(Points[i].tcorr) / (Points[i+dsr].s - Points[i-dsl].s) );
	}
}

void _curve::PrintProperties() {
	int i;
	for (i = dsl; i < Points.size()-dsr; i++) {
		printf("%f\t%f\t%f\n",Points[i].s,Points[i].tcorr,Points[i].k);
	}
}

void _curve::ExportPolyData(const char FileName[]) {
	int i, c = 0;;
	vtkSmartPointer<vtkPoints> vtkPts = vtkSmartPointer<vtkPoints>::New();
	for (i = dsl; i < Points.size()-dsr; i++) {
		vtkPts -> InsertNextPoint(Points[i].x,Points[i].y,Points[i].z);
		c++;
	}
	vtkSmartPointer<vtkFloatArray> Scalar = vtkSmartPointer<vtkFloatArray>::New();
	Scalar -> SetNumberOfComponents(1);
	Scalar -> SetNumberOfTuples(c);
	vtkSmartPointer<vtkCellArray> Line = vtkSmartPointer<vtkCellArray>::New();
	Line -> InsertNextCell(c);
	for (i = 0; i < c; i++) {
		Line -> InsertCellPoint(i);
		Scalar -> SetTuple1(i,Points[i+dsl].k);
	}
	vtkSmartPointer<vtkPolyData> PolyLine = vtkSmartPointer<vtkPolyData>::New();
	PolyLine -> SetPoints(vtkPts);
	PolyLine -> SetLines(Line);
	PolyLine -> GetPointData() -> SetScalars(Scalar);
	vtkSmartPointer<vtkPolyDataWriter> Writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	Writer -> SetFileName(FileName);
	Writer -> SetInputData(PolyLine);
	Writer -> Write();

}
*/
void NewFuckingRoutine(vtkCell *Line, vtkPoints *Points) {
	
	int i;
	double r[3];

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
/*
	vtkSmartPointer<vtkParametricSpline> Spline = vtkSmartPointer<vtkParametricSpline>::New();
	Spline -> SetPoints(LinePoints);
	Spline -> ClosedOff();
	Spline -> ParameterizeByLengthOn();
*/
	vtkSmartPointer<vtkKochanekSpline> SplineModel = vtkSmartPointer<vtkKochanekSpline>::New();

	vtkSmartPointer<vtkSplineFilter> SplineFilter = vtkSmartPointer<vtkSplineFilter>::New();
	SplineFilter->SetSpline(SplineModel);
	SplineFilter->SetInputData(PolyDataLine);
	SplineFilter->SetSubdivideToLength();
	SplineFilter->SetLength(0.01);
	SplineFilter->Update();

	vtkPolyData *Spline = SplineFilter -> GetOutput();

	printf("#Points 0: %d\n",(int)Line->GetNumberOfPoints());
	printf("#Points 1: %d\n",(int)Spline->GetNumberOfPoints());

	printf("#Length 0: %f\n",GetLength(Line,Points,false));
	printf("#Length 1: %f\n",GetLength(Spline->GetCell(0),Spline->GetPoints(),true));


}

void TestRun(const char FilePrefix[]) {
	
	char _FullPath[256];
	sprintf(_FullPath,"%s_skeleton.vtk",FilePrefix);

	vtkSmartPointer<vtkPolyDataReader> PolyReader = vtkSmartPointer<vtkPolyDataReader>::New();
	PolyReader -> SetFileName(_FullPath);
	PolyReader -> Update();

	vtkPolyData *Skell = PolyReader -> GetOutput();
	vtkPoints *VtkPts = Skell -> GetPoints();

	int NLines = Skell -> GetNumberOfCells();

	//std::vector<_curve> Curves;

	int line;
	for (line = 0; line < NLines; line++) {
		if ( GetLength(Skell->GetCell(line),VtkPts,false) > 3.0 ) {
			NewFuckingRoutine(Skell->GetCell(line),VtkPts);
		}
/*		_curve curve(Skell->GetCell(line),VtkPts);
		curve.SetOrder(2);
		if (curve.length>3.0) {
			curve.CalculateTangentVectors();
			curve.CalculateTangentCorrelation();
			curve.CalculateCurvature();
			Curves.push_back(curve);
		}
*/	}

/*	int curve, co = -1;
	double lmax = 0.0;
	for (curve = 0;  curve < Curves.size(); curve++) {
		//Curves[curve].PrintProperties();
		printf("%f\n",Curves[curve].length);
		if (Curves[curve].length > lmax) {
			lmax = Curves[curve].length;
			co = curve;
		}
	}
	if (co>=0) Curves[co].ExportPolyData("temp.vtk");
*/
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
		
		TestRun(_PrefixFile);
		printf("%s\n",_PrefixFile);
	
	}
	fclose(f);

	return 0;
}