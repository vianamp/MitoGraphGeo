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
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkBitArray.h>
#include <vtkFloatArray.h>
#include <vtkCellLocator.h>
#include <vtkTriangle.h>
#include <vtkKdTreePointLocator.h>
#include <vtkKdTree.h>
#include <vtkSplineFilter.h>
#include <vtkParametricSpline.h>
#include <vtkKochanekSpline.h>
#include <vtkParametricFunctionSource.h>


#define _min(a,b) ((a<b)?a:b)
#define _max(a,b) ((a>b)?a:b)
#define _eps 1E-6

class _point{
	private:
		double x, y, z;

	public:
		_point (double _x, double _y, double _z) {
			x = _x; y = _y; z = _z;
		}
		_point (double r[3]) {
			x = r[0]; y = r[1]; z = r[2];
		}
		void Get(double r[3]) {
			r[0] = x; r[1] = y; r[2] = z;
		}
		void Set(double r[3]) {
			x = r[0]; y = r[1]; z = r[2];
		}
};

class _points{
	public:
		_points() {}
		~_points() {}
		int Size() {
			return (int)Points.size();
		}
		void GetPoint(int i, double r[3]) {
			Points[i].Get(r);
		}
		void SetPoint(int i, double r[3]) {
			Points[i].Set(r);
		}		
		void AddPoint(double x, double y, double z) {
			_point point(x,y,z);
			Points.push_back(point);
		}
		void ProjectOnThis(vtkPolyData *Surface);
		void ShuffleOnThis(vtkPolyData *Surface);
		void Get1stNeighDist(double *NeighDist);
		void ExportAsPolyData(const char FileName[]);

	private:
		std::vector<_point> Points;
};

void _points::ProjectOnThis(vtkPolyData *Surface) {
	vtkSmartPointer<vtkKdTree> KdTree = vtkSmartPointer<vtkKdTree>::New();
	vtkPoints *EllPoints = Surface -> GetPoints();
	KdTree -> BuildLocatorFromPoints( EllPoints );

	int i, p;
	double r[3], d;
	for (i = 0; i < Points.size(); i++) {
		GetPoint(i,r);
		p = KdTree -> FindClosestPoint(r,d);
		EllPoints -> GetPoint(p,r);
		SetPoint(i,r);
		GetPoint(i,r);
	}
}

void _points::ShuffleOnThis(vtkPolyData *Surface) {
	int i, p;
	double r[3];
	for (i = 0; i < Points.size(); i++) {
		GetPoint(i,r);
		p = rand() % Surface -> GetNumberOfPoints();
		Surface -> GetPoint(p,r);
		SetPoint(i,r);
	}
}


void _points::Get1stNeighDist(double *NeighDist) {
	int i, j;
	double ri[3], rj[3], d;
	vtkSmartPointer<vtkPoints> TmpPoints = vtkSmartPointer<vtkPoints>::New();
	for (i = 0; i < Points.size(); i++) {
		GetPoint(i,ri);
		TmpPoints -> InsertPoint(i,ri);
	}
	vtkSmartPointer<vtkPolyData> TmpPolydt = vtkSmartPointer<vtkPolyData>::New();
	TmpPolydt -> SetPoints(TmpPoints);

    vtkSmartPointer<vtkKdTreePointLocator> Tree = vtkSmartPointer<vtkKdTreePointLocator>::New();
    Tree -> SetDataSet(TmpPolydt);
    Tree -> BuildLocator();

    vtkSmartPointer<vtkIdList> List = vtkSmartPointer<vtkIdList>::New();
    
    for (i = 0; i < Points.size(); i++) {
        GetPoint(i,ri);
        Tree -> FindClosestNPoints(2,ri,List);
        j = List -> GetId(1);
        TmpPoints -> GetPoint(j,rj);
		d = sqrt( pow(ri[0]-rj[0],2) + pow(ri[1]-rj[1],2) + pow(ri[2]-rj[2],2) );
		NeighDist[i] = d;
	}

}

void _points::ExportAsPolyData(const char FileName[]) {
	int i;
	double ri[3];
	vtkSmartPointer<vtkPoints> TmpPoints = vtkSmartPointer<vtkPoints>::New();
	for (i = 0; i < Points.size(); i++) {
		GetPoint(i,ri);
		TmpPoints -> InsertPoint(i,ri);
	}
	vtkSmartPointer<vtkPolyData> TmpPolydt = vtkSmartPointer<vtkPolyData>::New();
	TmpPolydt -> SetPoints(TmpPoints);

	vtkSmartPointer<vtkPolyDataWriter> Writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	Writer -> SetInputData(TmpPolydt);
	Writer -> SetFileName(FileName);
	Writer -> Write();
}

void TestRun(const char FilePrefix[]) {
	
	_points Nodes;

	float x, y, z;
	char _FullPath[256];
	sprintf(_FullPath,"%s.coo",FilePrefix);
	FILE *fc = fopen(_FullPath,"r");
	while ( fscanf(fc,"%f %f %f",&x,&y,&z) != EOF ) {
		Nodes.AddPoint(0.056*x,0.056*y,0.2*z);
	}
	fclose(fc);

	sprintf(_FullPath,"%s_ellipsoid.vtk",FilePrefix);

	vtkSmartPointer<vtkPolyDataReader> EllReader = vtkSmartPointer<vtkPolyDataReader>::New();
	EllReader -> SetFileName(_FullPath);
	EllReader -> Update();
	vtkPolyData *Ellipsoid = EllReader -> GetOutput();

	printf("Number of Nodes = %d\n",Nodes.Size());
	printf("Ellipsoid points = %d\n",(int)Ellipsoid->GetNumberOfPoints());

	int i, run = 0;
	double *NeighDist = new double[Nodes.Size()];
	
	sprintf(_FullPath,"%s.1stneighdist",FilePrefix);
	FILE *fd = fopen(_FullPath,"w");
	
	Nodes.ProjectOnThis(Ellipsoid);
	do {
		Nodes.Get1stNeighDist(NeighDist);
		Nodes.ShuffleOnThis(Ellipsoid);

		for (i = 0; i < Nodes.Size(); i++) {
			fprintf(fd,"%d\t%1.5f\n",(run==0)?1:0,NeighDist[i]);
		}

		run++;
	} while (run < 500);
	fclose(fd);

	delete[] NeighDist;

}

void GetSurfaceEntropy(const char FilePrefix[]) {

	char _FullPath[256];
	sprintf(_FullPath,"%s_ellipsoid.vtk",FilePrefix);

	vtkSmartPointer<vtkPolyDataReader> EllReader = vtkSmartPointer<vtkPolyDataReader>::New();
	EllReader -> SetFileName(_FullPath);
	EllReader -> Update();
	vtkPolyData *Ellipsoid = EllReader -> GetOutput();

	printf("Ellipsoid points = %d\n",(int)Ellipsoid->GetNumberOfPoints());

	sprintf(_FullPath,"%s_surface.vtk",FilePrefix);

	vtkSmartPointer<vtkPolyDataReader> SurReader = vtkSmartPointer<vtkPolyDataReader>::New();
	SurReader -> SetFileName(_FullPath);
	SurReader -> Update();
	vtkPolyData *Surface = SurReader -> GetOutput();

	printf("Surface points = %d\n",(int)Surface->GetNumberOfPoints());

    vtkSmartPointer<vtkCellLocator> CellLocator = vtkSmartPointer<vtkCellLocator>::New();
    CellLocator -> SetDataSet(Ellipsoid);
    CellLocator -> BuildLocator();

    vtkSmartPointer<vtkFloatArray> Density = vtkSmartPointer<vtkFloatArray>::New();
    Density -> SetNumberOfComponents(1);
    Density -> SetNumberOfTuples(Ellipsoid->GetNumberOfCells());
    Density -> FillComponent(0,0);

    double r[3], u[3], d, a;
    vtkIdType ido;
    int subido;
    for (vtkIdType id = 0; id < Surface->GetNumberOfPoints(); id++) {
    	Surface -> GetPoint(id,r);
        CellLocator -> FindClosestPoint(r,u,ido,subido,d);
        Density -> SetTuple1(ido,Density->GetTuple1(ido)+1.0);
    }
    Density -> Modified();

    double rho, mass = 0.0;
    for (vtkIdType id = 0; id < Ellipsoid->GetNumberOfCells(); id++) {
        vtkTriangle *T = dynamic_cast<vtkTriangle*>(Ellipsoid->GetCell(id));
		rho = Density -> GetTuple1(id) / (T->ComputeArea()*Surface->GetNumberOfPoints());
		mass += rho;
        Density -> SetTuple1(id,rho);
    }
    Density -> Modified();

    printf("Total mass: %1.5f\n",mass);
	sprintf(_FullPath,"%s.entropy",FilePrefix);
	FILE *f = fopen(_FullPath,"w");
    for (vtkIdType id = 0; id < Ellipsoid->GetNumberOfCells(); id++) {
        Density -> SetTuple1(id,Density->GetTuple1(id)/mass);
        fprintf(f,"%1.8f\n",Density->GetTuple1(id)/mass);
    }
    Density -> Modified();
    fclose(f);

    Ellipsoid -> GetCellData() -> SetScalars(Density);

	sprintf(_FullPath,"%s_ellipsoid2.vtk",FilePrefix);
	vtkSmartPointer<vtkPolyDataWriter> Writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	Writer -> SetInputData(Ellipsoid);
	Writer -> SetFileName(_FullPath);
	Writer -> Write();

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
		
		GetSurfaceEntropy(_PrefixFile);
		printf("%s\n",_PrefixFile);
	
	}
	fclose(f);

	return 0;
}