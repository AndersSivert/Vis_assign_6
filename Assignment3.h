//---------------------------------------------------------------------------
#pragma once
//---------------------------------------------------------------------------
#include "Experiment.h"
#include "LinearAlgebra.h"
#include "GLGeometryViewer.h"
//---------------------------------------------------------------------------

#include "Field2.hpp"

/// This is an example experiment.
///
/// The code is meant to demonstrate how
///  to use the GeoX framework
///
class Assignment3 : public Experiment
{
    GEOX_CLASS(Assignment3)

//Constructor / Destructor
public:
    Assignment3();
    virtual ~Assignment3();

//Methods
public:
    void DrawScalarField();
	void ReadData();
	void DrawData();
	void DrawIsoMid();
	void DrawIsoAsympt();
    virtual QWidget* createViewer();

protected:
	bool static PointSorter(const Point2D & p1,const Point2D & p2);

//Attributes
public:

    ///File name of the scalar field
    string ScalarfieldFilename;

	//Bool to show the grid
	bool DrawGrid;

	//Bool to show colored points for each vertex value
	bool DrawPoints;

	//Bool to draw the isocontours
	bool DrawIsoContours;

	//Bool for midpoint or asymptote decider
	bool AsymptoticDecider;

	bool AutoIsoValues;
	string manualIsoValues;
	int numberOfIsocontours;

	//Vector that contains all the iso-values
	vector<float> isoValues;
	//Vector that countains the colours for the iso-values
	vector<Vector4f> isoColours;

	ScalarField2 field;


protected:
    GLGeometryViewer* viewer;
};


