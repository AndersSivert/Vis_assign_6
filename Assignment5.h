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
class Assignment5 : public Experiment
{
    GEOX_CLASS(Assignment5)

//Constructor / Destructor
public:
    Assignment5();
    virtual ~Assignment5();

//Methods
public:
	void Comparison();
	void ReadFieldFromFile();
	void DrawStreamLines();
	Vector2f static GetFieldValue1(Vector2f pos);
	Vector2f RungeKuttaIntegration(Vector2f pos);
	Vector2f RungeKuttaIntegration(Vector2f pos, bool forwards);
	Vector2f EulerIntegration(Vector2f pos);
	Vector2f EulerIntegration(Vector2f pos, bool forwards);

    virtual QWidget* createViewer();

protected:
	

//Attributes
public:
	
	//Use RK for tasks after 1?
	bool RungeKutta;

	string fileName;
	string startingPoints;
	bool randomPoints;

	//Euler parameters
	int EulerSteps;
	float EulerStepSize;

	//RK parameters
	int RKSteps;
	float RKStepSize;
    
	//Streamline parameters
	float maxLength;
	bool directionField;
	bool grid;
	float minimumMagnitude;
	float maxAdjustment;
	bool showPoints;
	bool magnitudeColor;
	int adjustSteps;
	int adjustArea;

	vector<Vector4f> colors;
	int colorIndex;

	bool readField;

	VectorField2 field;

protected:
    GLGeometryViewer* viewer;
};


