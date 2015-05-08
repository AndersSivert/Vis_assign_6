//---------------------------------------------------------------------------
#pragma once
//---------------------------------------------------------------------------
#include "Experiment.h"
#include "LinearAlgebra.h"
#include "GLGeometryViewer.h"
//---------------------------------------------------------------------------

#include "Field2.hpp"

/// KTH DD2257 Visualization
///	Assignment 6
///
/// Authors: Jim Eriksson, Martin Rockstrom, Anders Sivertsson
///
class Assignment6 : public Experiment
{
    GEOX_CLASS(Assignment6)

//Constructor / Destructor
public:
    Assignment6();
    virtual ~Assignment6();

//Methods
public:
	void Comparison();
	void ReadFieldFromFile();
	void DrawStreamLines();

	void GenerateTexture();

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

	bool greyScale;
	ScalarField2 texture;


protected:
    GLGeometryViewer* viewer;
};


