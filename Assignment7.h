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
class Assignment7 : public Experiment
{
    GEOX_CLASS(Assignment7)

//Constructor / Destructor
public:
    Assignment7();
    virtual ~Assignment7();

//Methods
public:
	void ReadFieldFromFile();

	//Domain Decomposition and Change-Of-Sign test:
	bool IsZeroPossible(int x, int y);
	//Vector2f FindZero(int x, int y);
	bool IsZeroPossible(vector<Vector2f> points);
	void FindZero(vector<Vector2f> points);

	void FindCriticalPoints();
	void ClassifyCriticalPoints();
	void ComputeSeparatrices();


	//Legacy methods
	void GenerateTexture();
	void GenerateKernel();
	void ClassicLIC();
	void FastLIC();
	vector<Vector2f> GenerateStreamLines(Vector2f startPoint, int pointNumber);
	float convolveKernel(int startIndex, vector<Vector2f> line);
	void EnhanceContrast();
	
	//Basic methods
	Vector2f RungeKuttaIntegration(Vector2f pos);
	Vector2f RungeKuttaIntegration(Vector2f pos, bool forwards);


    virtual QWidget* createViewer();

protected:
	

//Attributes
public:
	//Filename to read from
	string fileName;

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
	//Storing the fields
	VectorField2 vField;
	ScalarField2 sField;

	//Zero-finding threshold

	//Storing the critical points
	vector<Vector2f> AllCriticals;
	
	vector<Vector2f> Source;
	vector<Vector2f> FocusRepelling;
	vector<Vector2f> Saddle;
	vector<Vector2f> Center;
	vector<Vector2f> Sink;
	vector<Vector2f> FocusAttracting;

	//Threshold value for finding zeroes
	float ZeroThreshold;

	//LIC texture parameters:
	//Greyscale or B&W
	bool greyScale;
	//Textures: input and output respectively
	ScalarField2 texture;
	ScalarField2 LICtexture;
	ScalarField2 TEMPtexture;
	//vField to keep track of the visited pixels in fastLIC
	ScalarField2 visited;

	//Values for the generated kernel
	int kernelLength;
	vector<float> kernelValues;

	//Texture resolutions
	int xPowerOfTwo;
	int yPowerOfTwo;
	//Seed for random
	int textureSeed;

	bool ContrastEnhancement;
	float DesiredMean;
	float DesiredDeviation;

protected:
    GLGeometryViewer* viewer;
};


