//---------------------------------------------------------------------------
#include "stdafx.h"
//---------------------------------------------------------------------------
#include "Assignment7.h"
//---------------------------------------------------------------------------
#include "Properties.h"
#include "GLGeometryViewer.h"
#include "GeoXOutput.h"
//---------------------------------------------------------------------------

#include <limits>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

IMPLEMENT_GEOX_CLASS( Assignment7, 0)
{
    BEGIN_CLASS_INIT( Assignment7 );
	ADD_SEPARATOR("Vectorfield file name")
	ADD_STRING_PROP(fileName, 0)
	
	ADD_SEPARATOR("Texture parameters")
	ADD_BOOLEAN_PROP(greyScale,0)
	ADD_INT32_PROP(xPowerOfTwo,0)
	ADD_INT32_PROP(yPowerOfTwo,0)
	ADD_INT32_PROP(textureSeed,0)
	ADD_BOOLEAN_PROP(ContrastEnhancement,0)
	ADD_FLOAT32_PROP(DesiredMean,0)
	ADD_FLOAT32_PROP(DesiredDeviation,0)


	ADD_SEPARATOR("Runge-Kutta parameters")
	ADD_INT32_PROP(RKSteps,0)
	//ADD_FLOAT32_PROP(RKStepSize,0)
	ADD_FLOAT32_PROP(minimumMagnitude,0)
	
	ADD_SEPARATOR("Kernel parameters")
	ADD_INT32_PROP(kernelLength,0)

	
	ADD_NOARGS_METHOD(Assignment7::ReadFieldFromFile)
	ADD_NOARGS_METHOD(Assignment7::GenerateTexture)
	ADD_NOARGS_METHOD(Assignment7::ClassicLIC)
	ADD_NOARGS_METHOD(Assignment7::FastLIC)
	ADD_NOARGS_METHOD(Assignment7::EnhanceContrast)

	
	

    
}

QWidget* Assignment7::createViewer()
{
    viewer = new GLGeometryViewer();
    return viewer;
}

Assignment7::Assignment7()
{
    viewer = NULL;

	fileName = "";		//Martin
	//fileName = "";	//Anders
	//fileName = "";	//Jim

	readField = false;
	
	xPowerOfTwo = 0;
	yPowerOfTwo = 0;
	
	textureSeed = 0;

	EulerSteps = 100;
	EulerStepSize = 0.1;
	
	RKSteps = 100;
	RKStepSize = 0.1;

	maxLength = 10;
	directionField = true;
	magnitudeColor = false;
	showPoints = true;
	grid = true;
	minimumMagnitude = 0.05;
	adjustArea = 0;
	adjustSteps = 5;
	maxAdjustment = 0.05;

	kernelLength = 10;

	ContrastEnhancement = false;
	DesiredMean = 0.5;
	DesiredDeviation = 0.1;
}

Assignment7::~Assignment7() {}

/// Returns the next power of two
///	Kudos to ExampleExperimentFields
///	Unchanged --> consider importing instead?
namespace
{
    ///Returns the next power-of-two
    int32 NextPOT(int32 n)
    {
        n--;
        n |= n >> 1;   // Divide by 2^k for consecutive doublings of k up to 32,
        n |= n >> 2;   // and then or the results.
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        n++;           // The result is a number of 1 bits equal to the number
                       // of bits in the original number, plus 1. That's the
                       // next highest power of 2.
        return n;
    }
}

Vector2f Assignment7::RungeKuttaIntegration(Vector2f pos) {
	Vector2f p1 =  vField.sample(pos[0],pos[1]);
	Vector2f p2 =  vField.sample(pos[0]+p1[0]*RKStepSize/2,pos[1]+p1[1]*RKStepSize/2);
	Vector2f p3 =  vField.sample(pos[0]+p2[0]*RKStepSize/2,pos[1]+p2[1]*RKStepSize/2);
	Vector2f p4 =  vField.sample(pos[0]+p3[0]*RKStepSize,pos[1]+p3[1]*RKStepSize);
	return pos + (p1 + p2*2 + p3*2 + p4)*RKStepSize/6;
}
Vector2f Assignment7::RungeKuttaIntegration(Vector2f pos, bool forwards) {
	float dt;
	if(forwards) {
		dt = RKStepSize;
	} else {
		dt = -RKStepSize;
	}

	if(directionField) {
		Vector2f p1 =  vField.sample(pos[0],pos[1]);
		p1.normalize();
		Vector2f p2 =  vField.sample(pos[0]+p1[0]*dt/2,pos[1]+p1[1]*dt/2);
		p2.normalize();
		Vector2f p3 =  vField.sample(pos[0]+p2[0]*dt/2,pos[1]+p2[1]*dt/2);
		p3.normalize();
		Vector2f p4 =  vField.sample(pos[0]+p3[0]*dt,pos[1]+p3[1]*dt);
		p4.normalize();
		return pos + (p1 + p2*2 + p3*2 + p4)*dt/6;
	} else {
		Vector2f p1 =  vField.sample(pos[0],pos[1]);
		Vector2f p2 =  vField.sample(pos[0]+p1[0]*dt/2,pos[1]+p1[1]*dt/2);
		Vector2f p3 =  vField.sample(pos[0]+p2[0]*dt/2,pos[1]+p2[1]*dt/2);
		Vector2f p4 =  vField.sample(pos[0]+p3[0]*dt,pos[1]+p3[1]*dt);
		return pos + (p1 + p2*2 + p3*2 + p4)*dt/6;
	}
}

void Assignment7::ReadFieldFromFile(){
	if (!vField.load(fileName))
    {
		if(!sField.load(fileName)) {
			output << "Error loading field file " << fileName << "\n";
			return;
		} else {
			// If we have a scalar field, use the gradient as a vector field.
			vField.clear();
			vField.init(sField.boundMin(),sField.boundMax(),sField.dims());
			for(int x=0;x<sField.dims()[0];x++) {
				for(int y=0;y<sField.dims()[1];y++) {
					vField.setNode(x,y,sField.sampleGradient(sField.nodePosition(x,y)));
				}
			}
			readField = true;
		}
    } else {
		readField = true;
	}

}

bool Assignment7::IsZeroPossible(int x, int y) {
	//TODO: THIS IS PLACEHOLDER
	return false;
}
Vector2f Assignment7::FindZero(int x, int y) {
	//TODO: THIS IS PLACEHOLDER
	Vector2f result;
	result.setZero();
	return result;
}

void Assignment7::FindCriticalPoints() {
	//TODO:Implement this method
}

void Assignment7::ClassifyCriticalPoints() {
	//TODO:Implement this method
}

void Assignment7::ComputeSeparatrices() {
	//TODO:Implement this method
}

//Legacy methods below this point
//Could be useful for bonus points

void Assignment7::GenerateTexture() {
	viewer->clear();
	texture.clear();
	visited.clear();
	int iWidth = NextPOT(vField.boundMax()[0]-vField.boundMin()[0]);
	if(xPowerOfTwo>0) {
		iWidth = pow(2.0,xPowerOfTwo);
	}
	int iHeight = NextPOT(vField.boundMax()[1]-vField.boundMin()[1]);
	if(yPowerOfTwo>0) {
		iHeight = pow(2.0,yPowerOfTwo);
	}

	srand(textureSeed);
	
	texture.init(vField.boundMin(),vField.boundMax(),makeVector2ui(iWidth,iHeight));
	LICtexture.init(vField.boundMin(),vField.boundMax(),makeVector2ui(iWidth,iHeight));
	visited.init(vField.boundMin(),vField.boundMax(),makeVector2ui(iWidth,iHeight));
	if(greyScale) {
		for(size_t i = 0; i < iWidth;i++) {
			for(size_t j = 0; j < iHeight; j++) {
				visited.setNodeScalar(i,j,0.0);
				texture.setNodeScalar(i,j,(float)rand()/RAND_MAX);
			}
		}
	} else {
		for(size_t i = 0; i < iWidth;i++) {
			for(size_t j = 0; j < iHeight; j++) {
				visited.setNodeScalar(i,j,0.0);
				float value = (float)rand()/RAND_MAX;
				texture.setNodeScalar(i,j, (value > 0.5) ? 1 : 0);
			}
		}
	}
	output << "segment length 1: " << min((texture.boundMax()[0] - texture.boundMin()[0])/iWidth,(texture.boundMax()[1] - texture.boundMin()[1])/iHeight) << "\n";
	output << "segment length 2: " << min((texture.nodePosition(0,0)-texture.nodePosition(1,0)).getSqrNorm(),(texture.nodePosition(0,0)-texture.nodePosition(0,1)).getSqrNorm()) << "\n";
	RKStepSize = min((texture.boundMax()[0] - texture.boundMin()[0])/iWidth,(texture.boundMax()[1] - texture.boundMin()[1])/iHeight);
	viewer->setTextureGray(texture.getData());
	viewer->refresh();
}

void Assignment7::GenerateKernel() {
	kernelValues.clear();
	//TODO: Implement different kernels, kernel sizes etc.
	//Currently box kernel of size 5
	for(int i=0;i<kernelLength;i++) {
		kernelValues.push_back(1.0/kernelLength);
	}
}

vector<Vector2f> Assignment7::GenerateStreamLines(Vector2f startPoint, int pointNumber) {
	vector<Vector2f> line;
	line.push_back(startPoint);
	bool done = false;
	std::vector<Vector2f>::iterator it;
	int steps = 0;
	
	//Forwards interpolation:
	while(!done && steps < pointNumber) {
		it = line.end();
		Vector2f pos = line.back();

		Vector2f newPos = RungeKuttaIntegration(pos,true);
		//Discard conditions: out of bounds, too-low or zero magnitude
		if(vField.insideBounds(newPos)) {
			line.insert(it,newPos);
			if(vField.sample(newPos).getSqrNorm()<minimumMagnitude) {
				done = true;
			}
		} else {
			done = true;
		}
		steps++;
	}
	//int stepsRemaining = RKSteps - steps;
	done = false;
	steps = 0;

	while(!done && steps < pointNumber) {
		it = line.begin();
		Vector2f pos = line.front();

		Vector2f newPos = RungeKuttaIntegration(pos,false);
		//Discard conditions: out of bounds, too-low or zero magnitude
		if(vField.insideBounds(newPos)) {
			line.insert(it,newPos);
			if(vField.sample(newPos).getSqrNorm()<minimumMagnitude) {
				done = true;
			}
		} else {
			done = true;
		}
		steps++;
	}

	return line;
}

float Assignment7::convolveKernel(int startIndex, vector<Vector2f> line) {
	//This function assumes you've used the equidistant streamline method
	int kernelBoundaries = (int)kernelLength/2;
	float result = 0;
	for(int i = 0; i<kernelLength; i++) {
		//Calculate the index for the line sample
		int lineIndex = startIndex + i - kernelBoundaries;
		//Correct for out-of-bounds
		lineIndex = lineIndex < 0 ? 0 : lineIndex;
		lineIndex = lineIndex > line.size()-1 ? line.size()-1 : lineIndex;

		//Figure out what pixel said point belongs to:
		Vector2ui pixelIndex = texture.closestNode(line[lineIndex]);
		
		//Perform the convolution and add to the result
		result += kernelValues[i]*texture.nodeScalar(pixelIndex[0],pixelIndex[1]);
	}
	return result;
}

void Assignment7::ClassicLIC() {
	
	viewer->clear();
	vector<vector<Vector2f>> foundLines;
	vector<vector<Vector2f>> foundSampleLines;
	vector<Vector2f> line;

	LICtexture = texture;

	float segmentLength = min((texture.nodePosition(0,0)-texture.nodePosition(1,0)).getSqrNorm(),(texture.nodePosition(0,0)-texture.nodePosition(0,1)).getSqrNorm());
	GenerateKernel();
	int kernelHalf = (int)kernelLength/2;

	for(int x = 0; x < texture.dims()[0]; x++) {
		for(int y = 0; y < texture.dims()[1]; y++) {
			Vector2f startPoint = texture.nodePosition(x,y);
			
			line = GenerateStreamLines(startPoint, kernelLength);
			foundLines.push_back(line);
			int startIndex = find(line.begin(),line.end(),startPoint) - line.begin();
			LICtexture.setNodeScalar(x,y,convolveKernel(startIndex,line));
			
			//viewer->setTextureGray(LICtexture.getData());
			//viewer->refresh();
		}
	}
	/*for(int i = 0; i<foundLines.size(); i++) {
		for(int j = 0; j<foundLines[i].size()-1;j++) {
				viewer->addLine(foundLines[i][j],foundLines[i][j+1], makeVector4f(0,0,1,1));
		}
	}*/
	if(ContrastEnhancement) {
		EnhanceContrast();
	} else {
		viewer->setTextureGray(LICtexture.getData());
		viewer->refresh();
	}
}

void Assignment7::FastLIC() {
	viewer->clear();
	vector<vector<Vector2f>> foundLines;
	vector<vector<Vector2f>> foundSampleLines;
	vector<Vector2f> line;

	LICtexture = texture;
	TEMPtexture= texture;
	float segmentLength = min((texture.nodePosition(0,0)-texture.nodePosition(1,0)).getSqrNorm(),(texture.nodePosition(0,0)-texture.nodePosition(0,1)).getSqrNorm());
	GenerateKernel();
	//int iterator = 0;

	for(int x = 0; x < texture.dims()[0]; x++) {
		for(int y = 0; y < texture.dims()[1]; y++) {
			//Do not calculate streamlines for visited pixels
			if(visited.nodeScalar(x,y)>0) {
				continue;
			}
			//calculate the streamline
			line = GenerateStreamLines(texture.nodePosition(x,y), RKSteps);
			foundLines.push_back(line);
		
			//move kernel along line
			for(int n = 0; n < line.size(); n++) {
				//Find out what pixel we're visiting
				Vector2ui pixel = texture.closestNode(line[n]);
				//add one to the visit-counter for this pixel
				visited.setNodeScalar(pixel[0],pixel[1],visited.nodeScalar(pixel[0],pixel[1])+1);
				//convolve the kernel for this part of the line, and add it to the total of that pixel
				LICtexture.setNodeScalar(pixel[0],pixel[1], 
					LICtexture.nodeScalar(pixel[0],pixel[1]) + convolveKernel(n,line));
				
				//TEMPtexture.setNodeScalar(pixel[0],pixel[1],LICtexture.nodeScalar(pixel[0],pixel[1])/visited.nodeScalar(pixel[0],pixel[1]));
					
				
			} 
		}
	}
	//Adjust the values
	float value;
	for(int x = 0; x < texture.dims()[0]; x++) {
		for(int y = 0; y < texture.dims()[1]; y++) {
			value = (LICtexture.nodeScalar(x,y)-texture.nodeScalar(x,y))/visited.nodeScalar(x,y);
			LICtexture.setNodeScalar(x,y,value);
			visited.setNodeScalar(x,y,0.0);
		}
	}
	if(ContrastEnhancement) {
		EnhanceContrast();
	} else {
		viewer->setTextureGray(LICtexture.getData());
		viewer->refresh();
	}
}

void Assignment7::EnhanceContrast() {
	int n = 0;
	float P = 0;
	float mu = 0;
	for(int x = 0; x < LICtexture.dims()[0];x++) {
		for(int y = 0; y < LICtexture.dims()[1];y++) {
			float p = LICtexture.nodeScalar(x,y);
			if(p>0){
				n++;
				P += p*p;
				mu += p;
			}
		}
	}
	if(n<=1) {
		return;
	}
	mu = mu/n;
	float sd = sqrt((P-n*mu*mu)/(n-1));

	float stretching = DesiredDeviation/sd;// < 1) ? DesiredDeviation/sd : 1;

	for(int x = 0; x < LICtexture.dims()[0];x++) {
		for(int y = 0; y < LICtexture.dims()[1];y++) {
			float p = LICtexture.nodeScalar(x,y);
			LICtexture.setNodeScalar(x,y,DesiredMean + stretching*(p-mu));
		}
	}

	viewer->setTextureGray(LICtexture.getData());
	viewer->refresh();

}