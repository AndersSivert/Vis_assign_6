//---------------------------------------------------------------------------
#include "stdafx.h"
//---------------------------------------------------------------------------
#include "Assignment6.h"
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

IMPLEMENT_GEOX_CLASS( Assignment6, 0)
{
    BEGIN_CLASS_INIT( Assignment6 );
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

	
	ADD_NOARGS_METHOD(Assignment6::ReadFieldFromFile)
	ADD_NOARGS_METHOD(Assignment6::GenerateTexture)
	ADD_NOARGS_METHOD(Assignment6::ClassicLIC)
	ADD_NOARGS_METHOD(Assignment6::FastLIC)
	ADD_NOARGS_METHOD(Assignment6::EnhanceContrast)
	ADD_NOARGS_METHOD(Assignment6::DrawStreamLines)

	
	

    
}

QWidget* Assignment6::createViewer()
{
    viewer = new GLGeometryViewer();
    return viewer;
}

Assignment6::Assignment6()
{
    viewer = NULL;
    RungeKutta = false;

	fileName = "C:\\Users\\Martin\\Desktop\\GeoX\\Assignment05\\Data\\ANoise2CT4.am";					//Martin
	//fileName = "C:\\Program Files\\GeoX\\experiments\\Visualization\\Assignment6\\Data\\ANoise2CT4.am";	//Anders
	//fileName = "";																					//Jim

	randomPoints=false;
	startingPoints="10";
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

Assignment6::~Assignment6() {}

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

Vector2f Assignment6::EulerIntegration(Vector2f pos) {
	return pos + field.sample(pos[0],pos[1])*EulerStepSize;
}

Vector2f Assignment6::EulerIntegration(Vector2f pos, bool forwards) {
	float dt;
	if(forwards) {
		dt = EulerStepSize;
	} else {
		dt = -EulerStepSize;
	}
	if(directionField){
		Vector2f sampled =field.sample(pos[0],pos[1]);
		sampled.normalize();
		return pos + sampled*dt;
	} else {
		return pos + field.sample(pos[0],pos[1])*dt;
	}
}

Vector2f Assignment6::RungeKuttaIntegration(Vector2f pos) {
	Vector2f p1 =  field.sample(pos[0],pos[1]);
	Vector2f p2 =  field.sample(pos[0]+p1[0]*RKStepSize/2,pos[1]+p1[1]*RKStepSize/2);
	Vector2f p3 =  field.sample(pos[0]+p2[0]*RKStepSize/2,pos[1]+p2[1]*RKStepSize/2);
	Vector2f p4 =  field.sample(pos[0]+p3[0]*RKStepSize,pos[1]+p3[1]*RKStepSize);
	return pos + (p1 + p2*2 + p3*2 + p4)*RKStepSize/6;
}
Vector2f Assignment6::RungeKuttaIntegration(Vector2f pos, bool forwards) {
	float dt;
	if(forwards) {
		dt = RKStepSize;
	} else {
		dt = -RKStepSize;
	}

	if(directionField) {
		Vector2f p1 =  field.sample(pos[0],pos[1]);
		p1.normalize();
		Vector2f p2 =  field.sample(pos[0]+p1[0]*dt/2,pos[1]+p1[1]*dt/2);
		p2.normalize();
		Vector2f p3 =  field.sample(pos[0]+p2[0]*dt/2,pos[1]+p2[1]*dt/2);
		p3.normalize();
		Vector2f p4 =  field.sample(pos[0]+p3[0]*dt,pos[1]+p3[1]*dt);
		p4.normalize();
		return pos + (p1 + p2*2 + p3*2 + p4)*dt/6;
	} else {
		Vector2f p1 =  field.sample(pos[0],pos[1]);
		Vector2f p2 =  field.sample(pos[0]+p1[0]*dt/2,pos[1]+p1[1]*dt/2);
		Vector2f p3 =  field.sample(pos[0]+p2[0]*dt/2,pos[1]+p2[1]*dt/2);
		Vector2f p4 =  field.sample(pos[0]+p3[0]*dt,pos[1]+p3[1]*dt);
		return pos + (p1 + p2*2 + p3*2 + p4)*dt/6;
	}
}

void Assignment6::ReadFieldFromFile(){
	if (!field.load(fileName))
    {
        output << "Error loading field file " << fileName << "\n";
        return;
    } else {
		readField = true;
	}

}

void Assignment6::DrawStreamLines()
{
	//viewer->clear();
	if(!readField){
		ReadFieldFromFile();
	}
	
	float xMin = field.boundMin()[0];
	float yMin = field.boundMin()[1];
	float xMax = field.boundMax()[0];
	float yMax = field.boundMax()[1];
	
	//Step 1: determine seed points
	vector<Vector2f> startPoints;
	if(randomPoints) {
		int N = stoi(startingPoints);
		for(int i=0;i<N;i++) {
			Vector2f newPos;
			newPos[0] = (float)rand()*(xMax-xMin)/RAND_MAX + xMin;
			newPos[1] = (float)rand()*(yMax-yMin)/RAND_MAX + yMin;
			startPoints.push_back(newPos);
		}
	} else if(grid) {
		string temp;
		stringstream ss(startingPoints);
		vector<float> gridDimensions;
		while(getline(ss,temp,',')){
			gridDimensions.push_back(stof(temp));
		}
		int xDim;
		int yDim;
		if(gridDimensions.size()==1) { 
			xDim = gridDimensions[0]; 
			yDim = gridDimensions[0];
		} else {
			xDim = gridDimensions[0]; 
			yDim = gridDimensions[1];
		}
		if(xDim==1&&yDim==1) {
			Vector2f newPos;
			newPos[0] = 0.5*(xMax+xMin);
			newPos[1] = 0.5*(yMax+yMin);
			startPoints.push_back(newPos);
		} else {
			for(int x=1; x<=xDim;x++) {
				for(int y=1; y<=yDim;y++) {
					Vector2f newPos;
					newPos[0] = x*(xMax-xMin)/(xDim+1) + xMin;
					newPos[1] = y*(yMax-yMin)/(yDim+1) + yMin;
					startPoints.push_back(newPos);
				}
			}
		}
	} else {
		Vector2f newPos;
		string temp;
		stringstream ss(startingPoints);
		if(getline(ss,temp,',')) {
			newPos[0] = stof(temp);
		} else {
			newPos[0] = 0.5*(xMin+xMax);
		}
		if(getline(ss,temp,',')) {
			newPos[1] = stof(temp);
		} else {
			newPos[1] = 0.5*(yMin+yMax);
		}
				
		startPoints.push_back(newPos);
	}
	
	if(adjustSteps>0 && maxAdjustment>0 && adjustArea > 0) {
		//float moveDistance = RKStepSize*max((xMax-xMin)/field.dims()[0],(yMax-yMin)/field.dims()[1]);
		//float moveDistance = 2*maxAdjustment/(adjustSteps*(adjustSteps+1));
		for(int n=0;n<adjustSteps;n++) {
			float stepLength = 2*maxAdjustment*(adjustSteps-n)/(adjustSteps*(adjustSteps+1));
			for(int i=0;i<startPoints.size();i++) {
				float ownMagn = field.sample(startPoints[i]).getSqrNorm();

				//lower-bounds node index:
				int xIndex = (int)(startPoints[i][0]-xMin)/(xMax-xMin);
				int xLower = (xIndex - adjustArea + 1 < 0) 
							? 0 : xIndex - adjustArea + 1;
				int xUpper = (xIndex + adjustArea > field.dims()[0]-1) 
							? field.dims()[0]-1 : xIndex + adjustArea;
				int yIndex = (int)(startPoints[i][1]-yMin)/(yMax-yMin);
				int yLower = (yIndex - adjustArea + 1 < 0) 
							? 0 : yIndex - adjustArea + 1;
				int yUpper = (yIndex + adjustArea > field.dims()[1]-1) 
							? field.dims()[1]-1 : yIndex + adjustArea;
				Vector2f movement;
				movement.setZero();
				
				//For all nodes in the cell
				for(int x = xLower; x <= xUpper; x++) {
					for(int y = yLower; y <= yUpper; y++) {
						//Find direction vector from position to node
						Vector2f temp = field.nodePosition(x,y)
										- startPoints[i];
						
						//temp[0] = x-0.5;
						//temp[1] = y-0.5;
						float dist = temp.getSqrNorm();
						dist = (dist == 0) ? 1 : dist;
						temp.normalize();
						
						//weigh it with the difference between the magnitude in the
						// point and the magnitude in the node. Move towards higher
						//magnitudes and away from lower.
						//output << "x: " << x << ", y " << y << "\n";
						movement += temp*(field.node(x,y).getSqrNorm()
										  - ownMagn)/dist;
					}
				}
				//output << "movement vector " << movement << "\n";
				//Now we have the weighted sum. Normalize this vector...
				movement.normalize();
				//And set its length to the above calculated one
				//movement = movement * moveDistance;//*(1.1*adjustSteps-n)/adjustSteps;
				movement = movement * stepLength;
				//Check if we're inside the bounds:
				if(field.insideBounds(startPoints[i]+movement)) {
					//If we are, move the point
					startPoints[i] = startPoints[i] + movement;
				}
			}
			
		}
	}
	
	vector<vector<Vector2f>> AllStreamLines;
	for(int i=0;i<startPoints.size();i++) {
		vector<Vector2f> StreamLine;
		StreamLine.push_back(startPoints[i]);
		bool dontSwap = false;
		bool forwards = true;
		int stepsTaken = 0;
		float length = 0;
		bool done = false;
		std::vector<Vector2f>::iterator it;

		while(!done && stepsTaken<RKSteps && length<maxLength) {
			Vector2f pos;
			if(forwards){
				pos = StreamLine.back();
				it = StreamLine.end();
			} else {
				pos = StreamLine.front();
				it = StreamLine.begin();
			}
			Vector2f newPos = RungeKuttaIntegration(pos,forwards);
			//Discard conditions: out of bounds, too-low or zero magnitude
			bool stopDirection = false;
			if(newPos[0]<xMin || newPos[0]>xMax || newPos[1]<yMin || newPos[1]>yMax) {
					stopDirection = true;
			} else {
				StreamLine.insert(it,newPos);
			
				//Too slow:
				if(field.sample(newPos).getSqrNorm() < minimumMagnitude) {
					stopDirection = true;
				}
				length += (newPos - pos).getSqrNorm();
				stepsTaken+=1;
			}
			if(dontSwap) {
				if(stopDirection) {
					done = true;
				} 
			} else {
				forwards = !forwards;
				if(stopDirection) {
					dontSwap = true;
				}
			} 

		}
		AllStreamLines.push_back(StreamLine);
	}
	if(magnitudeColor) {
		float maxMagnitude = 0;

		for(int i=0;i<field.dims()[0];i++) {
			for(int j=0;j<field.dims()[1];j++) {
				maxMagnitude = max(maxMagnitude,field.node(i,j).getSqrNorm());
			}
		}
		if (maxMagnitude == 0) {
			output << "all nodes are of magnitude 0" << "\n";
			maxMagnitude = 1;
		} 
		float c;
		float s = 0.3;
		float c1;
		float c2;
		float c3;
		for(int i=0;i<AllStreamLines.size();i++) {
			for(int n=1;n<AllStreamLines[i].size();n++) {
				c = 0.5*field.sample(AllStreamLines[i][n-1]).getSqrNorm();
				c += 0.5*field.sample(AllStreamLines[i][n]).getSqrNorm();
				c = c/maxMagnitude;
		
				c1 = exp(- c*c/(2*s*s))/(s*sqrt(2*3.1415));
				c2 = 0.9*exp(- (c-0.5)*(c-0.5)/(2*(s-0.05)*(s-0.05)))/(s*sqrt(2*3.1415));
				c3 = exp(- (c-1)*(c-1)/(2*(s-0.05)*(s-0.05)))/(s*sqrt(2*3.1415));
				viewer->addLine(AllStreamLines[i][n-1],AllStreamLines[i][n],makeVector4f(c3,c2,c1,1));
			}
		}
	} else {
		for(int i=0;i<AllStreamLines.size();i++) {
			for(int n=1;n<AllStreamLines[i].size();n++) {
				viewer->addLine(AllStreamLines[i][n-1],AllStreamLines[i][n],makeVector4f(1,1,1,1));
			}
		}
	}
	if(showPoints) {
		for(int i=0;i<AllStreamLines.size();i++) {
			for(int n=1;n<AllStreamLines[i].size();n++) {
				Point2D newPoint(startPoints[i]);
				newPoint.color = makeVector4f(0,1,0,1);
				viewer->addPoint(newPoint);
				viewer->addLine(AllStreamLines[i][n-1],AllStreamLines[i][n],makeVector4f(1,1,1,1));
			}
		}
		for(int i=0;i<startPoints.size();i++) {
			Point2D newPoint(startPoints[i]);
			newPoint.color = makeVector4f(1,0,0,1);
			viewer->addPoint(newPoint);
		}
		
	}
	viewer->refresh();
}

void Assignment6::GenerateTexture() {
	viewer->clear();
	texture.clear();
	visited.clear();
	int iWidth = NextPOT(field.boundMax()[0]-field.boundMin()[0]);
	if(xPowerOfTwo>0) {
		iWidth = pow(2.0,xPowerOfTwo);
	}
	int iHeight = NextPOT(field.boundMax()[1]-field.boundMin()[1]);
	if(yPowerOfTwo>0) {
		iHeight = pow(2.0,yPowerOfTwo);
	}

	srand(textureSeed);
	
	texture.init(field.boundMin(),field.boundMax(),makeVector2ui(iWidth,iHeight));
	LICtexture.init(field.boundMin(),field.boundMax(),makeVector2ui(iWidth,iHeight));
	visited.init(field.boundMin(),field.boundMax(),makeVector2ui(iWidth,iHeight));
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

void Assignment6::GenerateKernel() {
	kernelValues.clear();
	//TODO: Implement different kernels, kernel sizes etc.
	//Currently box kernel of size 5
	for(int i=0;i<kernelLength;i++) {
		kernelValues.push_back(1.0/kernelLength);
	}
}

vector<Vector2f> Assignment6::GenerateStreamLines(Vector2f startPoint, int pointNumber) {
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
		if(field.insideBounds(newPos)) {
			line.insert(it,newPos);
			if(field.sample(newPos).getSqrNorm()<minimumMagnitude) {
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
		if(field.insideBounds(newPos)) {
			line.insert(it,newPos);
			if(field.sample(newPos).getSqrNorm()<minimumMagnitude) {
				done = true;
			}
		} else {
			done = true;
		}
		steps++;
	}

	return line;
}

vector<Vector2f> Assignment6::GenerateStreamLineEquidistant(Vector2f startPoint, float segmentLength) {
	
	vector<Vector2f> sampleLine;
	sampleLine.push_back(startPoint);
	bool done = false;
	std::vector<Vector2f>::iterator it;
	float length = 0;
	float maxLength = kernelLength*segmentLength/1.9;
	float permanentLength = 0;
	int steps = 0;
	Vector2f pos = startPoint;
	//output << "new line \n";
	//Forwards interpolation:
	while(!done && steps < RKSteps && permanentLength<maxLength) {
		
		Vector2f newPos = RungeKuttaIntegration(pos,true);
		//output << "forwards length " << length << "\n";
		permanentLength += (newPos-pos).getSqrNorm();
		//Discard conditions: out of bounds, too-low or zero magnitude
		if(field.insideBounds(newPos)) {
			if(field.sample(newPos).getSqrNorm()<minimumMagnitude) {
				done = true;
			}
			
			//If we're passing the arc-length threshold, add a
			//point at the proper coordinates. the while is to handle
			//the stepping length being longer than the segment length
			int n = 1;
			while((newPos-pos).getSqrNorm()+length>=segmentLength) {
				it = sampleLine.end();
				Vector2f direction = newPos - pos;
				direction.normalize();
				pos = pos + direction*(segmentLength - length);
				//Linear interpolation of the point?
				sampleLine.insert(it,pos);
				//Reset the length to the next point
				length = 0;
			}
			//update the length we've traveled
			length += (newPos-pos).getSqrNorm();
			//store the new position
			pos = newPos;
		} else {
			done = true;
		}
		steps++;
	}
	done = false;
	maxLength = kernelLength*segmentLength;
	length = 0;
	permanentLength=0;
	
	pos = startPoint;

	while(!done && steps < RKSteps && permanentLength<maxLength) {
		
		Vector2f newPos = RungeKuttaIntegration(pos,false);
		//output << "forwards length " << length << "\n";
		permanentLength += (newPos-pos).getSqrNorm();
		//Discard conditions: out of bounds, too-low or zero magnitude
		if(field.insideBounds(newPos)) {
			if(field.sample(newPos).getSqrNorm()<minimumMagnitude) {
				done = true;
			}
			
			//If we're passing the arc-length threshold, add a
			//point at the proper coordinates. the while is to handle
			//the stepping length being longer than the segment length
			while((newPos-pos).getSqrNorm()+length>=segmentLength) {
				it = sampleLine.begin();
				Vector2f direction = newPos - pos;
				direction.normalize();
				pos = pos + direction*(segmentLength - length);
				//Linear interpolation of the point?
				sampleLine.insert(it,pos);
				//Reset the length to the next point
				length = 0;
			}
			//update the length we've traveled
			length += (newPos-pos).getSqrNorm();
			//store the new position
			pos = newPos;
		} else {
			done = true;
		}
		steps++;
	}
	//output << "line has length " << line.size() << "\n";

	return sampleLine;
}

vector<Vector2f> Assignment6::GenerateStreamLineEquidistantLong(Vector2f startPoint, float segmentLength) {
	
	vector<Vector2f> sampleLine;
	sampleLine.push_back(startPoint);
	bool done = false;
	std::vector<Vector2f>::iterator it;
	float length = 0;
	int steps = 0;
	int points = 0;
	Vector2f pos = startPoint;

	//Forwards interpolation:
	while(!done && steps < RKSteps && points < RKSteps) {
		
		Vector2f newPos = RungeKuttaIntegration(pos,true);
		
		//Discard conditions: out of bounds, too-low or zero magnitude
		if(field.insideBounds(newPos)) {
			if(field.sample(newPos).getSqrNorm()<minimumMagnitude) {
				done = true;
			}
			
			//If we're passing the arc-length threshold, add a
			//point at the proper coordinates. the while is to handle
			//the stepping length being longer than the segment length
			while((newPos-pos).getSqrNorm()+length>=segmentLength && points < RKSteps) {
				it = sampleLine.end();
				Vector2f direction = newPos - pos;
				direction.normalize();
				//Linear interpolation of the point
				pos = pos + direction*(segmentLength - length);
				sampleLine.insert(it,pos);
				//Reset the length to the next point
				length = 0;
				points += 1;
			}
			//update the length we've traveled
			length += (newPos-pos).getSqrNorm();
			//store the new position
			pos = newPos;
		} else {
			done = true;
		}
		steps++;
	}
	done = false;
	length = 0;
	points = 0;
	
	pos = startPoint;

	while(!done && steps < RKSteps && points < RKSteps) {
		
		Vector2f newPos = RungeKuttaIntegration(pos,false);
		//Discard conditions: out of bounds, too-low or zero magnitude
		if(field.insideBounds(newPos)) {
			if(field.sample(newPos).getSqrNorm()<minimumMagnitude) {
				done = true;
			}
			
			//If we're passing the arc-length threshold, add a
			//point at the proper coordinates. the while is to handle
			//the stepping length being longer than the segment length
			while((newPos-pos).getSqrNorm()+length>=segmentLength && points < RKSteps) {
				it = sampleLine.begin();
				Vector2f direction = newPos - pos;
				direction.normalize();
				//Linear interpolation of the point
				pos = pos + direction*(segmentLength - length);
				sampleLine.insert(it,pos);
				//Reset the length to the next point
				length = 0;
				points += 1;
			}
			//update the length we've traveled
			length += (newPos-pos).getSqrNorm();
			//store the new position
			pos = newPos;
		} else {
			done = true;
		}
		steps++;
	}
	//output << "line has length " << line.size() << "\n";

	return sampleLine;
}


float Assignment6::convolveKernel(int startIndex, vector<Vector2f> line) {
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

void Assignment6::ClassicLIC() {
	
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

void Assignment6::FastLIC() {
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

void Assignment6::EnhanceContrast() {
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