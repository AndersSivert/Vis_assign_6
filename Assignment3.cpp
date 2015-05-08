//---------------------------------------------------------------------------
#include "stdafx.h"
//---------------------------------------------------------------------------
#include "Assignment3.h"
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

IMPLEMENT_GEOX_CLASS( Assignment3, 0)
{
    BEGIN_CLASS_INIT( Assignment3 );

    ADD_SEPARATOR("Scalar field")
    ADD_STRING_PROP(ScalarfieldFilename, 0)
    ADD_NOARGS_METHOD(Assignment3::ReadData)
    ADD_NOARGS_METHOD(Assignment3::DrawData)
	ADD_BOOLEAN_PROP(DrawGrid,0);
	ADD_BOOLEAN_PROP(DrawPoints,0);
	
    ADD_SEPARATOR("Isocontour settings")
	ADD_BOOLEAN_PROP(DrawIsoContours,0)
	ADD_BOOLEAN_PROP(AutoIsoValues,0)
    ADD_STRING_PROP(manualIsoValues, 0)
	ADD_INT32_PROP(numberOfIsocontours,0)
	ADD_BOOLEAN_PROP(AsymptoticDecider,0)
	

    
}

QWidget* Assignment3::createViewer()
{
    viewer = new GLGeometryViewer();
    return viewer;
}

Assignment3::Assignment3()
{
    viewer = NULL;
    ScalarfieldFilename = "C:\\Users\\Martin\\Desktop\\GeoX\\Assignment03\\";
	field.setZero();
	DrawGrid = true;
	DrawPoints = false;
	DrawIsoContours = false;
	AsymptoticDecider = true;

	AutoIsoValues = false;
	manualIsoValues = "";
	numberOfIsocontours = 0;
}

Assignment3::~Assignment3() {}


void Assignment3::ReadData()
{
	

    if (!field.load(ScalarfieldFilename))
    {
        output << "Error loading field file " << ScalarfieldFilename << "\n";
        return;
    }
    ////Get the minimum/maximum value in that field
    //float32 min = std::numeric_limits<float32>::max();
    //float32 max = -std::numeric_limits<float32>::max();
    //for(size_t j=0; j<field.dims()[1]; j++)
    //{
    //    for(size_t i=0; i< field.dims()[0]; i++)
    //    {
    //        const float32 val = field.nodeScalar(i,j);
    //        min = val < min ? val : min;
    //        max = val > max ? val : max;
    //    }
    //}

    //Draw a point for each grid vertex.
    //for(size_t j=0; j<field.dims()[1]; j++)
    //{
    //    for(size_t i=0; i<field.dims()[0]; i++)
    //    {
    //        const float32 val = field.nodeScalar(i, j);
    //        const float32 c = (val - min) / (max - min);

    //        Point2D p;
    //        p.position  = field.nodePosition(i, j);
    //        p.size = 5;
    //        //Use a grayscale depending on the actual value
    //        p.color[0] = c; p.color[1] = c; p.color[2] = c;
    //        viewer->addPoint(p);
    //    }
    //}

}

void Assignment3::DrawData()
{
	
    viewer->clear();

	int32 yDim = field.dims()[0];
	int32 xDim = field.dims()[1];

	//Get the minimum/maximum value in that field
		float32 min = std::numeric_limits<float32>::max();
		float32 max = -std::numeric_limits<float32>::max();
		for(size_t j=0; j<field.dims()[1]; j++)
		{
		    for(size_t i=0; i< field.dims()[0]; i++)
		    {
		        const float32 val = field.nodeScalar(i,j);
		        min = val < min ? val : min;
		        max = val > max ? val : max;
		    }
		}

	if(DrawGrid)
	{
		for(size_t j=0; j<field.dims()[1];j++)
		{
			viewer->addLine(field.nodePosition(0,j),field.nodePosition(xDim-1,j),makeVector4f(1.0,1.0,1.0,1.0));
		}
		for(size_t i=0; i<field.dims()[0];i++)
		{
			viewer->addLine(field.nodePosition(i,0),field.nodePosition(i,yDim-1),makeVector4f(1.0,1.0,1.0,1.0));
		}
	}

	if(DrawPoints) {

		//Draw a point for each grid vertex.
		for(size_t j=0; j<field.dims()[1]; j++)
		{
		    for(size_t i=0; i<field.dims()[0]; i++)
		    {
		        const float32 val = field.nodeScalar(i, j);
		        const float32 c = (val - min) / (max - min);
				
				float s = 0.3;
				float c1 = exp(- c*c/(2*s*s))/(s*sqrt(2*3.1415));
				float c2 = exp(- (c-0.5)*(c-0.5)/(2*s*s))/(s*sqrt(2*3.1415));
				float c3 = exp(- (c-1)*(c-1)/(2*s*s))/(s*sqrt(2*3.1415));

		        Point2D p;
		        p.position  = field.nodePosition(i, j);
		        p.size = 5;
				
		        //Color things
		        //p.color[0] = c; p.color[1] = c*(1 - c); p.color[2] = 1-c;
				p.color[0] = c1; p.color[1] = c2; p.color[2] = c3;
		        viewer->addPoint(p);
		    }
		}
	}
	if(DrawIsoContours)
	{

		isoValues.clear();
		isoColours.clear();
		if(AutoIsoValues) {
			float interval = 1/(numberOfIsocontours+1);
			for(int i=1; i<=numberOfIsocontours;i++)
			{
				isoValues.push_back((i*(max-min)/(numberOfIsocontours+1)) + min);
			}
		} else {
			int n = 0;
			string temp;
			stringstream ss(manualIsoValues);
			while(getline(ss,temp,',')) {
				isoValues.push_back(stof(temp));
			}
			
		}
		for(int i=0;i<(int)isoValues.size();i++) {
			float c = (float)i/(float)isoValues.size();
			float s = 0.3;
			float c1 = exp(- c*c/(2*s*s))/(s*sqrt(2*3.1415));
			float c2 = exp(- (c-0.5)*(c-0.5)/(2*s*s))/(s*sqrt(2*3.1415));
			float c3 = exp(- (c-1)*(c-1)/(2*s*s))/(s*sqrt(2*3.1415));
			isoColours.push_back(makeVector4f(c1,c2,c3,1));
		}
		if(AsymptoticDecider)
		{
			DrawIsoAsympt();
		} else {
			DrawIsoMid();
		}
	}

    viewer->refresh();
}

void Assignment3::DrawIsoAsympt()
{
	vector< Point2D > CrossingPoints;
	for(size_t j=0; j<field.dims()[1]-1; j++)
		{
			for(size_t i=0; i< field.dims()[0]-1; i++)
		    {
				float values [4] = {field.nodeScalar(i, j),field.nodeScalar(i+1, j),
					field.nodeScalar(i+1, j+1),field.nodeScalar(i, j+1)};
				for(int k = 0; k<(int)isoValues.size();k++)
				{
					CrossingPoints.clear();
					//If all the scalar values in this grid are above or below the
					//iso value, skip evaluating for that value.
					if((isoValues[k] > *std::max_element(values, values+4))
						|| (isoValues[k] < *std::min_element(values,values+4)))
					{
						continue;
					}
					//Order of point checking
					int order [5][2];
				
					order[0][0] = i;   order[0][1] = j;
					order[1][0] = i+1; order[1][1] = j;
					order[2][0] = i+1; order[2][1] = j+1;
					order[3][0] = i;   order[3][1] = j+1;
					order[4][0] = i;   order[4][1] = j;
					
					
					for(int l = 0;l<4;l++)
					{
						float f0 = field.nodeScalar(order[l][0],order[l][1]);
						float f1 = field.nodeScalar(order[l+1][0],order[l+1][1]);
						if((f0 >= isoValues[k])	== (f1 >= isoValues[k]))
						{
							continue;
						}
						Point2D p;
						p.position = (field.nodePosition(order[l+1][0],order[l+1][1])
									 - field.nodePosition(order[l][0],order[l][1])).operator*
									((isoValues[k] - f0) / (f1 - f0));
						//We have now calculated the direction and length of the vector that goes from the first point to
						//the crossing point.
						p.position += field.nodePosition(order[l][0],order[l][1]);
						//Now we have moved the vector to its proper position.
						CrossingPoints.push_back(p);
						p.position[0];
					}
					if((int)CrossingPoints.size()==0)
					{
						continue;
					}
					std::sort(CrossingPoints.begin(),CrossingPoints.end(),PointSorter);
					for(int l=0; l < (int)CrossingPoints.size()-1; l += 2)
					{
						viewer->addLine(CrossingPoints[l].position,
										CrossingPoints[l+1].position,
										isoColours[k]);
					}
					
				
				}
			}
	}
}

void Assignment3::DrawIsoMid()
{
	
	//To get the isolines, we can simply connect the isoline-edge
	//intersections in the order we find them, going clockwise.
	//This should solve all cases; in the unambigous cases there's
	//only a single line, otherwise we just need to start at the
	//correct sign. also, in the ambigous case, if the first node
	//is the wrong sign, the next one in a clockwise direction
	//has to be the right sign, otherwise it wouldn't be an
	//ambigous case.

	vector< Point2D > CrossingPoints;
	for(size_t j=0; j<field.dims()[1]-1; j++)
		{
			for(size_t i=0; i< field.dims()[0]-1; i++)
		    {
				float values [4] = {field.nodeScalar(i, j),field.nodeScalar(i+1, j),
					field.nodeScalar(i+1, j+1),field.nodeScalar(i, j+1)};
				
				for(int k = 0; k<(int)isoValues.size();k++)
				{
					CrossingPoints.clear();
					//If all the scalar values in this grid are above or below the
					//iso value, skip evaluating for that value.
					if((isoValues[k] > *std::max_element(values, values+4))
						|| (isoValues[k] < *std::min_element(values,values+4)))
					{
						continue;
					}
					//Order of point checking
					int order [5][2];
					
					if((isoValues[k] >= values[0])
						!= (isoValues[k] >= 0.25*(values[0]+values[1]+values[2]+values[3])))
					{//if (i,j) is on the same side of the isovalue as the midpoint,
					//start there
						order[0][0] = i;   order[0][1] = j;
						order[1][0] = i+1; order[1][1] = j;
						order[2][0] = i+1; order[2][1] = j+1;
						order[3][0] = i;   order[3][1] = j+1;
						order[4][0] = i;   order[4][1] = j;
					} else
					{//Otherwise, start at the next point in the order of permutation.
					//If it's an unambigous case, where we start doesn't matter,
					//and if it is we're starting at the correct sign.
						order[0][0] = i+1; order[0][1] = j;
						order[1][0] = i+1; order[1][1] = j+1;
						order[2][0] = i;   order[2][1] = j+1;
						order[3][0] = i;   order[3][1] = j;
						order[4][0] = i+1; order[4][1] = j;

					}
					
					for(int l = 0;l<4;l++)
					{
						float f0 = field.nodeScalar(order[l][0],order[l][1]);
						float f1 = field.nodeScalar(order[l+1][0],order[l+1][1]);
						if((f0 >= isoValues[k])	== (f1 >= isoValues[k]))
						{
							continue;
						}
						Point2D p;
						p.position = (field.nodePosition(order[l+1][0],order[l+1][1])
									 - field.nodePosition(order[l][0],order[l][1])).operator*
									((isoValues[k] - f0) / (f1 - f0));
						//We have now calculated the direction and length of the vector that goes from the first point to
						//the crossing point.
						p.position += field.nodePosition(order[l][0],order[l][1]);
						//Now we have moved the vector to its proper position.
						CrossingPoints.push_back(p);
					}
					if((int)CrossingPoints.size()==0)
					{
						continue;
					}
					for(int l=0;l < (int)CrossingPoints.size()-1;l+=2)
					{
						viewer->addLine(CrossingPoints[l].position,
										CrossingPoints[l+1].position,
										isoColours[k]);
					}
					/*
					Point2D p1;
					Point2D p2;
					Point2D p3;
					
					p3.position= (p1.position + p2.position);*/
				
				}
			}
	}
		
}




bool Assignment3::PointSorter(const Point2D & p1,const Point2D & p2)
{
	return p1.position[0]<p2.position[0];
}