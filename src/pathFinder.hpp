#pragma once

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/opencv.hpp"

#include "../include/MST.h"
#include "../include/Matching/Matching.h"
#include "../include/Matching/Graph.h"
#include "../include/Christofides.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <math.h>

namespace cv{
namespace PathFinder {
	
bool pointInsideSegment(double x1, double x2, double y1, double y2, double x, double y) 
{ 
	
    if ((x >= x1) && (x < x2) && (y >= y1) && (y < y2)){
	return true;
	}
	else{
		return false; 
	}
} 

typedef std::map<int, Point2f> idPointMap;

class pathFinder{

private:
std::vector<Point2f> path;

public:
std::vector<Point2f> getPath(){return path;};

pathFinder(std::vector<Point2f> points, Mat stippledImage, int horSeg, int vertSeg){
	
Size size=stippledImage.size();
std::vector<idPointMap> imageSegmentation((horSeg*vertSeg));

int imageSegmentationPosition;

// The next loop distributes the stippled points into several regions (the segmentation of the image depends on horSeg and vertSeg).
// For each region an approximate optimal path by means of Christofides algorithm will be found. Also for each segment it includes one integer value for each point in the region which
//it will be used eventually as ID when running the Christofides algorithm.

for(int i=0; i<points.size();i++){
	bool assigned=false;
	imageSegmentationPosition=0;
	int horAux=0;
	int vertAux;
	
	while(horAux<horSeg && assigned==false){
	vertAux=0;

		while(vertAux<vertSeg && assigned==false){

			double horSegmentMin=horAux*(size.width/horSeg);
			double horSegmentMax=(horAux+1)*(size.width/horSeg);
			double vertSegmentMin=vertAux*(size.height/vertSeg);
			double vertSegmentMax=(vertAux+1)*(size.height/vertSeg);
			
//This lines of code prevent that since we are assigning the segment to the points doing intervals of the type [  )[  )...[  ) horizontally and the same rotated 90° clockwise vertically,
// the points lying exactly on the final edge not being considerated, hence we add 1 in order to also "catch" those points lying on the right or lower edge (so the last segement is "of type []" ).
			if((horAux+1)==horSeg){
				horSegmentMax=horSegmentMax+1;
			}
			if((vertAux+1)==vertSeg){
				vertSegmentMax=vertSegmentMax+1;
			}

			Point2f pointAux=points[i];
						
			assigned=pointInsideSegment(horSegmentMin, horSegmentMax, vertSegmentMin, vertSegmentMax, pointAux.x, pointAux.y);
			
			if(assigned){
			//idValueSegment gives an ID to the point within the map corresponding to the segment
				int idValueSegment;
				idPointMap* idPointMapSegmentAux= &(imageSegmentation[imageSegmentationPosition]);
//If the corresponding segment is not empty and there is already other points, it adds the new one
				if(idPointMapSegmentAux->rbegin() != idPointMapSegmentAux->rend()){
					idValueSegment = (idPointMapSegmentAux->rbegin()->first) + 1;
					idPointMapSegmentAux->insert({idValueSegment, pointAux});
				}
				else{
					idValueSegment=0;
					idPointMapSegmentAux->insert({idValueSegment, pointAux});
				}
			}
			vertAux++;
			imageSegmentationPosition++;
		}
		horAux++;
	}
}//End of the segmentation loop


//Here a vector of integers segmentTour is created. It controls in which order the image segments will be tour. This code tours the image as a "snake" from up to down and from left to right, so something like ununu...

std::vector<int> segmentTour;
int alternator=1;
int position=1;
segmentTour.push_back(position-1);
int multiplicator=1;
int counter=1;

while(segmentTour.size()<vertSeg*horSeg){
while(counter<vertSeg){
	position=position+1*alternator;
	segmentTour.push_back(position-1);
	counter++;
}
alternator=alternator*(-1);
counter=1;
position=position + vertSeg;

if (!(segmentTour.size()==vertSeg*horSeg)) {
segmentTour.push_back(position-1);
}

}




//Begin of the search of the pseudooptimal path
//For each element in the vector image segmentation, which contains a list of points with an associated id, runs the christofides algorithm
Graph G;
vector<double> cost;

std::cout << "Image Segmentation size: " << imageSegmentation.size() << std::endl;
std::cout << "\n " << std::endl;

for(int p=0; p<segmentTour.size(); p++){

int k= segmentTour[p];
 // std::cout << "Por aqui imageSegmentation[k] size"  << imageSegmentation[k].size() << std::endl;

	cost.clear();
		
//For the segment, it runs the Christofides algorithm in case there is points within that region
 if(imageSegmentation[k].size() > 2){
	 
	G = Graph(imageSegmentation[k].size());
	
idPointMap idPointMapAux = imageSegmentation[k];

		for(int i = 0; i < imageSegmentation[k].size(); i++)
			for(int j = i+1; j < imageSegmentation[k].size(); j++)
				G.AddEdge(i, j);

	for(int i = 0; i < G.GetNumEdges(); i++)
	{
	    pair<int, int> p = G.GetEdge(i);
	    int u = p.first, v = p.second;
		Point2f pointAuxFirst = idPointMapAux[u];
		Point2f pointAuxSecond = idPointMapAux[v];
	    cost.push_back( sqrt( pow(pointAuxFirst.x-pointAuxSecond.x, 2) + pow(pointAuxFirst.y-pointAuxSecond.y, 2) ) );
	}
		
		
//The results of the christofides algorithm used is a little bit strange, giving as output edges numbering and (non ordered) nodes numbers,
//so this part of the code just take that output and creates a simple ordered list of points which are the solution to that segement and 
//it adds this ordered list of points to the total solution list of points "path".
	pair< vector<int> , double > p = Christofides(G, cost);
	vector<int> solIDs = p.first;
	std::cout << "Segment number: " << k << endl;
	 std::cout << "Number of points in the segment: " << imageSegmentation[k].size() << std::endl;
	std::cout << "Number of edges in the solution: " << (solIDs.size()-1) << std::endl;
	std::cout << "Solution cost: " << p.second << endl;

	int previousEdge1, previousEdge2;
	
	previousEdge1 = G.GetEdge(solIDs[0]).first;
	previousEdge2 = G.GetEdge(solIDs[0]).second;
	
	for(int i = 0; i < (int)solIDs.size()-1; i++){
		Point2f solutionPathPointAux;
		if(i==0){
			i++;
			if(previousEdge1 == G.GetEdge(solIDs[i]).first){
				solutionPathPointAux = idPointMapAux[previousEdge2];
				path.push_back(solutionPathPointAux);
				std::cout << "Pushed back point: " << solutionPathPointAux.x << " " <<  solutionPathPointAux.y << std::endl;
				solutionPathPointAux = idPointMapAux[previousEdge1];
				path.push_back(solutionPathPointAux);
				std::cout << "Pushed back point: " << solutionPathPointAux.x << " " <<  solutionPathPointAux.y << std::endl;
				solutionPathPointAux = idPointMapAux[G.GetEdge(solIDs[i]).second];
				path.push_back(solutionPathPointAux);
				std::cout << "Pushed back point: " << solutionPathPointAux.x << " " <<  solutionPathPointAux.y << std::endl;
			}
			else if(previousEdge1 == G.GetEdge(solIDs[i]).second){
				solutionPathPointAux = idPointMapAux[previousEdge2];
				path.push_back(solutionPathPointAux);
				std::cout << "Pushed back point: " << solutionPathPointAux.x << " " <<  solutionPathPointAux.y << std::endl;
				solutionPathPointAux = idPointMapAux[previousEdge1];
				path.push_back(solutionPathPointAux);
				std::cout << "Pushed back point: " << solutionPathPointAux.x << " " <<  solutionPathPointAux.y << std::endl;
				solutionPathPointAux = idPointMapAux[G.GetEdge(solIDs[i]).first];
				path.push_back(solutionPathPointAux);
				std::cout << "Pushed back point: " << solutionPathPointAux.x << " " <<  solutionPathPointAux.y << std::endl;
			}
			else if (previousEdge2 == G.GetEdge(solIDs[i]).first){
				solutionPathPointAux = idPointMapAux[previousEdge1];
				path.push_back(solutionPathPointAux);
				std::cout << "Pushed back point: " << solutionPathPointAux.x << " " <<  solutionPathPointAux.y << std::endl;
				solutionPathPointAux = idPointMapAux[previousEdge2];
				path.push_back(solutionPathPointAux);
				std::cout << "Pushed back point: " << solutionPathPointAux.x << " " <<  solutionPathPointAux.y << std::endl;
				solutionPathPointAux = idPointMapAux[G.GetEdge(solIDs[i]).second];
				path.push_back(solutionPathPointAux);
				std::cout << "Pushed back point: " << solutionPathPointAux.x << " " <<  solutionPathPointAux.y << std::endl;
			}
			else {
				solutionPathPointAux = idPointMapAux[previousEdge1];
				path.push_back(solutionPathPointAux);
				std::cout << "Pushed back point: " << solutionPathPointAux.x << " " <<  solutionPathPointAux.y << std::endl;
				solutionPathPointAux = idPointMapAux[previousEdge2];
				path.push_back(solutionPathPointAux);
				std::cout << "Pushed back point: " << solutionPathPointAux.x << " " <<  solutionPathPointAux.y << std::endl;
				solutionPathPointAux = idPointMapAux[G.GetEdge(solIDs[i]).first];
				path.push_back(solutionPathPointAux);
				std::cout << "Pushed back point: " << solutionPathPointAux.x << " " <<  solutionPathPointAux.y << std::endl;
			}
			previousEdge1=G.GetEdge(solIDs[i]).first;
			previousEdge2=G.GetEdge(solIDs[i]).second;
		}
		else{
			if(previousEdge1 == G.GetEdge(solIDs[i]).first){
				solutionPathPointAux = idPointMapAux[G.GetEdge(solIDs[i]).second];
				path.push_back(solutionPathPointAux);
				std::cout << "Pushed back point: " << solutionPathPointAux.x << " " <<  solutionPathPointAux.y << std::endl;
			}
			else if(previousEdge2== G.GetEdge(solIDs[i]).first){
				solutionPathPointAux = idPointMapAux[G.GetEdge(solIDs[i]).second];
				path.push_back(solutionPathPointAux);
				std::cout << "Pushed back point: " << solutionPathPointAux.x << " " <<  solutionPathPointAux.y << std::endl;
			}
			else{
				solutionPathPointAux = idPointMapAux[G.GetEdge(solIDs[i]).first];
				path.push_back(solutionPathPointAux);
				std::cout << "Pushed back point: " << solutionPathPointAux.x << " " <<  solutionPathPointAux.y << std::endl;
			}
			previousEdge1=G.GetEdge(solIDs[i]).first;
			previousEdge2=G.GetEdge(solIDs[i]).second;
		}		
	}
	std::cout << "\n " << endl;
 }//End of the if for segments containing two points or more within them
 else if (imageSegmentation[k].size() == 1){
	Point2f solutionPathPointAux;
	idPointMap idPointMapAux = imageSegmentation[k];
	solutionPathPointAux = idPointMapAux[0];
	path.push_back(solutionPathPointAux);
	
	std::cout << "Segment number: " << k << endl;
	std::cout << "Single point in the segment. Point added: " << solutionPathPointAux.x << " " << solutionPathPointAux.y << endl;
	std::cout << "\n" << endl;

}//End of the else if in case there is only one point in the segment
 else if (imageSegmentation[k].size() == 2){
	Point2f solutionPathPointAux1;
	Point2f solutionPathPointAux2;
	idPointMap idPointMapAux = imageSegmentation[k];
	solutionPathPointAux1 = idPointMapAux[0];
	solutionPathPointAux2 = idPointMapAux[1];
	path.push_back(solutionPathPointAux1);
	path.push_back(solutionPathPointAux2);
	
	std::cout << "Segment number: " << k << endl;
	std::cout << "Two points in the segment. First point added: " << solutionPathPointAux1.x << " " << solutionPathPointAux1.y << endl;
	std::cout << "Two points in the segment. Second point added: " << solutionPathPointAux2.x << " " << solutionPathPointAux2.y << endl;
	std::cout << "\n" << endl;

}//End of the else if in case there is only one point in the segment

}//End of Christofides algorithm for each segment

}

}; //End of the class

} //End of pathFinder namespace
} //End of cv namespace
