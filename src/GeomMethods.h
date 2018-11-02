#pragma once

#include "ofMain.h"
#include <algorithm>

#include <vector>
#include "Pt.h"

class GeomMethods
{

public:
	GeomMethods();
	~GeomMethods();

	vector<Pt> genSitePts(int, int); // set a number of site - points 
	vector<Pt> genGridPts(int, int, int); //set points inside a boundry - plane
	
	//voronoi
	vector<Pt> setSiteToGrid(vector<Pt>, vector<Pt>); //voronoi	
	vector<vector<Pt>> gridPtsInSite(vector<Pt>, vector<Pt>); //get grid points in convex hull
	
	//convexhull
	vector<vector<Pt>> initConvexHull(vector<Pt>, vector<Pt>); //initialize and set up convex hull	
	vector<Pt> genConvexHull(vector<Pt>); //convex hull
	vector<Pt> offsetedHull(vector<Pt>, int);
	
	//utility functions
	int orientationPts(Pt, Pt, Pt);
	int ptInPoly(vector<Pt>, Pt);
	int doesintx(Pt, Pt, Pt, Pt);
	Pt intxPt(Pt, Pt, Pt, Pt);
	vector<Pt> simplify(vector<Pt>);
	Pt proj(Pt, Pt, Pt);


	//occupation
	vector<vector<Pt>> spine(vector<Pt>, float, float);

	// spline
	vector<Pt> cubicSubdiv(vector<Pt>);
	Pt cubicSplit(Pt, Pt);
	vector<Pt>  cubicAvg(vector<Pt>);
	
};

