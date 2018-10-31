#pragma once

#include "Pt.h"
#include <vector>
#include <algorithm>

class RegionMethods
{
public:
	vector<Pt> hull;
	vector<vector<Pt>> segs;
	vector<vector<Pt>> bays;

	RegionMethods(vector<Pt>);
	~RegionMethods();

	void setHull(vector<Pt>);
	void setSegs(vector<vector<Pt>>);
	void setBays(vector<vector<Pt>>);

	vector<Pt> getHull();
	vector<vector<Pt>> getSegs();
	vector<vector<Pt>> getBays();
};


