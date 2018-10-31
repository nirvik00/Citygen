#pragma once

#include "ofMain.h"
#include <cmath>

class Pt
{
private:
	int siteId;
public:
	float x, y, z;
	int lock;
	Pt() {}
	Pt(float, float, float);
	~Pt();
	
	float Dis(Pt);

	void setId(int t);
	int getId();
};

