#include "Pt.h"


Pt::Pt(float x_, float y_, float z_)
{
	x = x_;
	y = y_;
	z = z_;
}

Pt::~Pt()
{
}


float Pt::Dis(Pt b) {
	float dx = (x - b.x);
	float dy = (y - b.y);
	float dz = (z - b.z);
	return sqrt(dx*dx + dy*dy + dz*dz);
}

void Pt::setId(int t) {
	siteId = t;
}

int Pt::getId() {
	return siteId;
}
