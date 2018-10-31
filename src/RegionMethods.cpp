#include "RegionMethods.h"



RegionMethods::RegionMethods(vector<Pt> p){
	hull.clear();
	for (int i = 0; i < p.size(); i++) {
		hull.push_back(p[i]);
	}
}

RegionMethods::~RegionMethods(){}

void RegionMethods::setHull(vector<Pt> p) {
	hull.clear();
	for (int i = 0; i < p.size(); i++) {
		hull.push_back(p[i]);
	}
}

void RegionMethods::setBays(vector<vector<Pt>> p) {
	bays.clear();
	for (int i = 0; i < p.size(); i++) {
		bays.push_back(p[i]);
	}
}

void RegionMethods::setSegs(vector<vector<Pt>> p) {
	segs.clear();
	for (int i = 0; i < p.size(); i++) {
		segs.push_back(p[i]);
	}
}

vector<Pt> RegionMethods::getHull() {
	return hull;
}

vector<vector<Pt>> RegionMethods::getSegs() {
	return segs;
}

vector<vector<Pt>> RegionMethods::getBays() {
	return bays;
}