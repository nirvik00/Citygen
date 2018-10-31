#pragma once

#include "ofMain.h"
#include "GeomMethods.h"
#include "Pt.h"
#include "RegionMethods.h"

#include <vector>
#include <algorithm>

using namespace std;

struct Colr {
	int re, gr, bl;
	Colr(int a, int b, int c) {
		re = a;
		gr = b;
		bl = c;
	}
};

struct Seg {
public:
	Pt a, b;
	Seg(){}
	Seg(Pt a_, Pt b_) {
		a = a_; b = b_;
	}	
	void display() {
		ofDrawLine(a.x, a.y, a.z, b.x, b.y, b.z);
	}
	void display(Seg s) {
		Pt p = s.a; Pt q = s.b;
		float e = 0;
		ofSetLineWidth(3);
		ofSetColor(0, 0, 0, 150);
		ofDrawLine(a.x, e, a.z, b.x, e, b.z);
		ofDrawLine(b.x, e, b.z, q.x, e, q.z);
		ofDrawLine(q.x, e, q.z, p.x, e, p.z);
		ofDrawLine(p.x, e, p.z, a.x, e, a.z);
	}
};

struct Quad {
	Seg A, B;
	Pt p, q, r, s;
	Quad(Seg a, Seg b) {
		A = a; B = b;
		p = A.a; q = A.b; r = B.b; s = B.a;
	}
	void display() {
		ofSetLineWidth(1);
		ofSetColor(0, 0, 0, 150);
		ofDrawLine(p.x, p.y, p.z, q.x, q.y, q.z);
		ofDrawLine(q.x, q.y, q.z, r.x, r.y, r.z);
		ofDrawLine(r.x, r.y, r.z, s.x, s.y, s.z);
		ofDrawLine(s.x, s.y, s.z, p.x, p.y, p.z);
	}
};

struct Cell {
private:
	//int FILL, HEIGHT;

public:
	Pt p, q, r, s;
	int FILL, HEIGHT;
	int I, J;		
	Cell(Pt a, Pt b, Pt c, Pt d, int i, int j) {
		p = a; q = b; r = c; s = d;
		I = i; J = j;
	}
	void display() {
		ofSetLineWidth(1);
		ofSetColor(0, 0, 0, 150);
		ofDrawLine(p.x, p.y, p.z, q.x, q.y, q.z);
		ofDrawLine(q.x, q.y, q.z, r.x, r.y, r.z);
		ofDrawLine(r.x, r.y, r.z, s.x, s.y, s.z);
		ofDrawLine(s.x, s.y, s.z, p.x, p.y, p.z);
	}
	void setHeight(int t) {
		HEIGHT = t;
	}
	void setFill(int t) {
		FILL = t;
	}
	int getHeight() {
		return HEIGHT;
	}
	int getFill() {
		return FILL;
	}
	void display2() {
		if (FILL == 1) {
			ofSetLineWidth(1);
			ofSetColor(0, 0, 0, 150);
			for (int e = 0; e < HEIGHT; e++) {
				float i = e * 10;
				ofDrawLine(p.x, i, p.z, q.x, i, q.z);
				ofDrawLine(q.x, i, q.z, r.x, i, r.z);
				ofDrawLine(r.x, i, r.z, s.x, i, s.z);
				ofDrawLine(s.x, i, s.z, p.x, i, p.z);
			}			
		}
	}
};

struct Block {
public:
	vector<Seg> seg;
	vector<Pt> hullpts;
	vector<Quad> blockquadvec;
	vector<Cell> cellvec;
	Block(vector<Seg> s, vector<Pt>p_) {
		seg = s; hullpts = p_;
	}
	void setQuadVec(vector<Quad> quads) {
		blockquadvec.clear();
		vector<Quad>().swap(blockquadvec);
		blockquadvec = quads;
	}
	void setCellVec(vector<Cell> cells) {
		cellvec.clear();
		vector<Cell>().swap(cellvec);
		cellvec = cells;
	}
};



class ofApp : public ofBaseApp {
public:
	void nsInit();
	void nsOccupy();
	void nsGenCell();
	void nsInteract();

	void setup();
	void update();
	void draw();

	void keyPressed(int key);
	void keyReleased(int key);
	void mouseMoved(int x, int y);
	void mouseDragged(int x, int y, int button);
	void mousePressed(int x, int y, int button);
	void mouseReleased(int x, int y, int button);
	void mouseEntered(int x, int y);
	void mouseExited(int x, int y);
	void windowResized(int w, int h);
	void dragEvent(ofDragInfo dragInfo);
	void gotMessage(ofMessage msg);

	ofEasyCam cam;
	
	vector<Colr> colorVec;
	vector<Pt> sitePtVec;
	vector<Pt> gridPtVec;
	vector<vector<Pt>> hullPts;
	vector<vector<Pt>> scaleHullPts;
	vector<vector<Pt>> intHullPts;
	vector<vector<Pt>> splinePts;

	vector<RegionMethods> regVec;
	vector<Seg> spinevec;
	vector<Block> blockvec;
	

	int ITERATION = 1;
	int BOARD_DIMENSION = 300;
	int NUM_SITES = 3;
	int CELL_LE = BOARD_DIMENSION / 100;
	int CELL_WI = BOARD_DIMENSION / 100;
	int SCALE_HULL = (int)NUM_SITES/5;
	float SPINE_DEPTH = 50;
	float BAY_DEPTH = 50;
	GeomMethods geomMethods;

	
};