#pragma once

#include "ofMain.h"
#include "GeomMethods.h"
#include "Pt.h"

#include <vector>
#include <algorithm>

using namespace std;

struct Colr {
	int re, gr, bl;
	Colr() {}
	Colr(int a, int b, int c) {
		re = a;
		gr = b;
		bl = c;
	}
	void display() {
		std::cout << "--- Color check : " << re << ", " << gr << ", " << bl << endl;
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
	Quad(Pt a, Pt b, Pt c, Pt d) {
		p = a; q = b; r = c; s = d;
	}
	void display() {
		ofSetLineWidth(1);
		ofSetColor(255, 0, 0);
		p.y = -5; q.y = -5; r.y = -5, s.y = -5;
		ofDrawLine(p.x, p.y, p.z, q.x, q.y, q.z);
		
		//ofDrawLine(q.x, q.y, q.z, r.x, r.y, r.z);

		ofSetColor(0, 0, 255);
		ofSetLineWidth(5);
		ofDrawLine(r.x, r.y, r.z, s.x, s.y, s.z);
		//ofDrawLine(s.x, s.y, s.z, p.x, p.y, p.z);

		ofSetColor(0, 0, 255);
		ofDrawSphere(p.x, p.y, p.z, 5);
		ofSetColor(0, 255, 0);
		ofDrawSphere(q.x, -10, q.z, 5);
	}
};

struct Cell {
private:
	//int FILL, HEIGHT;

public:
	Colr colr;
	Pt p, q, r, s;
	int CELL_FILL = 0, CELL_HEIGHT = 0.0;
	float CELL_ANGLE=0;
	int I, J, Num;		
	Cell() {}
	Cell(Pt a, Pt b, Pt c, Pt d, int j, int i, int n, float cellangle, Colr cl) {
		p = a; q = b; r = c; s = d;
		I = i; J = j; Num = n;
		CELL_ANGLE = cellangle;
	}
	void display() {
		ofSetLineWidth(1);
		ofSetColor(0, 0, 0, 150);
		ofDrawLine(p.x, p.y, p.z, q.x, q.y, q.z);
		ofDrawLine(q.x, q.y, q.z, r.x, r.y, r.z);
		ofDrawLine(r.x, r.y, r.z, s.x, s.y, s.z);
		ofDrawLine(s.x, s.y, s.z, p.x, p.y, p.z);
	}
	void setMinAngle(float t) {
		CELL_ANGLE = t;
	}
	float getMinAngle(float t) {
		return CELL_ANGLE;
	}
	void setHeight(int t) {
		CELL_HEIGHT = t;
	}
	void setFill(int t) {
		CELL_FILL = t;
	}
	int getHeight() {
		return CELL_HEIGHT;
	}
	int getFill() {
		return CELL_FILL;
	}
	void display2() {
		if (CELL_FILL == 1) {
			for (int e = 1; e < CELL_HEIGHT; e++) {
				float i = e * 10;
				ofSetLineWidth(1);
				ofSetColor(0, 0, 0);
				ofMesh mesh;
				ofPoint p0(p.x, (e - 1) * 10, p.z); ofPoint p1(p.x, e * 10, p.z);
				ofPoint q0(q.x, (e - 1) * 10, q.z); ofPoint q1(q.x, e * 10, q.z);
				ofPoint r0(r.x, (e - 1) * 10, r.z); ofPoint r1(r.x, e * 10, r.z);
				ofPoint s0(s.x, (e - 1) * 10, s.z); ofPoint s1(s.x, e * 10, s.z);

				if (e > 0 && (e < CELL_HEIGHT - 1)) {
					// avoid boolean operations by eliminating the horizontal lines
					ofDrawLine(p0, q0); ofDrawLine(q0, r0); ofDrawLine(r0, s0); ofDrawLine(s0, p0);//bottom horizontal
					ofDrawLine(p1, q1); ofDrawLine(q1, r1); ofDrawLine(r1, s1); ofDrawLine(s1, p1);//top horizontal
					ofSetLineWidth(2);
					ofDrawLine(p0, p1); ofDrawLine(q0, q1); ofDrawLine(r0, r1); ofDrawLine(s0, s1);//vertical 
					
					//front faces
					mesh.addVertex(p0); mesh.addVertex(q1); mesh.addVertex(q0);
					mesh.addVertex(p0); mesh.addVertex(q1); mesh.addVertex(p1);

					mesh.addVertex(q0); mesh.addVertex(r1); mesh.addVertex(r0);
					mesh.addVertex(q0); mesh.addVertex(r1); mesh.addVertex(q1);

					mesh.addVertex(r0); mesh.addVertex(s1); mesh.addVertex(r1);
					mesh.addVertex(r0); mesh.addVertex(s1); mesh.addVertex(s0);

					mesh.addVertex(p0); mesh.addVertex(s1); mesh.addVertex(p1);
					mesh.addVertex(p0); mesh.addVertex(s1); mesh.addVertex(s0);
					
					//bottom face
					mesh.addVertex(p0); mesh.addVertex(q0); mesh.addVertex(r0);
					mesh.addVertex(p0); mesh.addVertex(s0); mesh.addVertex(r0);

					//top face 
					mesh.addVertex(p1); mesh.addVertex(q1); mesh.addVertex(r1);
					mesh.addVertex(p1); mesh.addVertex(s1); mesh.addVertex(r1);
					
				}
				else {
					ofSetLineWidth(2);
					//show the roof line for clarity
					//ofDrawLine(p0, q0); ofDrawLine(q0, r0); ofDrawLine(r0, s0); ofDrawLine(s0, p0);//bottom horizontal
				}
				mesh.setupIndicesAuto();
				//ofSetColor(colr.re, colr.gr, colr.bl);
				ofSetColor(150, 150, 150);
				mesh.draw();
			}	
		}
		else {
			ofMesh mesh;
			ofPoint p0(p.x, 0, p.z);
			ofPoint q0(q.x, 0, q.z);
			ofPoint r0(r.x, 0, r.z);
			ofPoint s0(s.x, 0, s.z);
			//bottom face
			mesh.addVertex(p0); mesh.addVertex(q0); mesh.addVertex(r0);
			mesh.addVertex(p0); mesh.addVertex(s0); mesh.addVertex(r0);
			mesh.setupIndicesAuto();
			ofSetColor(0, 250, 0);
			mesh.draw();
		}
	}
	void display3() {
		p.y = -10;
		q.y = -10;
		r.y = -10;
		s.y = -10;
		if (J < 3) {
			ofSetLineWidth(1);
			ofSetColor(255, 0, 0);
			ofDrawLine(p.x, p.y, p.z, r.x, r.y, r.z);
			ofDrawLine(q.x, q.y, q.z, s.x, s.y, s.z);
		}
		else if (J > Num - 3) {
			ofSetLineWidth(1);
			ofSetColor(0, 0, 255);
			ofDrawLine(p.x, p.y, p.z, r.x, r.y, r.z);
			ofDrawLine(q.x, q.y, q.z, s.x, s.y, s.z);
		}
		else {
			ofSetLineWidth(1);
			ofSetColor(0, 255, 0);
			ofDrawLine(p.x, p.y, p.z, r.x, r.y, r.z);
			ofDrawLine(q.x, q.y, q.z, s.x, s.y, s.z);
		}
	}
};

struct Block {
public:
	Colr colr;
	vector<Seg> seg;
	vector<Pt> hullpts;
	vector<Quad> blockquadvec;
	vector<Cell> cellvec;
	vector<vector<Cell>> rowcellvec; //vector of cells in a row
	float areaOccupied;
	Block(vector<Seg> s, vector<Pt>p_, Colr cl) {
		seg = s; hullpts = p_; colr = cl;
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
	// NS METHODS
	void nsInit();
	void nsOccupy();
	void nsGenCell();
	void nsRules();
	//curve skeleton
	void nsgenCrvSkeleton(); //generate curve skeleton from block attributes:1.segments 2.hull
	float nscomputeArea();


	// OF METHODS
	void setup();
	void update();
	void draw();

	void keyPressed(int key);

	// UNUSED OF METHODS
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
	ofLight light;
	
	vector<Colr> colorVec;
	vector<Pt> sitePtVec;
	vector<Pt> gridPtVec;
	vector<vector<Pt>> hullPts;
	vector<vector<Pt>> scaledHullPts;
	vector<vector<Pt>> intHullPts;
	vector<vector<Pt>> splinePts;


	vector<vector<Pt>> secHullPts;

	vector<Seg> spinevec;
	vector<Block> blockvec;
	
	int ITERATION = 0;
	// ratio of BOARD_DIMENSION : NUM_SITES = 500 : 5; 1500 : 45
	int BOARD_DIMENSION = 1000;
	int NUM_SITES = 5;
	int CELL_LE = BOARD_DIMENSION / 100;
	int CELL_WI = BOARD_DIMENSION / 100;
	int SCALE_HULL = 10; //20 // max(10, (int)NUM_SITES / 15);
	float CELL_LENGTH = 50;
	float BAY_DEPTH = 50;
	float MIN_CELL_ANGLE = PI/10;
	GeomMethods geomMethods;
	float REQ_FAR = 2.5;
	
};	