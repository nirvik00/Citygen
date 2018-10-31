#include "ofApp.h"

void ofApp::nsInit() {
	ITERATION++;
	std::cout << "1->2 Command received at nsINIT" << endl;
	colorVec.clear();
	sitePtVec.clear();
	gridPtVec.clear();
	hullPts.clear();
	scaleHullPts.clear();
	spinevec.clear();
	blockvec.clear();
	intHullPts.clear();

	vector<Colr>().swap(colorVec);
	vector<Pt>().swap(sitePtVec);
	vector<Pt>().swap(gridPtVec);
	vector<vector<Pt>>().swap(hullPts);
	vector<vector<Pt>>().swap(scaleHullPts);
	vector<Seg>().swap(spinevec);
	vector<Block>().swap(blockvec);
	vector<vector<Pt>>().swap(intHullPts);


	for (int i = 0; i < NUM_SITES; i++) {
		int re = (int)ofRandom(0, 255);
		int gr = (int)ofRandom(0, 255);
		int bl = (int)ofRandom(0, 255);
		Colr col(re, gr, bl);
		colorVec.push_back(col);
	}
		
	sitePtVec = geomMethods.genSitePts(NUM_SITES, BOARD_DIMENSION);
	vector<Pt>tempGridPtVec = geomMethods.genGridPts(CELL_LE, CELL_WI, BOARD_DIMENSION);
	gridPtVec = geomMethods.setSiteToGrid(sitePtVec, tempGridPtVec);
		
	
	hullPts = geomMethods.initConvexHull(sitePtVec, gridPtVec);
	cout << "\t REGIONS constructed...going to next step" << endl;

	//scale hull down	
	for (int i = 0; i < hullPts.size(); i++) {
		vector<Pt>temp = hullPts[i];
		vector<Pt> newPts;
		for (int j = 0; j < temp.size(); j++) {
			Pt a, b;
			if (j == 0) {
				a = temp[temp.size() - 1];
				b = temp[0];
			}
			else {
				a = temp[j - 1];
				b = temp[j];
			}			
			Pt c((a.x + b.x) / 2, (a.y + b.y) / 2, (a.z + b.z) / 2);
			float d = a.Dis(b);
			Pt u((b.x - a.x) / d, (b.y - a.y) / d, (b.z - a.z) / d);
			float sc = SCALE_HULL;
			Pt v(c.x - u.z*sc, c.y + u.y, c.z + u.x*sc);
			if (geomMethods.ptInPoly(temp, v) == 1) {
				newPts.push_back(v);
			}			
		}
		vector<Pt> tempHullPts = geomMethods.genConvexHull(newPts);
		scaleHullPts.push_back(tempHullPts);
	}
	
	//scaleHullPts = hullPts;
	cout << "\t Hulls constructed...going to next step" << endl;
	nsOccupy();
}

void ofApp::nsOccupy() {
	cout << "2/2-> Command received at nsOCCUPY...2 objectives here " << endl;
	for (int i = 0; i < scaleHullPts.size(); i++) {
		vector<Pt> temp = scaleHullPts[i];
		Pt spA, spB; float maxD = -1;
		for (int j = 0; j < temp.size(); j++) {
			for (int k = 0; k < temp.size(); k++) {
				float d = temp[j].Dis(temp[k]);
				if (d > maxD) {
					maxD = d;
					spA = temp[j];
					spB = temp[k];
				}
			}
		}
		spinevec.push_back(Seg(spA, spB));
	}
	cout << "\t Spine constructed" << endl;
	for (int i = 0; i < spinevec.size(); i++) {
		vector<Seg>segvec; vector<Seg> segvecL; vector<Seg> segvecR;
		vector<Pt>hullpts = scaleHullPts[i];

		Pt A = spinevec[i].a; Pt B = spinevec[i].b;

		float D = A.Dis(B);
		Pt u((B.x - A.x) / D, (B.y - A.y) / D, (B.z - A.z) / D);
		Pt v(-u.z, u.y, u.x); Pt w(u.z, u.y, -u.x);

		for (int j = 0; j < D; j += BAY_DEPTH) {
			Pt p(A.x + u.x*j, A.y + u.y*j, A.z + u.z*j);
			float sc = 1000; float sc1 = 30;
			Pt q(p.x + v.x*sc, p.y + v.y*sc, p.z + v.z*sc);
			Pt r(p.x + w.x*sc, p.y + w.y*sc, p.z + w.z*sc);

			Pt q_(p.x + v.x*sc1, p.y + v.y*sc1, p.z + v.z*sc1);

			Pt I(-1, -1, -1); Pt J(-1, -1, -1);
			for (int k = 0; k < hullpts.size(); k++) {
				Pt h0, h1;
				if (k == 0) {
					h0 = hullpts[hullpts.size() - 1];
					h1 = hullpts[0];
				}
				else {
					h0 = hullpts[k - 1];
					h1 = hullpts[k];
				}
				Pt intx = geomMethods.intxPt(h0, h1, p, q);
				Pt intx2 = geomMethods.intxPt(h0, h1, p, r);
				if (intx.x != -1 && intx.y != -1 && intx.z != -1) {
					I = intx;
				}
				if (intx2.x != -1 && intx2.y != -1 && intx2.z != -1) {
					J = intx2;
				}
			}

			int T0 = (I.x != -1 && I.y != -1 && I.z != -1 && I.Dis(p) > SPINE_DEPTH);
			int T1 = (J.x != -1 && J.y != -1 && J.z != -1 && J.Dis(p) > SPINE_DEPTH);
			
			/*
			if (T0 == 1) {
				segvecL.push_back(Seg(p, I));
			}
			if (T1 == 1) {
				segvecR.push_back(Seg(p, J));
			}
			*/
			if (T0 == 1 && T1 == 1) { segvec.push_back(Seg(I, J)); }
			else if (T0 == 0 && T1 == 1) { segvec.push_back(Seg(p, J)); }
			else if (T0 == 1 && T1 == 0) { segvec.push_back(Seg(I, p)); }
		}
		blockvec.push_back(Block(segvec, hullpts));
	}
	cout << "\t Bays constructed" << endl;
	for (int i = 0; i < blockvec.size(); i++) {
		vector<Quad> qv;
		vector<Seg>segvec = blockvec[i].seg;
		for (int j = 1; j < segvec.size(); j++) {
			Seg s0, s1;
			s0 = segvec[j - 1];
			s1 = segvec[j];
			Quad quad(s0, s1);
			qv.push_back(quad);
		}
		blockvec[i].setQuadVec(qv);
	}
	cout << "Objective complete at nsOCCUPY\nITERATION NUMBER= " << ITERATION << "\n\n\n" << endl;
	nsGenCell();
}

void ofApp::nsGenCell() {
	for (int i = 0; i < blockvec.size(); i++) {
		vector<Quad>qvec = blockvec[i].blockquadvec;
		vector<Cell>cellvec;
		for (int j = 0; j < qvec.size(); j++) {
			qvec[j].display();
			vector<Pt>bottom; vector<Pt>top;
			Pt p = qvec[j].p; Pt q = qvec[j].q; Pt r = qvec[j].r; Pt s = qvec[j].s;
			int itr = 0;
			for (float k = 0.f; k < 1.0; k += 0.1) {
				Pt a(p.x + (q.x - p.x)*k, p.y + (q.y - p.y)*k, p.z + (q.z - p.z)*k);
				Pt b(s.x + (r.x - s.x)*k, s.y + (r.y - s.y)*k, s.z + (r.z - s.z)*k);
				//ofDrawLine(a.x, a.y, a.z, b.x, b.y, b.z);
				bottom.push_back(a); top.push_back(b);
				Pt p0, p1, p2, p3;
				if (itr > 0) {
					p0 = bottom[itr - 1]; p1 = bottom[itr]; p2 = top[itr]; p3 = top[itr - 1];
					cellvec.push_back(Cell(p0, p1, p2, p3, itr, j));
				}
				itr++;
			}
		}
		blockvec[i].setCellVec(cellvec);
	}
	cout << "Objective complete at nsGENCELL" << endl;
	nsInteract();
}

void ofApp::nsInteract() {
	//interact cells and display updated  form
	for (int i = 0; i < blockvec.size(); i++) {
		vector<Cell>::iterator j = blockvec[i].cellvec.begin();
		while(j!= blockvec[i].cellvec.end()){
			int J = (*j).J;
			int ht, fill;
			if (J < 3 || J>7) {
				int R = ofRandom(0, 10);
				if (R > 3) {
					fill = 1;
				}
				else {
					fill = 0;
				}
				int S = ofRandom(0, 10);
				if (S < 4) {
					ht = 2;
				}
				else if (S > 3 && S < 5) {
					ht = 5;
				}
				else {
					ht = 10;
				}
				(*j).HEIGHT=ht;
				(*j).FILL=fill;
			}	
			j++;
		}
	}
	//cout << "Objective complete at nsINTERACTION" << endl;
}

void ofApp::setup(){
	ofEnableDepthTest();
	ofSetVerticalSync(true);
	cam.lookAt(ofVec3f(-500, 0, -500));
	cam.setDistance(1000);
	cam.setPosition(ofVec3f(1250, 250, 1250));
	nsInit();
}

void ofApp::update(){	
	//nsInit();
}

void ofApp::draw(){
	cam.begin();
	
	ofSetBackgroundColor(255, 255, 255); ofFill(); ofDrawAxis(50);	

	for (int i = 0; i < hullPts.size(); i++) {
		vector<Pt> hull = hullPts[i];
		for (int j = 0; j < hull.size(); j++) {
			Pt a, b;
			if (j == 0) {
				a = hull[hull.size() - 1];
				b = hull[j];
			}
			else {
				a = hull[j - 1];
				b = hull[j];
			}
			ofSetLineWidth(1);
			ofSetColor(255, 100, 0, 100);
			ofDrawLine(a.x, a.y, a.z, b.x, b.y, b.z);
		}
	}

	for (int i = 0; i < scaleHullPts.size(); i++) {
		vector<Pt> hull = scaleHullPts[i];
		for (int j = 0; j < hull.size(); j++) {
			Pt a, b;
			if (j == 0) {
				a = hull[hull.size() - 1];
				b = hull[j];
			}
			else {
				a = hull[j - 1];
				b = hull[j];
			}
			ofSetLineWidth(1);
			ofSetColor(255, 0, 255, 100);
			ofDrawLine(a.x, a.y, a.z, b.x, b.y, b.z);
		}
	}

	//display CElls
	for (int i = 0; i < blockvec.size(); i++) {
		vector<Cell> cellvec = blockvec[i].cellvec;
		for (int j = 0; j < cellvec.size(); j++) {
			cellvec[j].display();
			cellvec[j].display2();
		}
	}

	ofSetLineWidth(1);
	cam.end();
}

void ofApp::keyPressed(int key) {
	if (key == 'r') { 		
		cout << "\n\n\nUSER Command: "<< ITERATION <<", for new board" << endl;
		nsInit(); 
	}

	if (key == 'i') {
		cout << "\n\n INTERACTion command" << endl;
		nsInteract();
	}

	if (key == 's') {
		for (int i = 0; i < blockvec.size(); i++) {
			vector<Cell> cellvec = blockvec[i].cellvec;
			for (int j = 0; j < cellvec.size(); j++) {
				cout << " height" << cellvec[j].getHeight() << " Fill: " << cellvec[j].getFill() << endl;
			}
		}
	}
}

void ofApp::keyReleased(int key){

}

void ofApp::mouseMoved(int x, int y ){

}

void ofApp::mouseDragged(int x, int y, int button){

}

void ofApp::mousePressed(int x, int y, int button){

}

void ofApp::mouseReleased(int x, int y, int button){

}

void ofApp::mouseEntered(int x, int y){

}

void ofApp::mouseExited(int x, int y){
}

void ofApp::windowResized(int w, int h){

}

void ofApp::gotMessage(ofMessage msg){

}

void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
