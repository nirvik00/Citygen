#include "ofApp.h"


/* ns methods  */
void ofApp::nsInit() {
	ITERATION++;
	std::cout << "1->2 Command received at nsINIT" << endl;
	colorVec.clear();
	sitePtVec.clear();
	gridPtVec.clear();
	hullPts.clear();
	scaledHullPts.clear();
	spinevec.clear();
	blockvec.clear();
	intHullPts.clear();

	vector<Colr>().swap(colorVec);
	vector<Pt>().swap(sitePtVec);
	vector<Pt>().swap(gridPtVec);
	vector<vector<Pt>>().swap(hullPts);
	vector<vector<Pt>>().swap(scaledHullPts);
	vector<Seg>().swap(spinevec);
	vector<Block>().swap(blockvec);
	vector<vector<Pt>>().swap(intHullPts);


	for (int i = 0; i < NUM_SITES; i++) {
		int re = (int)ofRandom(0, 255);
		int gr = (int)ofRandom(0, 255);
		int bl = (int)ofRandom(0, 255);
		Colr col(re, gr, bl);
		col.display();
		colorVec.push_back(col);
	}
		
	sitePtVec = geomMethods.genSitePts(NUM_SITES, BOARD_DIMENSION);
	vector<Pt>tempGridPtVec = geomMethods.genGridPts(CELL_LE, CELL_WI, BOARD_DIMENSION);
	gridPtVec = geomMethods.setSiteToGrid(sitePtVec, tempGridPtVec);
		
	
	hullPts = geomMethods.initConvexHull(sitePtVec, gridPtVec);
	std::cout << "\t REGIONS constructed...going to next step" << endl;

	scaledHullPts = hullPts;
	std::cout << "\t Hulls constructed...going to next step" << endl;
	nsOccupy();
	nsgenCrvSkeleton();
	nsGenCell();
	nsRules();
	nscomputeArea();
}

void ofApp::nsOccupy() {
	vector<vector<Pt>>().swap(secHullPts);
	std::cout << "2/2-> Command received at nsOCCUPY...2 objectives here " << endl;
	for (int i = 0; i < scaledHullPts.size(); i++) {
		vector<Pt> temp = scaledHullPts[i];
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
		vector<Pt>().swap(temp);
		spinevec.push_back(Seg(spA, spB));		
	}
	std::cout << "\t Spine constructed" << endl;

	for (int i = 0; i < spinevec.size(); i++) {
		vector<Seg>segvec; vector<Seg> segvecL; vector<Seg> segvecR;
		vector<Pt>hullpts = scaledHullPts[i];

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

			int T0 = (I.x != -1 && I.y != -1 && I.z != -1 && I.Dis(p) > CELL_LENGTH);
			int T1 = (J.x != -1 && J.y != -1 && J.z != -1 && J.Dis(p) > CELL_LENGTH);

			if (T0 == 1 && T1 == 1) { segvec.push_back(Seg(I, J));  }
			else if (T0 == 0 && T1 == 1) { segvec.push_back(Seg(p,J)); }
			else if (T0 == 1 && T1 == 0) { segvec.push_back(Seg(I, p)); }
		}
		blockvec.push_back(Block(segvec, hullpts, colorVec[i]));
	}
	std::cout << "\t Segments constructed @ nsOccupy" << endl;
}

void ofApp::nsgenCrvSkeleton() {
	for (int i = 0; i < blockvec.size(); i++) {
		vector<Quad> qv;
		vector<Pt>csSpine;
		vector<Pt> Le; vector<Pt> Ri;
		vector<Seg> seg = blockvec[i].seg;
		if (seg.size() < 3) { continue; }
		for (int j = 0; j < seg.size(); j++) {
			Pt a = seg[j].a; Le.push_back(a);
			Pt b = seg[j].b; Ri.push_back(b);
			Pt c((a.x + b.x) / 2, (a.y + b.y) / 2, (a.z + b.z) / 2);
			ofSetColor(0);
			csSpine.push_back(c);
		}
		vector<Seg>xseg; vector<Seg>fseg;
		for (int j = 1; j < Le.size() - 1; j++) {
			Pt a = Le[j - 1];
			Pt b = Le[j];
			Pt c = Le[j + 1];
			Pt u((c.x - a.x) / a.Dis(c), (c.y - a.y) / a.Dis(c), (c.z - a.z) / a.Dis(c));
			Pt v(-u.z, u.y, u.x); Pt w(u.z, u.y, -u.x);
			float sc = 10000;
			Pt p(b.x + v.x*sc, b.y + v.y*sc, b.z + v.z*sc);
			Pt q(b.x + w.x*sc, b.y + w.y*sc, b.z + w.z*sc);
			for (int k = 1; k < csSpine.size(); k++) {
				Pt m, n;
				m = csSpine[k - 1];
				n = csSpine[k];
				Pt I = geomMethods.intxPt(p, b, m, n);
				Pt J = geomMethods.intxPt(q, b, m, n);
				if (I.x != -1 && I.y != -1 && I.z != -1) {
					xseg.push_back(Seg(b, I));
				}
				if (J.x != -1 && J.y != -1 && J.z != -1) {
					xseg.push_back(Seg(b, J));
				}
			}
		}
		for (int j = 1; j < Ri.size() - 1; j++) {
			Pt a = Ri[j - 1];
			Pt b = Ri[j]; 
			Pt c = Ri[j + 1];
			Pt u((c.x - a.x) / a.Dis(c), (c.y - a.y) / a.Dis(c), (c.z - a.z) / a.Dis(c));
			Pt v(-u.z, u.y, u.x); Pt w(u.z, u.y, -u.x);
			float sc = 10000;
			Pt p(b.x + v.x*sc, b.y + v.y*sc, b.z + v.z*sc);
			Pt q(b.x + w.x*sc, b.y + w.y*sc, b.z + w.z*sc);
			for (int k = 1; k < csSpine.size(); k++) {
				Pt m, n;
				m = csSpine[k - 1];
				n = csSpine[k];
				Pt I = geomMethods.intxPt(p, b, m, n);
				Pt J = geomMethods.intxPt(q, b, m, n);
				if (I.x != -1 && I.y != -1 && I.z != -1) {
					xseg.push_back(Seg(b, I));
				}
				if (J.x != -1 && J.y != -1 && J.z != -1) {
					xseg.push_back(Seg(b, J));
				}
			}
		}
		for (int j = 0; j < xseg.size(); j++) {
			//xseg.a=point on hull; xseg.b= spine point; 
			Pt a = xseg[j].a; Pt b = xseg[j].b;
			float minD = BOARD_DIMENSION * 2; int sum = 0;
			Pt m(-1, -1, -1); Pt n(-1, -1, -1);
			//if this seg is cut by LONGER segment ignore it.
			for (int k = 0; k < xseg.size(); k++) {
				Pt p = xseg[k].a; Pt q = xseg[k].b;
				if (a.Dis(b) < p.Dis(q)) { continue; }
				if (a.Dis(p) < 1 && b.Dis(q)<1) { continue; }
				Pt I = geomMethods.intxPt(a, b, p, q);
				if (I.x != -1 && I.y != -1 && I.z != -1) {
					float d = I.Dis(a);
					if (d < minD && d>1) {
						m = a; n = I; minD = d;
						sum++;
					}
				}
			}
			if (sum == 0) {
				fseg.push_back(Seg(a, b));
			}
			else {
				xseg[j].a = m; xseg[j].b = n;
				fseg.push_back(Seg(m, n));
			}
		}
		vector<Seg> gseg;
		for (int j = 0; j < fseg.size(); j++) {
			Pt a = xseg[j].a; Pt b = xseg[j].b;
			float minD = BOARD_DIMENSION * 2; int sum = 0;
			Pt m(-1, -1, -1); Pt n(-1, -1, -1);
			//dont ignore any intersection
			for (int k = 0; k < fseg.size(); k++) {
				Pt p = fseg[k].a; Pt q = fseg[k].b;
				if (a.Dis(p) < 1 && b.Dis(q) < 1) { continue; }
				Pt I = geomMethods.intxPt(a, b, p, q);
				if (I.x != -1 && I.y != -1 && I.z != -1) {
					float d = I.Dis(a);
					if (d < minD && d>1) {
						m = a; n = I; minD = d;
						sum++;
					}
				}
			}
			if (sum == 0) {
				gseg.push_back(Seg(a, b));
			}
			else {
				fseg[j].a = m; fseg[j].b = n;
				gseg.push_back(Seg(m, n));
			}
		}

		for (int j = 1; j < gseg.size() - 1; j++) {
			//gseg.a=point on hull; fseg.b= spine point; 
			Pt a, b, c, d;
			a = gseg[j - 1].a; b = gseg[j - 1].b;
			d = gseg[j].a; c = gseg[j].b;
			Pt I = geomMethods.intxPt(a, c, b, d);
			if (I.x != -1 && I.y != -1 && I.z != -1) {
				if (max(d.Dis(a), b.Dis(c)) < max(c.Dis(d), b.Dis(a))){
					qv.push_back(Quad(a, b, c, d));
				}				
			}
		}
		csSpine.clear(); vector<Pt>().swap(csSpine);
		blockvec[i].setQuadVec(qv);
	}
	std::cout << "Objective complete at genCrvSkeleton\nITERATION NUMBER= " << ITERATION << "\n\n\n" << endl;
}

void ofApp::nsGenCell() {
		
	//FOCUS ON ORTHO INTERIOR
	for (int i = 0; i < blockvec.size(); i++) {
		vector<Cell>cellvec;
		vector<vector<Cell>> rowcellvec;
		vector<Quad>qvec = blockvec[i].blockquadvec;
		for (int j = 0; j < qvec.size(); j++) {
			vector<Pt>bottom; vector<Pt>top;

			Pt P = qvec[j].p; Pt Q = qvec[j].q; float dpq = P.Dis(Q);
			Pt u((Q.x - P.x) / dpq, (Q.y - P.y) / dpq, (Q.z - P.z) / dpq);
			Pt p(P.x + u.x*SCALE_HULL, P.y + u.y*SCALE_HULL, P.z + u.z*SCALE_HULL);
			Pt q(Q.x - u.x*SCALE_HULL, Q.y - u.y*SCALE_HULL, Q.z - u.z*SCALE_HULL);

			Pt R = qvec[j].r; Pt S = qvec[j].s; float drs = R.Dis(S);
			Pt v((R.x - S.x) / drs, (R.y - S.y) / drs, (R.z - S.z) / drs);			
			Pt s(S.x + v.x*SCALE_HULL, S.y + v.y*SCALE_HULL, S.z + v.z*SCALE_HULL);
			Pt r(R.x - v.x*SCALE_HULL, R.y - v.y*SCALE_HULL, R.z - v.z*SCALE_HULL);
			
			//it is possible that a cell is too long on one side 
			//because the first segment goes unchecked - pq < rs: forms distorted trapezoidal form
			
			bottom.push_back(p); top.push_back(s);

			//handle the mid section
			int itr = CELL_LENGTH;
			float DIS = 0;
			while (DIS < p.Dis(q)) {
				Pt e;
				if (p.Dis(q) > r.Dis(s)) {
					e=Pt(p.x + u.x*itr, p.y + u.y*itr, p.z + itr * u.z);
					Pt f = geomMethods.proj(s, e, r);				
					//when projected, then it may not sit inside 
					if (abs(e.Dis(q) + e.Dis(p) - p.Dis(q)) < 1) {
						if (e.Dis(q) > CELL_LENGTH / 2) {
							if (abs(f.Dis(r) + f.Dis(s) - r.Dis(s)) < 1) {//if the point is outside segemtn but it can be above- handle intersection
								bottom.push_back(e); top.push_back(f);
							}
							else {//handle the intersection -> proj can be above then this dcase is wrong
								//intersection with the joiing edge from r-q NOT r-s
								//also handle simple case r-s
								Pt I = geomMethods.intxPt(e, f, r, q);
								Pt J = geomMethods.intxPt(e, f, r, s);
								if (I.x != -1 && I.y != -1 && I.z != -1) {
									if (I.Dis(r) > CELL_LENGTH / 2) {
										bottom.push_back(e); top.push_back(I);
									}
								}
								else if (J.x != -1 && J.y != -1 && J.z != -1) {
									if (J.Dis(r) > CELL_LENGTH / 2) {
										bottom.push_back(e); top.push_back(J);
									}
								}
							}
						}
					}
				}
				else {
					e=Pt(s.x + v.x*itr, s.y + v.y*itr, s.z + itr * v.z);
					Pt f = geomMethods.proj(p, e, q);
					//when projected, then it may not sit inside 
					if (abs(e.Dis(s) + e.Dis(r) - s.Dis(r)) < 1) {
						if (e.Dis(r) > CELL_LENGTH / 2) {
							if (abs(f.Dis(p) + f.Dis(q) - p.Dis(q)) < 1) {//if the point is outside segemtn but it can be above- handle intersection
								bottom.push_back(f); top.push_back(e);
							}
							else {//handle the intersection -> proj can be above then this dcase is wrong
								//intersection with the joiing edge from p-q NOT r-q
								//also handle simple case p-q
								Pt I = geomMethods.intxPt(e, f, r, q);
								Pt J = geomMethods.intxPt(e, f, p, q);
								if (I.x != -1 && I.y != -1 && I.z != -1) {
									if (I.Dis(q) > CELL_LENGTH / 2) {
										bottom.push_back(I); top.push_back(e);
									}									
								}
								else if (J.x != -1 && J.y != -1 && J.z != -1) {
									if (J.Dis(q) > CELL_LENGTH / 2){
										bottom.push_back(J); top.push_back(e);
									}
								}
							}
						}
					}
				}
				DIS = e.Dis(q);
				itr += CELL_LENGTH;
			}
			bottom.push_back(q); top.push_back(r); 

			//cout << "number of cells in quad (col)= " << num << endl;
			for (int k = 1; k < bottom.size(); k++) {
				Pt p0 = bottom[k - 1];
				Pt p1 = bottom[k];
				Pt p2 = top[k];
				Pt p3 = top[k - 1];
				float angArr[4];
				angArr[0] = geomMethods.computeAngle(p0, p1, p2);
				angArr[1] = geomMethods.computeAngle(p1, p2, p3);
				angArr[2] = geomMethods.computeAngle(p2, p3, p0);
				angArr[3] = geomMethods.computeAngle(p3, p0, p1);
				float minAngle = 1000;
				for (int ite = 0; ite < 4; ite ++) {
					if (angArr[ite] < minAngle) {
						minAngle = angArr[ite];
					}
				}
				float d0 = p0.Dis(p1); 
				float d1 = p0.Dis(p3); 
				float d2 = p1.Dis(p2);
				float d3 = p2.Dis(p3);
				float d4 = p0.Dis(p2);
				float d5 = p1.Dis(p3);
				int o0 = (d0 > 0.35*CELL_LENGTH && d0 < 3 * CELL_LENGTH);
				int o1 = (d1 > 0.35*CELL_LENGTH && d1 < 3 * CELL_LENGTH);
				int o2 = (d2 > 0.35*CELL_LENGTH && d2 < 3 * CELL_LENGTH);
				int o3 = (d3 > 0.35*CELL_LENGTH && d3 < 3 * CELL_LENGTH);
				int o4 = (d4 > 0.35*CELL_LENGTH && d4 < 3 * CELL_LENGTH);
				int o5 = (d5 > 0.35*CELL_LENGTH && d5 < 3 * CELL_LENGTH);
				if (o1 == 1 && o1 == 1 && o2 == 1 && o3 == 1 && o4==1 && o5==1 && minAngle>MIN_CELL_ANGLE) {
					
					Cell cell(p0, p1, p2, p3, k, j, bottom.size(), minAngle, colorVec[i]);
					cell.setMinAngle(minAngle);
					cellvec.push_back(cell);
				}			
			}		
			vector<Pt>().swap(bottom);
			vector<Pt>().swap(top);
		}		
		colorVec[i].display();
		blockvec[i].setCellVec(cellvec);
		vector<Cell>().swap(cellvec);
	}
	
	std::cout << "Objective complete at nsGENCELL" << endl;
	nsRules();
}

void ofApp::nsRules() {
	//interact cells and display updated  form
	for (int i = 0; i < blockvec.size(); i++) {
		vector<Cell>::iterator itr = blockvec[i].cellvec.begin();

		float req_angle = PI / 10;
		while(itr != blockvec[i].cellvec.end()){
			int J = (*itr).J; int num = (*itr).Num; float cellang = (*itr).CELL_ANGLE;
			//cout <<"check : "<< cellang << ", " << req_angle << endl;
			(*itr).CELL_FILL = 0; (*itr).CELL_HEIGHT = 0; 
			int ht, fill;
			if (J < 2 || J > num - 2) {
				int R = ofRandom(0, 10); 
				if (R > 0 && cellang>req_angle) {
					fill = 1;
				}
				else {
					fill = 0;
				}
			}
			else {
				int interior_hull = ofRandom(0, 10);
				if (interior_hull > 5 && cellang > req_angle) {
					fill = 1;
				}
				else {
					fill = 0;
				}				
			}
			if(fill==1){
				int S = ofRandom(0, 10);
				if (J<3 || J>num - 3) {
					if (S < 4) {
						ht = 3;
					}
					else if (S > 3 && S < 7) {
						ht = 5;
					}
					else {
						ht = 8;
					}
				}
				else {//interiro hull
					if (S < 4) {
						ht = 5;
					}
					else if (S > 3 && S < 5) {
						ht = 7;
					}
					else {
						ht = 18;
					}
				}				
				(*itr).CELL_HEIGHT=ht;
				(*itr).CELL_FILL=fill;
			}	
			itr++;
		}
	}
	std::cout << "Objective complete at nsINTERACTION" << endl;
	/* post process cells - merging */
}

float ofApp::nscomputeArea() {
	std::cout << "\n\n========================" << endl;
	std::cout << "\n\n\tAREA OCCUPIED" << endl;
	std::cout << "\n\n========================" << endl;
	vector<float> blockcellarea;
	float totalAreaOccupied=0;
	for (int i = 0; i < blockvec.size(); i++) {
		std::vector<Cell>cellvec = blockvec[i].cellvec;
		float sum_cells=0;
		for (int j = 0; j < cellvec.size(); j++) {
			Cell cell = cellvec[j];
			Pt p = cell.p; Pt q = cell.q; Pt r = cell.r; Pt s = cell.s;
			//2 triangles p,q,r + p,r,s
			float a1 = p.Dis(q), b1 = q.Dis(r), c1 = p.Dis(r);
			float a2 = p.Dis(r), b2 = r.Dis(s), c2 = p.Dis(s);
			float s1 = (a1 + b1 + c1) / 2; 
			float s2 = (a2 + b2 + c2) / 2;
			float ar1 = sqrt((s1 - a1)*(s1 - b1)*(s1 - c1));
			float ar2 = sqrt((s2 - a2)*(s2 - b2)*(s2 - c2));
			sum_cells += (ar1 + ar2)*cell.CELL_HEIGHT;
			totalAreaOccupied+= (ar1 + ar2)*cell.CELL_HEIGHT;
		}
		std::cout << "area occupied by block: <" << i << "> = " << sum_cells << endl;
		blockcellarea.push_back(sum_cells);
	}
	float site_area = BOARD_DIMENSION * BOARD_DIMENSION;
	float fsr = totalAreaOccupied / site_area;
	std::cout << "\nTotal area occupied by cells: " << totalAreaOccupied << endl;
	std::cout << "Total Area of region = " << site_area << endl;
	std::cout << "FSR achieved = " << fsr << endl;
	return fsr;
}


/* frameworks */

void ofApp::setup(){
	//light.enable();
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
			ofSetColor(255, 100, 0, 200);
			ofDrawLine(a.x, a.y, a.z, b.x, b.y, b.z);
		}
	}

	for (int i = 0; i < scaledHullPts.size(); i++) {
		vector<Pt> hull = scaledHullPts[i];
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
			ofSetLineWidth(3);
			ofSetColor(255, 0, 255, 100);
			ofDrawLine(a.x, a.y, a.z, b.x, b.y, b.z);
		}
	}


	//display CElls
	for (int i = 0; i < blockvec.size(); i++) {
		vector<Cell> cellvec = blockvec[i].cellvec;
		vector<Quad> quads = blockvec[i].blockquadvec;
		for (int j = 0; j < quads.size(); j++) {
			//quads[j].display();
		}
		for (int j = 0; j < cellvec.size(); j++) {
			//cellvec[j].display();
			cellvec[j].display2();//discrete cells unmerged MESH
			//cellvec[j].display3(); //plot diagonals
		}
	}
	
	ofSetLineWidth(1);
	cam.end();
	ofDrawBitmapStringHighlight("Press 'r' or 'R' to reconfigure organization\nPress 'i' or 'I' to reconfigure buildings", 10, 20);
	ofDrawBitmapStringHighlight("student: Nirvik Saha\nadviser: Dennis R Shelden, John R Haymaker", 10, ofGetWindowHeight() - 20);
}

void ofApp::keyPressed(int key) {
	if (key == 'r') { 		
		cout << "\n\n\nUSER Command: "<< ITERATION <<", for new board" << endl;
		nsInit(); 
	}

	if (key == 'i') {
		std::cout << "\n\n INTERACTion command" << endl;
		nsRules();
	}

	if (key == 's') {
		for (int i = 0; i < blockvec.size(); i++) {
			vector<Cell> cellvec = blockvec[i].cellvec;
			for (int j = 0; j < cellvec.size(); j++) {
				std::cout << " height" << cellvec[j].getHeight() << " Fill: " << cellvec[j].getFill() << endl;
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

