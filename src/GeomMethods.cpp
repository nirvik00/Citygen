#include "GeomMethods.h"



GeomMethods::GeomMethods()
{
}


GeomMethods::~GeomMethods()
{
}

vector<Pt> GeomMethods::genSitePts(int t, int DI) {
	vector<Pt> ptVec;
	for (int i = 0; i < t; i++) {
		float x = ofRandomf() * DI;
		float y = 0;
		float z = ofRandomf() * DI;
		ptVec.push_back(Pt(x, 0, z));
	}
	return ptVec;
}

vector<Pt> GeomMethods::genGridPts(int le , int wi, int DI) {
	vector<Pt> ptVec;
	for (int i= -DI; i < DI; i+=le) {
		for (int j = -DI; j < DI; j += wi) {
			ptVec.push_back(Pt(i, 0, j));
		}
	}
	return ptVec;
}

vector<Pt> GeomMethods::setSiteToGrid(vector<Pt>sitePts, vector<Pt>gridPts) {
	for (int i = 0; i < gridPts.size(); i++) {
		Pt g = gridPts[i];
		float minD = 100000000; 
		float idx = -1;
		for (int j = 0; j < sitePts.size(); j++) {
			Pt s2 = sitePts[j];
			float d2 = s2.Dis(g);
			if (d2 < minD) {
				idx = j;
				minD = d2;
			}
		}
		if (idx > -1) {
			gridPts[i].setId(idx);
		}		
	}
	return gridPts;
}

vector<vector<Pt>> GeomMethods::gridPtsInSite(vector<Pt> sitePtVec, vector<Pt>gridPtVec) {
	vector<vector<Pt>> gridPtsInHull;
	for (int i = 0; i < sitePtVec.size(); i++) {
		vector<Pt> temp;
		for (int j = 0; j < gridPtVec.size(); j++) {
			int t = gridPtVec[j].getId();
			if (t == i) {
				temp.push_back(gridPtVec[j]);
			}
		}
		gridPtsInHull.push_back(temp);
		temp.clear(); vector<Pt>().swap(temp);
	}
	return gridPtsInHull;
}


//convex hull
vector<vector<Pt>> GeomMethods::initConvexHull(vector<Pt> sitePtVec, vector<Pt>gridPtVec) {
	vector<vector<Pt>> hullPts;
	for (int i = 0; i < sitePtVec.size(); i++) {
		vector<Pt> temp;
		for (int j = 0; j < gridPtVec.size(); j++) {
			int t = gridPtVec[j].getId();
			if (t == i) {
				temp.push_back(gridPtVec[j]);
			}
		}
		vector<Pt>tempHullPts = genConvexHull(temp);
		if (tempHullPts.size() > 3) {
			hullPts.push_back(tempHullPts);
		}
		temp.clear(); vector<Pt>().swap(temp);
		tempHullPts.clear(); vector<Pt>().swap(tempHullPts);
	}
	return hullPts;
}

int GeomMethods::orientationPts(Pt p, Pt q, Pt r) {
	float val = (q.z - p.z)*(r.x - q.x) - (r.z - q.z)*(q.x - p.x);
	if (val == 0) return 0;
	return (val > 0) ? 1 : 2;
}

vector<Pt> GeomMethods::genConvexHull(vector<Pt> pts) {
	vector<Pt> hull;
	int n = pts.size();
	if (n < 3) { 
		return hull; 
	}
	int r = 0;
	/*
	*  find the least x coordinate - this is the first point 
	*/
	for (int i = 1; i < n; i++) {
		if (pts[i].x < pts[r].x) { 
			r = i; 
		}
	}
	/* 
	*  initialize the first point
	*  until the iteration reaches the last point
	*  take 3 points such that their orientation is always ccw
	*/
	int p = r, q;
	do {
		hull.push_back(pts[p]);
		q = (p + 1) % n;
		for (int i = 0; i < n; i++) {
			if (orientationPts(pts[p], pts[i], pts[q]) == 2) {
				q = i;
			}
		}
		p = q;
	} while (p != r); 
	hull.push_back(hull[0]);
	return hull;
}


// utility
int GeomMethods::ptInPoly(vector<Pt> pts, Pt P) {
	Pt extreme(P.x + 100000, P.y, P.z);
	int count = 0;
	for (int i = 0; i < pts.size(); i++) {
		Pt a, b;
		if (i == 0) {
			a = pts[pts.size() - 1];
			b = pts[0];
		}
		else {
			a = pts[i - 1];
			b = pts[i];
		}
		int T = doesintx(a, b, P, extreme);
		if (T == 1) {
			count++;
		}
	}
	if (count % 2 == 1) {
		return 1;
	}
	else {
		return 0;
	}
}

vector<Pt> GeomMethods::simplify(vector<Pt> pts) {
	vector<Pt> new_pts;
	for (int i = 0; i < pts.size(); i++) {
		Pt p(pts[i].x, pts[i].y, pts[i].z);
		new_pts.push_back(p);
	}
	//for angles
	for (int i = 0; i < pts.size(); i++) {
		Pt a = pts[i];
		for (int j = i; j < pts.size(); j++) {
			Pt b = pts[j];
			for (int k = i; k < j; k++) {
				Pt c = pts[k];
				double dab = a.Dis(b);
				double dbc = b.Dis(c);
				double dca = c.Dis(a);
				if (dab > 1 && dbc > 1 && dca > 1) {
					double den_ab = b.x - a.x;
					if (abs(b.x - a.x) < 0.0001) {
						den_ab = 0.00001;
					}
					double den_bc = b.x - c.x;
					if (abs(b.x - c.x) < 0.0001) {
						den_bc = 0.00001;
					}
					double slope_ab = atan((b.z - a.z) / den_ab)*180/PI;
					double slope_bc = atan((b.z - c.z) / den_bc)*180/PI;

					bool o1 = c.x > min(a.x, b.x);
					bool o2 = c.x < max(a.x, b.x);
					bool o3 = c.y > min(a.y, b.y);
					bool o4 = c.y < max(a.y, b.y);
					int dir = 0;
					if (o1 == true && o2 == true && o3 == true && o4 == true) {
						dir = 1;
					}

					if ((abs(slope_ab - slope_bc) > 179 || abs(slope_ab - slope_bc) < 1) && dir==1) {
						int req=-1;
						for (int ite = 0; ite < new_pts.size(); ite++) {
							Pt d = new_pts[ite];
							if (c.x == d.x && c.y == d.y && c.z == d.z) {
								req = ite;
								break;
							}
							ite++;
						}
						if (req > -1) {
							new_pts.erase(new_pts.begin() + req);
						}
					}
				}
			}
		}
	}
	//for distance
	/*
	vector<Pt> rev;
	for (int i = 0; i < new_pts.size(); i++) {
		Pt p(new_pts[i].x, new_pts[i].y, new_pts[i].z);
		rev.push_back(p);
	}
	for (int i = 0; i < new_pts.size(); i++) {
		Pt a = pts[i];
		for (int j = i; j < new_pts.size(); j++) {
			Pt b = pts[j];
			double dab = a.Dis(b);
			//cout << "distances " << dab << endl;
			if (dab < 10 && dab>1) {
				int req = -1;
				for (int ite = 0; ite < rev.size(); ite++) {
					Pt d = rev[ite];
					if (b.x == d.x && b.y == d.y && b.z == d.z) {
						req = ite;
						break;
					}
				}
				if (req > -1) {
					rev.erase(rev.begin() + req);
				}				
			}
		}
	}
	*/
	return new_pts;
}

Pt GeomMethods::proj(Pt a, Pt b, Pt c) {
	//proj b on ac
	Pt u(a.x - c.x, a.y - c.y, a.z - c.z); Pt v(b.x - c.x, b.y - c.y, b.z - c.z);
	float e = (u.x*v.x + u.y*v.y + u.z*v.z) / (c.Dis(a)*c.Dis(a));
	Pt r(c.x + (e*u.x), c.y + (e*u.y), c.z + (e*u.z));
	return r;
}

int GeomMethods::doesintx(Pt p, Pt q, Pt r, Pt s) {
	float a1 = q.z - p.z; float b1 = p.x - q.x;	float c1 = a1 * q.x + b1 * q.z;
	float a2 = s.z - r.z; float b2 = r.x - s.x; float c2 = a2 * s.x + b2 * s.z;
	float det = a1 * b2 - a2 * b1; if (det == 0) { return 0; }
	float mx = (c1*b2 - c2 * b1) / det; float my = (c2*a1 - c1 * a2) / det;
	Pt I(mx, 0, my);
	if ((abs(I.Dis(r) + I.Dis(s) - r.Dis(s)) < 1) && (abs(I.Dis(p) + I.Dis(q) - p.Dis(q)) < 1)) { return 1; }
	else { return 0; }
}

Pt GeomMethods::intxPt(Pt p, Pt q, Pt r, Pt s) {
	float a1 = q.z - p.z; float b1 = p.x - q.x;	float c1 = a1 * q.x + b1 * q.z;
	float a2 = s.z - r.z; float b2 = r.x - s.x; float c2 = a2 * s.x + b2 * s.z;
	float det = a1 * b2 - a2 * b1; if (det == 0) { return Pt(-1,-1,-1); }
	float mx = (c1*b2 - c2 * b1) / det; float my = (c2*a1 - c1 * a2) / det;
	Pt I(mx, 0, my);
	if ((abs(I.Dis(r) + I.Dis(s) - r.Dis(s)) < 1) && (abs(I.Dis(p) + I.Dis(q) - p.Dis(q)) < 1)) { return I; }
	else { return Pt(-1,-1,-1); }
}


//occupation
vector<vector<Pt>> GeomMethods::spine(vector<Pt> pts, float DE, float SP) {
	float d0=-1;
	Pt A, B;
	for (int i = 0; i < pts.size(); i++) {
		for (int j = 0; j < pts.size(); j++) {
			Pt a = pts[i]; Pt b = pts[j];
			if (a.Dis(b) > d0) {
				A = a; B = b; d0 = a.Dis(b);
			}
		}
	}
	vector<vector<Pt>> PTS; //R.push_back(A); R.push_back(B); //spine
	
	float d = A.Dis(B);
	for (float i = 0; i < d; i += DE) {
		vector<Pt>temp;
		float sc = 1000000; //scale factor for intersection

		Pt g((B.x - A.x) / d, (B.y - A.y) / d, (B.z - A.z) / d); //unit vector
		Pt v(-g.z, g.y, g.x); //normal 1 
		Pt w(g.z, g.y, -g.x); //normal 2

		Pt u((B.x - A.x)*i / d, (B.y - A.y)*i / d, (B.z - A.z)*i / d); //increment spine points - bays depths


		Pt p(A.x + u.x, A.y + u.y, A.z + u.z);
		Pt q(p.x + v.x*sc, p.y + v.y*sc, p.z + v.z*sc); //along normal 1
		Pt r(p.x + w.x*sc, p.y + w.y*sc, p.z + w.z*sc); //along normal 2

		//Pt pL(p.x + v.x*SP, p.y + v.y*SP, p.z + v.z*SP); //offset 1
		//Pt pR(p.x + w.x*SP, p.y + w.y*SP, p.z + w.z*SP); //offset 2
		//temp.push_back(p); //spine point always index 0
		int prevOr;
		for (int j = 0; j < pts.size(); j++) {
			Pt a, b;
			if (j == 0) {
				a = pts[pts.size() - 1];
				b = pts[0];
			}
			else {
				a = pts[j - 1];
				b = pts[j];
			}
			//orientation
			int o1 = orientationPts(A, B, v); int o2 = orientationPts(A, B, w);
			Pt I, J;
			if (j == 0) {
				prevOr = o1;
				I = intxPt(a, b, p, q);
				J = intxPt(a, b, p, r);
			}
			else {
				if (prevOr == o1) {
					I = intxPt(a, b, p, q);
					J = intxPt(a, b, p, r);
				}
				else {
					I = intxPt(a, b, p, r);
					J = intxPt(a, b, p, q);
				}
			}
			if (I.x != -1 && I.y != -1 && I.z != -1) {
				temp.push_back(I);
			}
			if (J.x != -1 && J.y != -1 && J.z != -1) {
				temp.push_back(J);
			}
		}
		PTS.push_back(temp);
		temp.clear(); vector<Pt>().swap(temp);
	}
	return PTS;
}






//spline
vector<Pt> GeomMethods::cubicSubdiv(vector<Pt> pts) {
	vector<Pt> temp;
	for (int i = 0; i < pts.size(); i++) {
		Pt p, q;
		if (i < pts.size() - 1) {
			p = pts[i];
			q = pts[i + 1];
		}
		else {
			p = pts[i];
			q = pts[0];
		}
		temp.push_back(p);
		Pt r = cubicSplit(p, q);
		temp.push_back(r);
	}
	vector<Pt> spline = cubicAvg(temp);
	temp.clear(); vector<Pt>().swap(temp);
	return spline;
}

Pt GeomMethods::cubicSplit(Pt p, Pt q) {
	return Pt((p.x + q.x) / 2, (p.y + q.y) / 2, (p.z + q.z) / 2);
}

vector<Pt> GeomMethods::cubicAvg(vector<Pt> pts) {
	vector<Pt> temp;
	for (int i = 0; i < pts.size(); i++) {
		Pt p, q, r;
		if (i == 0) {
			p = pts[pts.size() - 1];
			q = pts[0];
			r = pts[1];
		}
		else if(i>0 && i<pts.size()-1) {
			p = pts[i - 1];
			q = pts[i];
			r = pts[i + 1];
		}
		else {
			p = pts[i];
			q = pts[0];
			r = pts[1];
		}
		float d = p.Dis(r);
		float sc = 1;
		Pt u(r.x - p.x, r.y - p.y, r.z - p.z);
		Pt v(-u.z*sc/d, u.y*sc/d, u.x*sc/d);
		Pt s;
		if (i < pts.size() - 1) {
			s = cubicSplit(q, r);
		}
		else {
			s = cubicSplit(p, q);
		}
		Pt w(s.x + v.x, s.y + v.y, s.z + v.x);
		temp.push_back(w);
	}
	return temp;
}