
// for Voronoi
#include <iostream>
#include <algorithm>
#include <set>

#include "Voronoi.h"  // THE HEADER

namespace VoronoiNamespace{

Voronoi::Voronoi(){
	edges = 0;
}

// ==== GetEdges - main algorithm

Edges * Voronoi::GetEdges(Vertices * v, double w, double h){
//Edges * Voronoi::GetEdges(){
	places = v;
	width  = w;
	height = h;
	root   = 0;

	if( !edges ) {
		 edges = new Edges();
	}else{ // clear old edges and vertices
		for(Vertices::iterator	i = points.begin(); i != points.end(); ++i) delete (*i);
		for(Edges::iterator		i = edges->begin(); i != edges->end(); ++i) delete (*i);
		points.clear();
		edges->clear();
	}

	for(Vertices::iterator i = places->begin(); i!=places->end(); ++i)	{  // put all places to event-que
		queue.push(new VEvent( *i, true));
	}

	VEvent * e;
	while( !queue.empty() )	{ // pop events from que by priority
		e = queue.top();
		queue.pop();
		ly = e->point->y;
		if(deleted.find(e) != deleted.end()) { delete(e); deleted.erase(e); continue;} // ?
		if(e->pe) InsertParabola(e->point);
		else RemoveParabola(e);
		delete(e);
	}

	FinishEdge(root);

	for(Edges::iterator i = edges->begin(); i != edges->end(); ++i)	{ // ?
		if( (*i)->neighbour) {
			(*i)->start = (*i)->neighbour->end;
			delete (*i)->neighbour;
		}
	}

	return edges;
}

// ==== InsertParabola

void	Voronoi::InsertParabola(Vec2d * p){
	if(!root){root = new VParabola(p); return;}

	if(root->isLeaf && root->site->y - p->y < 1){ // degenerovaný pripad - obì spodní místa ve stejné výce
		Vec2d * fp = root->site;
		root->isLeaf = false;
		root->SetLeft( new VParabola(fp) );
		root->SetRight(new VParabola(p)  );
		//Vec2d * s = new Vec2d( (p->x + fp->x)/2, height ); // zaèátek hrany uprostøed míst
		Vec2d * s = new Vec2d();   s->set( (p->x + fp->x)/2, height ); // zaèátek hrany uprostøed míst
		points.push_back(s);
		if(p->x > fp->x) root->edge = new VEdge(s, fp, p); // rozhodnu, který vlevo, který vpravo
		else root->edge = new VEdge(s, p, fp);
		edges->push_back(root->edge);
		return;
	}

	VParabola * par = GetParabolaByX(p->x);

	if(par->cEvent)	{
		deleted.insert(par->cEvent);
		par->cEvent = 0;
	}

	//Vec2d * start = new Vec2d(p->x, GetY(par->site, p->x));
    Vec2d * start = new Vec2d(); start->set( p->x, GetY(par->site, p->x) );
	points.push_back(start);

	VEdge * el = new VEdge(start, par->site, p);
	VEdge * er = new VEdge(start, p, par->site);

	el->neighbour = er;
	edges->push_back(el);

	// pøestavuju strom .. vkládám novou parabolu
	par->edge = er;
	par->isLeaf = false;

	VParabola * p0 = new VParabola(par->site);
	VParabola * p1 = new VParabola(p);
	VParabola * p2 = new VParabola(par->site);

	par->SetRight(p2);
	par->SetLeft(new VParabola());
	par->Left()->edge = el;

	par->Left()->SetLeft(p0);
	par->Left()->SetRight(p1);

	CheckCircle(p0);
	CheckCircle(p2);
}

// ==== RemoveParabola

void	Voronoi::RemoveParabola(VEvent * e){
	VParabola * p1 = e->arch;

	VParabola * xl = VParabola::GetLeftParent(p1);
	VParabola * xr = VParabola::GetRightParent(p1);

	VParabola * p0 = VParabola::GetLeftChild(xl);
	VParabola * p2 = VParabola::GetRightChild(xr);

	if(p0 == p2) std::cout << "chyba - pravá a levá parabola má stejné ohnisko!\n";

	if(p0->cEvent){ deleted.insert(p0->cEvent); p0->cEvent = 0; }
	if(p2->cEvent){ deleted.insert(p2->cEvent); p2->cEvent = 0; }

	//Vec2d * p = new Vec2d(e->point->x, GetY(p1->site, e->point->x));
	Vec2d * p = new Vec2d(); p->set( e->point->x, GetY(p1->site, e->point->x) );
	points.push_back(p);

	xl->edge->end = p;
	xr->edge->end = p;

	VParabola * higher;
	VParabola * par = p1;
	while(par != root)	{
		par = par->parent;
		if(par == xl) higher = xl;
		if(par == xr) higher = xr;
	}
	higher->edge = new VEdge(p, p0->site, p2->site);
	edges->push_back(higher->edge);

	VParabola * gparent = p1->parent->parent;
	if(p1->parent->Left() == p1)	{
		if(gparent->Left()  == p1->parent) gparent->SetLeft ( p1->parent->Right() );
		if(gparent->Right() == p1->parent) gparent->SetRight( p1->parent->Right() );
	}else{
		if(gparent->Left()  == p1->parent) gparent->SetLeft ( p1->parent->Left()  );
		if(gparent->Right() == p1->parent) gparent->SetRight( p1->parent->Left()  );
	}

	delete p1->parent;
	delete p1;

	CheckCircle(p0);
	CheckCircle(p2);
}

// ==== FinishEdge

void	Voronoi::FinishEdge(VParabola * n){
	if(n->isLeaf) {delete n; return;}
	double mx;
	if(n->edge->direction->x > 0.0){ mx = std::max(width, n->edge->start->x + 10 ); }
	else                           { mx = std::min(0.0,   n->edge->start->x - 10 ); }

	//Vec2d * end = new Vec2d(mx, mx * n->edge->f + n->edge->g);
	Vec2d * end = new Vec2d(); end->set( mx,  mx * n->edge->f + n->edge->g );
	n->edge->end = end;
	points.push_back(end);

	FinishEdge(n->Left() );
	FinishEdge(n->Right());
	delete n;
}

// ==== GetXOfEdge

double	Voronoi::GetXOfEdge( VParabola * par, double y){

	VParabola * left = VParabola::GetLeftChild(par);
	VParabola * right= VParabola::GetRightChild(par);

	Vec2d * p = left->site;
	Vec2d * r = right->site;

	double dp = 2.0 * (p->y - y);
	double a1 = 1.0 / dp;
	double b1 = -2.0 * p->x / dp;
	double c1 = y + dp / 4 + p->x * p->x / dp;

		   dp = 2.0 * (r->y - y);
	double a2 = 1.0 / dp;
	double b2 = -2.0 * r->x/dp;
	double c2 = ly + dp / 4 + r->x * r->x / dp;

	double a = a1 - a2;
	double b = b1 - b2;
	double c = c1 - c2;

	// solve quadratic equation
	double disc = b*b - 4 * a * c;
	double x1 = (-b + std::sqrt(disc)) / (2*a);
	double x2 = (-b - std::sqrt(disc)) / (2*a);

	if( p->y < r->y ) {
		return std::max(x1, x2);
	}else{
		return std::min(x1, x2);
	}
}

// ==== GetParabolaByX

VParabola * Voronoi::GetParabolaByX(double xx){
	VParabola * par = root;
	double x = 0.0;

	while(!par->isLeaf) {
		// projdu stromem dokud nenarazím na vhodný list
		x = GetXOfEdge(par, ly);
		if(x>xx){
			par = par->Left();
		}else{
			par = par->Right();
		}
	}
	return par;
}

// ==== GetY

double	Voronoi::GetY(Vec2d * p, double x) { // ohnisko, x-souøadnice
	double dp = 2 * (p->y - ly);
	double a1 = 1 / dp;
	double b1 = -2 * p->x / dp;
	double c1 = ly + dp / 4 + p->x * p->x / dp;
	return(a1*x*x + b1*x + c1);
}

// ==== CheckCircle

void	Voronoi::CheckCircle(VParabola * b){

	VParabola * lp = VParabola::GetLeftParent (b);
	VParabola * rp = VParabola::GetRightParent(b);

	VParabola * a  = VParabola::GetLeftChild (lp);
	VParabola * c  = VParabola::GetRightChild(rp);

	if(!a || !c || a->site == c->site) return;

	Vec2d * s = 0;
	s = GetEdgeIntersection(lp->edge, rp->edge);
	if(s == 0) return;

	double dx = a->site->x - s->x;
	double dy = a->site->y - s->y;

	double d = std::sqrt( (dx * dx) + (dy * dy) );

	if(s->y - d >= ly) { return;}

    Vec2d * newv = new Vec2d(); newv->set( s->x,  s->y - d );
    VEvent * e = new VEvent( newv, false);
	//VEvent * e = new VEvent( new Vec2d(s->x, s->y - d), false);
	points.push_back(e->point);
	b->cEvent = e;
	e->arch = b;
	queue.push(e);

}

// ==== GetEdgeIntersection

Vec2d * Voronoi::GetEdgeIntersection(VEdge * a, VEdge * b){

	double x = (b->g-a->g) / (a->f - b->f);
	double y = a->f * x + a->g;

	if((x - a->start->x)/a->direction->x < 0) return 0;
	if((y - a->start->y)/a->direction->y < 0) return 0;

	if((x - b->start->x)/b->direction->x < 0) return 0;
	if((y - b->start->y)/b->direction->y < 0) return 0;

	//Vec2d * p = new Vec2d(x, y);
	Vec2d * p = new Vec2d(); p->set(x, y);
	points.push_back(p);
	return p;

}



}

