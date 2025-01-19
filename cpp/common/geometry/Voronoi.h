#ifndef Voronoi_h
#define Voronoi_h

#include <list>
#include <queue>
#include <set>

#include "Vec2.h"

// from here: http://blog.ivank.net/fortunes-algorithm-and-implementation.html

namespace VoronoiNamespace{

	class VEvent;
	class VEdge;

	//===============
    //   VParabola
	//===============

	class VParabola{
		public:
		bool		isLeaf;   //  flag whether the node is Leaf or internal node
		Vec2d      * site;     //  pointer to the focus point of parabola (when it is parabola)
		VEdge     * edge;     //  pointer to the edge (when it is an edge)
		VEvent    * cEvent;   //  pointer to the event, when the arch disappears (circle event)
		VParabola * parent;   //  pointer to the parent node in tree

		VParabola * _left;
		VParabola * _right;

		void		SetLeft (VParabola * p) {_left  = p; p->parent = this;}
		void		SetRight(VParabola * p) {_right = p; p->parent = this;}

		VParabola *	Left () { return _left;  }
		VParabola * Right() { return _right; }

/*
		static VParabola * GetLeft			(VParabola * p);   //  returns the closest left leave of the tree
		static VParabola * GetRight			(VParabola * p);   //  returns the closest right leafe of the tree
		static VParabola * GetLeftParent	(VParabola * p);   //  returns the closest parent which is on the left
		static VParabola * GetRightParent	(VParabola * p);   //  returns the closest parent which is on the right
		static VParabola * GetLeftChild		(VParabola * p);   //  returns the closest leave which is on the left of current node
		static VParabola * GetRightChild	(VParabola * p);   //  returns the closest leave which is on the right of current node

*/

		VParabola(){
			site	= 0;
			isLeaf	= false;
			cEvent	= 0;
			edge	= 0;
			parent	= 0;
		}

		VParabola(Vec2d * s){
			site	= s;
			isLeaf	= true;
			cEvent	= 0;
			edge	= 0;
			parent	= 0;
		}

		// ==== Tree operations

		static VParabola * GetLeft	(VParabola * p){
			return GetLeftChild(GetLeftParent(p));
		}

		static VParabola * GetRight	(VParabola * p){
			return GetRightChild(GetRightParent(p));
		}

		static VParabola * GetLeftParent	(VParabola * p){
			VParabola * par		= p->parent;
			VParabola * pLast	= p;
			while(par->Left() == pLast) {
				if(!par->parent) return 0;
				pLast = par;
				par = par->parent;
			}
			return par;
		}

		static VParabola * GetRightParent (VParabola * p){
			VParabola * par		= p->parent;
			VParabola * pLast	= p;
			while(par->Right() == pLast) {
				if(!par->parent) return 0;
				pLast = par; par = par->parent;
			}
			return par;
		}

		static VParabola * GetLeftChild (VParabola * p){
			if(!p) return 0;
			VParabola * par = p->Left();
			while(!par->isLeaf) par = par->Right();
			return par;
		}

		static VParabola * GetRightChild (VParabola * p){
			if(!p) return 0;
			VParabola * par = p->Right();
			while(!par->isLeaf) par = par->Left();
			return par;
		}


	};

	//===============
    //   VEdge
	//===============

	class VEdge{
		public:

		Vec2d *	start;      //  pointer to start point
		Vec2d *	end;        //  pointer to end point
		Vec2d *	direction;  //  irectional vector, from "start", points to "end", normal of |left, right|
		Vec2d *	left;       //  pointer to Voronoi place on the left side of edge
		Vec2d *	right;      //  pointer to Voronoi place on the right side of edge

		double	f;          //  directional coeffitients satisfying equation y = f*x + g (edge lies on this line)
		double	g;          //

		VEdge * neighbour;      // some edges consist of two parts, so we add the pointer to another part to connect them at the end of an algorithm

		VEdge(Vec2d * s, Vec2d * a, Vec2d * b)	{
			start		= s;  //
			left		= a;
			right		= b;
			neighbour	= 0;
			end			= 0;

			f = (b->x - a->x) / (a->y - b->y) ;
			g = s->y - f * s->x ;
			//direction = new Vec2d(b->y - a->y, -(b->x - a->x));
			direction = new Vec2d();  direction->set( b->y - a->y, -(b->x - a->x) );
		}

		~VEdge(){
			delete direction ;
		}

	};

	//===============
    //   VEvent
	//===============

	class VEvent{
		public:
		Vec2d *	point;      // the point at which current event occurs (top circle point for circle event, focus point for place event)
		bool		pe;     // whether it is a place event or not
		double		y;      // y coordinate of "point", events are sorted by this "y"
		VParabola * arch;   // if "pe", it is an arch above which the event occurs

		VEvent(Vec2d * p, bool pev){
			point	= p;
			pe		= pev;  // point, at which the event occurs
			y		= p->y;
			arch	= 0;
		}

		struct CompareEvent : public std::binary_function<VEvent*, VEvent*, bool>{
			bool operator()(const VEvent* l, const VEvent* r) const { return (l->y < r->y); }
		};

	};

	typedef std::list<Vec2d *>		Vertices;
	typedef std::list<VEdge *>		Edges;

	//===============
    //   Voronoi
	//===============

	class Voronoi{ 
        public:
		Voronoi();
		Edges *			GetEdges( Vertices * v, double w, double h ); // main algorithm
		//Edges *		GetEdges( );     // main algorithm

		Vertices  *	places;          //  container of places with which we work
		Edges     *	edges;           //  container of edges which will be the result
		double		width, height;   //  width of the diagram
		VParabola *	root;            //  height of the diagram
		double		ly;              //  the root of the tree, that represents a beachline sequence

        private:

		std::set<VEvent *> deleted;                                                         //  set  of deleted (false) Events (since we can not delete from PriorityQueue
		std::list<Vec2d *> points;                                                          //  list of all new points that were created during the algorithm
		std::priority_queue<VEvent *, std::vector<VEvent *>, VEvent::CompareEvent> queue;   //  priority queue with events to process

		void        InsertParabola(Vec2d * p);                  //   processing the place event
		void        RemoveParabola(VEvent * e);                 //   processing the circle event
		void        FinishEdge(VParabola * n);                  //   recursively finishes all infinite edges in the tree
		double      GetXOfEdge(VParabola * par, double y);      //   returns the current x position of an intersection point of left and right parabolas
		VParabola * GetParabolaByX(double xx);                  //   returns the current x position of an intersection point of left and right parabolas
		double      GetY(Vec2d * p, double x);                  //
		void        CheckCircle(VParabola * b);                 //   checks the circle event (disappearing) of this parabola
		Vec2d     * GetEdgeIntersection(VEdge * a, VEdge * b);  //
	};


}


#endif
