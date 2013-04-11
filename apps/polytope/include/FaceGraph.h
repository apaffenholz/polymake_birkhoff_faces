#include <polymake/Array.h>
#include <polymake/Bitset.h>
#include <polymake/Vector.h>
#include <polymake/Graph.h>

#ifndef POLYMAKE_POLYTOPE_FACE_GRAPH_H
#define POLYMAKE_POLYTOPE_FACE_GRAPH_H

namespace polymake {  namespace polytope {

    class FaceGraph {
      
    private :
      int _n;
      int _n_missing;
      int _n_degtwo;
      Array<int> _udeg, _ldeg;
      Array<int> _lower_min, _upper_min;
      Array<Bitset> _ul, _ll;
      Array<bool> _u_deg2, _l_deg2;
      Array<int> _min_u_attach;
      Array<int> _min_l_attach;
      
      bool _deg2complete;
      bool _deg3complete;

      FaceGraph() {};
      
    public : 

      friend std::ostream & operator<< ( std::ostream & o, const FaceGraph& eg );
      
    FaceGraph(int n, int d) : _n(n), _n_missing(d+2*n-1), _udeg(n), _ldeg(n), _lower_min(n), _upper_min(n), 
	_ul(n,Bitset(n)), _ll(n,Bitset(n)), _u_deg2(n,false), _l_deg2(n,false), _min_u_attach(n), _min_l_attach(n)
	{}

    FaceGraph(int n, int d, int k) : _n(n), _n_missing(d+2*n-1), _n_degtwo(k), _udeg(n), _ldeg(n), _lower_min(n), _upper_min(n), 
	_ul(n,Bitset(n)), _ll(n,Bitset(n)), _u_deg2(n,false), _l_deg2(n,false), _min_u_attach(n), _min_l_attach(n)
	{}


    FaceGraph(const FaceGraph& eg) : 
        _n(eg.size()), 
	_n_missing(eg.n_missing_edges()), 
	_n_degtwo(eg.nodes_of_deg2()), 
	_udeg(eg.upper_degrees()), 
	_ldeg(eg.lower_degrees()), 
	_lower_min(eg.u_minimal_lower_index()), 
	_upper_min(eg.l_minimal_upper_index()),
	_ul(eg.size(),entire(eg.upper())),  
	_ll(eg.size(),entire(eg.lower())),
	_u_deg2(eg.u_deg2()),
	_l_deg2(eg.l_deg2()),
	_deg2complete(eg.deg2complete()), 
	_deg3complete(eg.deg3complete()) {}
      
      FaceGraph & operator= ( FaceGraph& eg ) {
	_n = eg.size();
	_n_degtwo = eg.nodes_of_deg2();
	_udeg = eg.upper_degrees();
	_ldeg = eg.lower_degrees();
	_lower_min = eg.u_minimal_lower_index();
	_upper_min = eg.l_minimal_upper_index();
	_ul = eg.upper();
	_ll = eg.lower();
	_u_deg2 = eg.u_deg2();
	_l_deg2 = eg.l_deg2();
	_deg2complete = eg.deg2complete();
	_deg3complete = eg.deg3complete();
	_n_missing = eg.n_missing_edges();
	return *this;
      }
      
      const int size() const {return _n; }
      const int nodes_of_deg2() const {return _n_degtwo; }

      void set_n_degtwo (int k ) { _n_degtwo = k; }

      const Array<int> upper_degrees() const { return _udeg; }
      const Array<int> lower_degrees() const { return _ldeg; }
      const Array<Bitset> upper() const { return _ul; }
      const Array<Bitset> lower() const { return _ll; }
      const int upper_degree(int i) const { return _udeg[i]; }
      const int lower_degree(int i) const { return _ldeg[i]; }
      const Bitset upper(int i) const { return _ul[i]; }
      const Bitset lower(int i) const { return _ll[i]; }
      const bool deg2complete() const { return _deg2complete; }
      const bool deg3complete() const { return _deg3complete; }
      const int n_missing_edges() const { return _n_missing; }

      const int u_minimal_lower_index ( int i ) const { return _lower_min[i]; }
      const int l_minimal_upper_index ( int i ) const { return _upper_min[i]; }
      const Array<int> u_minimal_lower_index() const { return _lower_min; }
      const Array<int> l_minimal_upper_index() const { return _upper_min; }

      void set_u_min_lower_index ( int i, int j ) { _lower_min[i] = j; }
      void set_l_min_upper_index ( int i, int j ) { _upper_min[i] = j; }

      const Array<bool> u_deg2() const { return _u_deg2; }
      const Array<bool> l_deg2() const { return _l_deg2; }
      bool u_deg2(int i) const { return _u_deg2[i]; }
      bool l_deg2(int i) const { return _l_deg2[i]; }
      
      void set_u_deg2 ( int i, bool b ) { _u_deg2[i]=b; }
      void set_l_deg2 ( int i, bool b ) { _l_deg2[i]=b; }
      
      const int upper_first_unused() const { int k = 0; while ( k<_n && !_ul[k].empty() ) k++; return k; }
      const int lower_first_unused() const { int k = 0; while ( k<_n && !_ll[k].empty() ) k++; return k; }

      const int min_u_attach ( int i ) {
	if ( i == 0 || u_deg2(i) != u_deg2(i-1) )
	  return 0;
	
	int a = *(upper(i-1).begin());
	if ( a < _n-_n_degtwo ) 
	  return a;
	else 
	  return 0;
      }

      const int min_l_attach ( int i ) {
	if ( i == 0 || l_deg2(i) != l_deg2(i-1) )
	  return 0;

	int a = *(lower(i-1).begin());
	if ( a < _n-_n_degtwo ) 
	  return a;
	else 
	  return 0;
      }


      const bool is_connected() const;
      const bool is_elementary() const;
      const Vector<Vector<int> > matchings() const;
      const IncidenceMatrix<> as_matrix() const;

      const Graph<> as_graph() const;
      
      void add_edge ( int u, int l ); 
      const bool edge ( int u, int l ) const { return _ul[u].contains(l); }

      void add_ear ( int u, int l, int additional_nodes = 0 );
      
    };
        
  }}

#endif
