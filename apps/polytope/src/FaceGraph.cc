#include <polymake/Array.h>
#include <polymake/IncidenceMatrix.h>
#include "polymake/polytope/FaceGraph.h"
#include "polymake/Set.h"
#include "polymake/PowerSet.h"


namespace polymake {  namespace polytope {
    
    Vector<Vector<int> > matchings_rec ( const Array<Bitset>& ul, const Array<Bitset>& ll ) {
      
      Vector<Vector<int> > mlist(0);  // initialize empty list for the matchings generated in this step
      if ( ul.size() == 1 ) {         // not much to do anymore: graph has two nodes that are either connected or not.
	if ( ul[0].size() != 0 ) {
	  Vector<int> v(1); 
	  v[0]=0;
	  mlist |= v;
	}
      } else {                       // now we have a graph with at least two nodes on each side
	int l = (*ul[0].begin());   // pick an edge. We split into matching containing that edge and all others
	if ( ul[0].size() > 1 && ll[l].size() > 1 ) {   // otherwise there is no matching *not* containing the edge (0,l)
	  Array<Bitset> copy_ul(ul);                    // copy the graph and remove the edge.
	  Array<Bitset> copy_ll(ll);
	  copy_ul[0] -= l;
	  copy_ll[l] -= 0;
	  mlist |= matchings_rec ( copy_ul, copy_ll );  // get the matchings in the new graph
	}
	// now we deal with the matchings containing edge (0,l)
	// this is more complicated: we have to remove the two nodes and later reinsert them
	// this requires renumbering the nodes in the graph with 0,l deleted, recursing, and then again renumbering the nodes
	int d = ul.size()-1;
	Array<Bitset> copy_ul(d);  
	Array<Bitset> copy_ll(d);
	bool other_node_of_deg0 = false;
	for ( int i = 0; i < d && !other_node_of_deg0; ++i ) {
	  Bitset un(d), ln(d);
	  for ( Entire<Bitset>::const_iterator uliIt = entire( ul[i+1]); !uliIt.at_end(); ++uliIt )
	    if ( (*uliIt) < l ) 
	      un+= *uliIt;
	    else 
	      if ( (*uliIt) > l )
		un+= (*uliIt)-1;
	  for ( Entire<Bitset>::const_iterator lliIt = entire( ll[ i<l ? i : i+1]); !lliIt.at_end(); ++lliIt )
	    if ( *lliIt != 0 ) 
	      ln+= (*lliIt)-1;  
	  if ( un.empty() || ln.empty() )
	    other_node_of_deg0 = true;
	  copy_ul[i] = un;
	  copy_ll[i] = ln;
	}
	if ( !other_node_of_deg0 ) {
	  Vector<Vector<int> > mlist_temp = matchings_rec ( copy_ul, copy_ll );
	  for ( Entire<Vector<Vector<int> > >::iterator vit = entire(mlist_temp); !vit.at_end(); ++vit ) {
	    Vector<int> v(d+1); 
	    v[0]=l;
	    for ( int j = 0; j < d; ++j )
	      if ( (*vit)[j] < l ) 
		v[j+1] = (*vit)[j];
	      else 
		v[j+1] = (*vit)[j]+1;
	    *vit = v;
	  }
	  mlist |= mlist_temp;
	}
      }
	
      return mlist;
    }


    const Vector<Vector<int> > FaceGraph::matchings () const {
      return matchings_rec( _ul, _ll );
    }

    Vector<Vector<int> > perfect_matchings (Array<Set<int> > ul, Array<Set<int> > ll) {
      return matchings_rec( ul, ll );
    }
    
    // Check whether the matrix M is fully indecomposable
    const bool FaceGraph::is_connected() const {
      
      Bitset U(sequence(1,_n-1)), L(_ll[0]);
      bool took_union = false;
      int i = 1;
      do {
	if ( U.contains(i) && !(L*_ul[i]).empty() ) {
	  took_union = true;
	  L += _ll[i];
	  U -= i;
	  i = *(U.begin());
	  continue;
	} else {
	  ++i;
	}
      } while ( took_union && i < _n );

      return U.empty();
    }


    const IncidenceMatrix<> FaceGraph::as_matrix() const {

      IncidenceMatrix<> G (_n,_n);
      for ( int i = 0; i < _n; ++i )
	for ( Entire<Bitset>::const_iterator bit = entire(_ul[i]); !bit.at_end(); ++bit )
	  G(i,*bit) = 1;

      return G;
    }

    const Graph<> FaceGraph::as_graph() const {

      Graph<> G(2*_n);
      for ( int i = 0; i < _n; i++ ) 
	for ( Entire<Bitset>::const_iterator adj = entire(_ul[i]); !adj.at_end(); ++adj )
	  G.edge(i,*adj+_n);

      return G;
    }
  
  // Test, whether a matrix M defines an elementary graph
  // test criterion is G=(A,B), N(S) > |S| for all subsets S of A
  // see Carvalho, ... Thm 1.2
    const bool FaceGraph::is_elementary() const {
      
      const Set<int> S();
      
      for( int k = 1; k < _n; k++ )  {
	for (Subsets_of_k<const sequence&>::iterator pit=entire(all_subsets_of_k(sequence(0,_n),k)); !pit.at_end(); ++pit) {
	  Set<int> NS;
	  Set<int> upper_level_choice(*pit);
	  
	  for ( Entire<Set<int> >::const_iterator sit = entire(upper_level_choice); !sit.at_end(); ++sit )
	    NS += _ul[*sit];
	  
	  if ( NS.size() <= k )
	    return false;
	}   
      } 
      
      return true;
    }


      
    void FaceGraph::add_ear ( int u, int l, int additional_nodes ) {

      if ( _udeg[u] == 0 || _ldeg[l] == 0 )
	throw std::runtime_error("attaching nodes for ear have degree 0");

      if ( additional_nodes > 0 ) {
	int up = upper_first_unused();
	int lo = lower_first_unused();

	if ( _n-up < additional_nodes || _n-lo < additional_nodes )
	  throw std::runtime_error("not enough space for ear");

	add_edge( u, lo );
	add_edge( up, l );
	for ( int j = 0; j < additional_nodes-1; ++j ) {
	  add_edge( up++, lo+1);
	  add_edge( up, lo++);
	}
	add_edge(up, lo);
      } else {
	add_edge ( u,l );
      }
    }

    void FaceGraph::add_edge ( int u, int l ) { 
      if ( _ul[u].contains(l) || _ll[l].contains(u) ) {
	cout << "upper: " << u << " lower: " << l << endl << "containment: " << _ul[u].contains(l) << " " << _ll[l].contains(u) 
	     << endl << "in: " << _ul[u] << " and " << _ll[l] << endl;;
	throw std::runtime_error("attempt to attach existing edge");
      }
      --_n_missing;
      ++_udeg[u]; ++_ldeg[l]; 
      _ul[u] += l; _ll[l] += u; 
    }


    std::ostream & operator << ( std::ostream & o, const FaceGraph& eg ) {
      wrap(o) << "************************************************************************" << endl;
      wrap(o) << "upper level: " << endl << eg.upper() << endl << "-------" << endl;
      wrap(o) << "lower level: " << endl << eg.lower() << endl << "-------" << endl;
      wrap(o) << "upper degrees: " << eg.upper_degrees() << endl;
      wrap(o) << "lower degrees: " << eg.lower_degrees() << endl << "------" << endl;
      wrap(o) << "first unused nodes: " << eg.upper_first_unused() << "   " << eg.lower_first_unused() << endl;
      wrap(o) << "************************************************************************" << endl;
      return o;
    }

    UserFunction4perl(" ", &perfect_matchings, "perfect_matchings($,$)");

  }
}
