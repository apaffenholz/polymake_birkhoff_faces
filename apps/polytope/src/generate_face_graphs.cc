#include  <algorithm>
#include <polymake/client.h>
#include <polymake/Array.h>
#include <polymake/IncidenceMatrix.h>
#include "polymake/polytope/FaceGraph.h"
#include <polymake/FacetList.h>
#include <../../ext/nauty/apps/graph/include/graph_compare.h>  // currently assumes that polymake runs with nauty

//# define BIRKHOFF_DEBUG

namespace polymake { namespace polytope {

    IncidenceMatrix<> verticesInFacets ( const Array<Bitset>& ul, const Vector<Vector<int> >& matchings ) {

      FacetList FL(matchings.dim());
      for ( int i = 0; i < ul.size(); ++i ) {
	for ( Entire<Bitset>::const_iterator bit = entire(ul[i]); !bit.at_end(); ++bit ) {
	  Set<int> F;
	  for ( int j = 0; j < matchings.dim(); ++j )
	    if ( matchings[j][i] != *bit )
	      F += j;
	  FL.insertMax(F);
	}
      }

      IncidenceMatrix<> I(FL.size(), matchings.dim(), entire(FL));
      return I;
    }

    void AddNewGraph ( FaceGraph & fg, std::vector<std::pair<Graph<>,IncidenceMatrix<> > > & result ) {

      if ( fg.is_elementary() && fg.is_connected() ) {
	IncidenceMatrix<> I = verticesInFacets(fg.upper(), fg.matchings() );
	bool newPoly = true;
	for ( Entire<std::vector<std::pair<Graph<>, IncidenceMatrix<> > > >::const_iterator aIt = entire(result); !aIt.at_end(); ++aIt )
	  //	  if ( CallPolymakeFunction("isomorphic", aIt->second, I ) ) {  # FIXME no direct isomorphism check for incidence matrices implemented
	  if ( graph::isomorphic( aIt->second, I ) ) {
	    newPoly = false;
	    break;
	  }
	if ( newPoly ) {
	  std::pair<Graph<>, IncidenceMatrix<> > p;
	  p.first = fg.as_graph();
	  p.second = I;
	  result.push_back(p);
#ifdef BIRKHOFF_DEBUG
	  cout << "final graph: " << fg << endl;
#endif
	}
      }
    }

    void AddEar (FaceGraph& eg, int remainingFreeNodes, int remainingEdges, int max_length, std::vector<std::pair<Graph<>,IncidenceMatrix<> > >& result ) {
      if ( remainingFreeNodes > 0 ) 
	for ( int i = 1; i <= std::min(max_length,remainingFreeNodes); ++i ) 
	  for ( int u = 0; u < eg.upper_first_unused(); ++u ) 
	    for ( int l = 0; l < eg.lower_first_unused(); ++l ) {
	      FaceGraph eg_copy(eg);
	      eg_copy.add_ear(u,l,i);
	      AddEar(eg_copy, remainingFreeNodes-i, remainingEdges, i, result);
	    }
      else 
	AddNewGraph(eg,result);
    }

    void AddEdge_level3 (FaceGraph eg, int freeEdges, int un, int ln, std::vector<std::pair<Graph<>,IncidenceMatrix<> > >& result, bool fromBelow = 0 ) {
      
#ifdef BIRKHOFF_DEBUG
      if ( fromBelow )
	cout << "[AddEdge_level3] Entering from [AddEdge_level2], have to add " << freeEdges << " edges to graph " << endl << eg.upper() << endl;
#endif

      if ( freeEdges == 0 ) {
	AddNewGraph(eg,result);
      } else {
	for ( int i = un; i < eg.size(); ++i ) {
	  int n = (*static_cast<Set<int> >(eg.lower(i)).rbegin());
	  for ( int j = n+1; j < eg.size(); ++j ) {
	    if ( !eg.upper(i).contains(j) ) {
	      FaceGraph copy_eg(eg);
	      copy_eg.add_edge(i,j);	      
	      AddEdge_level3(copy_eg,freeEdges-1, i, j, result, 1);
	    }
	  }
	}
      }
    }

    void AddEdge_level2 (FaceGraph eg, Set<int> unfinished, int nExcessEdges, int freeEdges, std::vector<std::pair<Graph<>,IncidenceMatrix<> > >& result, bool fromBelow = 0 ) {

#ifdef BIRKHOFF_DEBUG
      static int counter = 0;
#endif
      static Vector<IncidenceMatrix<> > VG;

#ifdef BIRKHOFF_DEBUG
      if ( fromBelow )
	cout << "[AddEdge_level2] Entering from [AddEdge_level1], have to add " << freeEdges << " edges to graph " << endl << eg.upper() << endl;
#endif

      if ( unfinished.size() == 0 ) {
	IncidenceMatrix<> G = eg.as_matrix();
	bool newGraph = true;
	for ( Entire<Vector<IncidenceMatrix<> > >::const_iterator egIt = entire(VG); !egIt.at_end(); ++egIt )
	  //	  if ( CallPolymakeFunction("isomorphic", *egIt, G) ) {  # FIXME no direct isomorphism check for incidence matrices implemented
	  if ( graph::isomorphic(*egIt, G) ) {
	    newGraph = false;
	    break;
	  }
	if ( newGraph ) {
	  VG |= G;
	  if ( freeEdges > 0 ) {
#ifdef BIRKHOFF_DEBUG
	    cout << "[AddEdge_level2] recursion number : " << ++counter << endl;
#endif
	    AddEdge_level3(eg,freeEdges, 0, 0, result, 1);
	  } else {
	    AddNewGraph(eg,result);
	  }
	}
      } else {
	int i = *(unfinished.begin());
	Set<int> copy_unfinished(eg.lower(i).size() < 2 ? unfinished : unfinished-i );
	for ( int j = 0; j < eg.size(); ++j ) {
	  if ( !eg.lower(i).contains(j) ) {
	    FaceGraph copy_eg(eg);
	    copy_eg.add_edge(j,i);
	    AddEdge_level2( copy_eg, copy_unfinished, nExcessEdges, freeEdges, result);
	  }
	}
      }
    }

    // add three edges to all nodes in the upper level
    // we assume that adjacencies in the lower level come in lex order, 
    // and at no point there are gaps in the vector of nodes with positive degree in the lower level
    // un is the current node in the upper level
    // ln1, ln2, ln3 are the neighbors of the nodes added in the previous recursion
    // freeEdges is the remaining number of edges that can be distributed freely (in level3).
    // lastUsedNode is the index of the last node with deg > 0 in the lower level
    void AddEdge_level1 (FaceGraph eg, int nExcessEdges, int freeEdges, int lastUsedNode, 
			 int un, int ln1, int ln2, int ln3, std::vector<std::pair<Graph<>,IncidenceMatrix<> > >& result) {
#ifdef BIRKHOFF_DEBUG
      cout << "[AddEdge_level1] entering" << endl;
      static int counter = 0;
#endif

      if ( un == eg.size() ) {
#ifdef BIRKHOFF_DEBUG
	cout << "[AddEdge_level1] recursion number : " << ++counter << endl;
#endif
	for ( int i=ln1; i < eg.size();  ++i )
	  if ( eg.lower(i).size() < 3 )
	    freeEdges -= 3-eg.lower(i).size(); 
	if ( freeEdges >= 0 ) {
	  Set<int> unfinished;
	  for ( int i=0; i < eg.size(); ++i )
	    if ( eg.lower(i).size() < 3 )
	      unfinished+=i;
	  if ( unfinished.size() > 0 ) {
	    AddEdge_level2( eg, unfinished, nExcessEdges, freeEdges, result, 1);
	  } else {
	    AddEdge_level3(eg, freeEdges, 0, 0, result);
	  }
	}
      } else {
	for ( int i = ln1; i <= std::min ( lastUsedNode+1, eg.size()-1 ); ++i )
	  for ( int j = i==ln1 ? ln2 : i+1; j <= std::min ( std::max ( i+1, lastUsedNode+1 ), eg.size()-1 ); ++j )
	    for ( int k = ( i== ln1 && j==ln2 ) ? ln3 : j+1; k <= std::min ( std::max ( lastUsedNode+1, j+1), eg.size()-1 ); ++k ) {
	      int copy_freeEdges = freeEdges;
	      for ( int l = ln1; l < i; ++l )
		if ( eg.lower(l).size() < 3 )
		  copy_freeEdges -= 3-eg.lower(l).size();
	      if ( copy_freeEdges >= 0 ) {
		FaceGraph copy_eg(eg);
		copy_eg.add_edge(un,i);
		copy_eg.add_edge(un,j);
		copy_eg.add_edge(un,k);
		AddEdge_level1(copy_eg, nExcessEdges, copy_freeEdges, std::max(lastUsedNode, k), un+1, i, j, k, result);
	      }
	    }
      }
    }


    // ----------------------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------
    // construction, if we have a node of deg 2
    // ------------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------------

    void AddRemainingEdges ( FaceGraph & fg,                // the current graph 
			      int n,                         // its number of nodes per layer
			      int d,                         // dimension of the face
			      int k,                         // number of nodes of deg 2 per layer (so the last k cannot be touched)
			      int u_start , int l_start,     // the first node in upper/lower layer we can add to
			      std::vector<std::pair<Graph<>,IncidenceMatrix<> > >& result // collect the completed graphs
			      ) {

#ifdef BIRKHOFF_DEBUG
      cout << "[AddRemainingEdges] entering" << endl;
#endif

      if ( fg.n_missing_edges() > 0 ) {
	for ( int i = u_start; i < n-k; ++i ) {
	  int ls = i == u_start ? l_start : 0;
	  if ( ls < fg.u_minimal_lower_index(i) ) ls = fg.u_minimal_lower_index(i)+1;
	  for ( int j = ls; j < n-k; ++ j ) {
	    if ( !fg.edge(i,j) ) { 
	      FaceGraph copy_fg(fg);
	      copy_fg.add_edge(i,j);
	      AddRemainingEdges(copy_fg, n, d, k, i, j, result);
	    }
	  }
	}
      } else {
#ifdef BIRKHOFF_DEBUG
	cout << "found a new graph: " << endl << fg << endl;
#endif
	AddNewGraph(fg,result);
      }

    }


    bool FillDeg3NodesLower ( FaceGraph & fg,                // the current graph 
			      int n,                         // its number of nodes per layer
			      int d,                         // dimension of the face
			      int k,                         // number of nodes of deg 2 per layer (so the last k cannot be touched)
			      int u_start , int l_start,     // the first node in upper/lower layer we can add to
			      int missing_edges_on_lower,    // how many we should still add
			      std::vector<std::pair<Graph<>,IncidenceMatrix<> > >& result // collect the completed graphs
			      ) {

#ifdef BIRKHOFF_DEBUG
      cout << "[FillDeg3NodesLower] entering with data" << endl
	   << "n_deg_2: " << k << endl 
	   << "u_start: " << u_start << endl << "l_start: " << l_start << endl
	   << "missing_edges_on_lower: " << missing_edges_on_lower << endl
	   << "total missing edges: " << fg.n_missing_edges() << endl;
#endif

      if ( fg.n_missing_edges() < missing_edges_on_lower )
	return false;

      // fill deg 3 nodes in the lower layer
      for ( int j = l_start; j < n-k; ++j ) {
	if ( fg.lower_degree(j) < 3 ) {
	  int us = j == l_start ? u_start : 0;
	  for ( int i = us; i < n-k; ++i ) {
	    if ( !fg.edge(i,j)  && j >= fg.u_minimal_lower_index(i) ) {
	      FaceGraph copy_fg(fg);
	      copy_fg.add_edge(i,j);
	      copy_fg.set_l_min_upper_index(j,i);
	      FillDeg3NodesLower(copy_fg, n, d, k, i, j, missing_edges_on_lower-1, result);
	    }
	  }
	}
      }

      if  ( missing_edges_on_lower == 0 ) 
	AddRemainingEdges(fg, n, d, k, 0, 0, result);

      return true;      
    }



    bool FillDeg3NodesUpper ( FaceGraph & fg,              // the current graph 
			      int n,                       // its number of nodes per layer
			      int d,                       // dimension of the face
			      int k,                       // number of nodes of deg 2 per layer (so the last k cannot be touched)
			      int u_start , int l_start,   // the first node in upper/lower layer we can add to
			      int missing_edges_on_upper, int missing_edges_on_lower,  // how many we should still add
			      std::vector<std::pair<Graph<>,IncidenceMatrix<> > >& result                           // collect the completed graphs
			      ) {

#ifdef BIRKHOFF_DEBUG
      cout << "[FillDeg3NodesUpper] entering with data" << endl
	   << "n_deg_2: " << k << endl 
	   << "u_start: " << u_start << endl << "l_start: " << l_start << endl
	   << "missing_edges_on_upper: " << missing_edges_on_upper << endl
	   << "missing_edges_on_lower: " << missing_edges_on_lower << endl
	   << "total missing edges: " << fg.n_missing_edges() << endl;
#endif

      // check whether we can still have enough edges to bring all non-deg2-nodes to degree at least 3
      if ( fg.n_missing_edges() < missing_edges_on_upper || fg.n_missing_edges() < missing_edges_on_lower ) 
	return false;

      // fill deg 3 nodes in the upper layer
      for ( int i = u_start; i < n-k; ++i ) {
#ifdef BIRKHOFF_DEBUG
	    cout << "testing upper node " << i  << endl;
#endif
	    if ( fg.upper_degree(i) < 3 ) {           // are we at a node that is not yet filled?
	      int ls = i == u_start ? l_start : 0;    // if we are at the same node we have considered in the 
                                                      // previous step, then we can assume that we add on a lex greater node
	      if ( ls < fg.min_u_attach(i) ) {        // if nodes have the same number of deg2 neighbors, then we can assume
                                                      // that the one with the smaller neighbor comes first
	    ls = fg.min_u_attach(i);
#ifdef BIRKHOFF_DEBUG
	    cout << "replacing lower start for upper " << i << " by " << fg.min_u_attach(i) << endl;
#endif
	  }

	      // now we loop over the allowed nodes in the lower layer
	  for ( int j = ls; j < n-k; ++ j ) {
#ifdef BIRKHOFF_DEBUG
	    cout << "testing edge (" << i << ", " << j << ")" << endl;
#endif
	    if ( !fg.edge(i,j) ) {             // if not present, prepare a new recursion step
	      FaceGraph copy_fg(fg);
	      copy_fg.add_edge(i,j);
	      copy_fg.set_u_min_lower_index(i,j);
	      if ( fg.lower_degree(j) < 3 )    // do we also fill a node on the lower layer?
		FillDeg3NodesUpper(copy_fg, n, d, k, i, j, missing_edges_on_upper-1, missing_edges_on_lower-1, result);
	      else
		FillDeg3NodesUpper(copy_fg, n, d, k, i, j, missing_edges_on_upper-1, missing_edges_on_lower, result);
#ifdef BIRKHOFF_DEBUG
	      cout << "returned to " << fg.n_missing_edges() << " missing edges" << endl;
#endif	      
	    }
	  }
	}
      }

      // return if we did not manage to fill all upper nodes
      if ( missing_edges_on_upper > 0 ) 
	return false;

      // otherwise start filling the lower layer
      FillDeg3NodesLower(fg, n, d, k, 0, 0, missing_edges_on_lower, result);

      return true;
    }




    void AddDeg2Loop ( FaceGraph & fg, int n, int d,
		       int k,                        // total number of nodes of deg 2
		       int c,                        // first free node for deg 2
		       int u_attach, int l_attach,   // nodes where the last deg2-loop has been added
		       std::vector<std::pair<Graph<>,IncidenceMatrix<> > >& result) {

#ifdef BIRKHOFF_DEBUG
      cout << "[AddDeg2Loop] entering with data" << endl << "total n_deg_2: " << k << endl << "first free deg 2: " << c << endl
	   << "last used upper: " << u_attach << endl << "last_used_lower: " << l_attach << endl
	   << "and graph" << endl << fg << endl;
#endif

      for ( int i = u_attach; i < n-k; ++i ) {    // loop over all nodes in the upper layer 
	                                          // where the next deg-2-loop can be attached

	int lstart = i == u_attach ? l_attach : 0;
	for ( int j = lstart; j < n-k; ++j ) {  // loop over all nodes in the lower layer
	                                        // where the next deg-2-loop can be attached
	                                        // if the upper node has already been used, we can
                                                // assume that the next loop does not end on a node with smaller index


	  FaceGraph copy_fg(fg);
	  copy_fg.add_edge(i,c);
	  copy_fg.add_edge(c,j);
	  copy_fg.add_edge(c,c);
	  copy_fg.set_u_deg2(i,true);
	  copy_fg.set_u_deg2(j,true);
	  if ( !copy_fg.edge(i,j) ) 
	    copy_fg.add_edge(i,j);

	  if ( c < n-1 ) {             // still a free deg-2-node
	    AddDeg2Loop(copy_fg,n,d, k, c+1, i, j, result);
	  } else {
#ifdef BIRKHOFF_DEBUG
	    cout << "new face graph after deg2: " << endl << copy_fg << endl;
#endif

	    // check that we can actually still fill it
	    int ud = 0;
	    int ld = 0;
	    for ( int h = 0; h < n-k; ++h ) {
	      ud += std::max(3-copy_fg.upper_degree(h),0);
	      ld += std::max(3-copy_fg.lower_degree(h),0);
	    }

	    if ( ud <= copy_fg.n_missing_edges() && ld <= copy_fg.n_missing_edges() ) {
	      FillDeg3NodesUpper ( copy_fg, n, d, k, 0, 0 , ud, ld, result );
	    } else {
#ifdef BIRKHOFF_DEBUG
	      cout << "no continuation possible: " << endl << "missing edges: " << copy_fg.n_missing_edges() << endl
		   << "missing in upper level: " << ud << endl 
		   << "missing in lower level: " << ld << endl;
#endif
	    }
	  }
	}      
      }
    }


    Array<std::pair<Graph<>, IncidenceMatrix<> > > generate_face_graphs(int n, int d) { 

      std::vector<std::pair<Graph<>,IncidenceMatrix<> > > result;
      FaceGraph G(n,d); 

      // if d=1, return the cycle of length 4 if n=2 and nothing otherwise
      if ( d == 1 ) {
	if ( n != 2 )
	  return result;
	else {
	  G.add_edge(0,0);
	  G.add_edge(0,1);
	  G.add_edge(1,0);
	  G.add_edge(1,1);
	  
	  IncidenceMatrix<> I(2,2);
	  I(0,0)=1; I(1,1)=1;

	  std::pair<Graph<>, IncidenceMatrix<> > p(G.as_graph(),I);
	  result.push_back(p);
	  return result;
	}
      }

      int min_nodes_of_deg2 = n-d+1 > 0 ? n-d+1 : 0;
      int max_nodes_of_deg2 = n-1;

      // we deal first with graphs that do not have a node of degree 2
      if ( min_nodes_of_deg2 == 0 ) {
	min_nodes_of_deg2++;

	G.add_edge(0,0);
	G.add_edge(0,1);
	G.add_edge(0,2);
	AddEdge_level1(G,d-1-n,d-1-n,2,1,0,1,2,result);
	
      }

      // now we look at graphs with at least one pair of nodes of degree 2

      for ( int k = min_nodes_of_deg2; k <= max_nodes_of_deg2; ++k ) {  // fix the number of deg 2 nodes	
	FaceGraph H(n,d,k); 
#ifdef BIRKHOFF_DEBUG
	cout << "[generate_face_graphs] working on deg-2-graphs with " << k << " nodes of deg 2" << endl;
#endif	

	H.add_edge(0,0);      // we add the first path of length 3 to the graph
	H.add_edge(0,n-k);
	H.add_edge(n-k,0);
	H.add_edge(n-k,n-k);
	H.set_u_deg2(0,true);
	H.set_l_deg2(0,true);

#ifdef BIRKHOFF_DEBUG
	cout << "[generate_face_graphs] currently missing edges " << H.n_missing_edges() << endl;
#endif	
	if ( k > 1 ) {
	  AddDeg2Loop ( H, n, d, k,   // the number of nodes of deg 2
			n-k+1,        // the first one of the last k we can add to 
			0,            // the first possible node in the upper layer we can attach to
			0,            // the first possible node in the lower layer we can attach to
			result);     
	} else {
	  int ud = 0;
	  int ld = 0;
	  for ( int h = 0; h < n-k; ++h ) {
	    ud += std::max(3-H.upper_degree(h),0);
	    ld += std::max(3-H.lower_degree(h),0);
	  }

	  if ( ud <= H.n_missing_edges() && ld <= H.n_missing_edges() ) {
	    FillDeg3NodesUpper ( H, n, d, k, 0, 0 , ud, ld, result );
	  }
	}
      }

      return result;
    }

    UserFunction4perl(" ", &generate_face_graphs, "generate_face_graphs($,$)");

  }
}
