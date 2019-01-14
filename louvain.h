#ifndef __LOUVAIN_H
#define __LOUVAIN_H

// The code below implements the algorithm described in the paper
// "Fast unfolding of communities in large networks"
// Vincent D Blondel et al J. Stat. Mech. (2008) P10008
//
// Implemented by Hyun Min Kang on July 10, 2018

#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <cstdio>

#include "Error.h"

class lnode {
public:
  int32_t id;

  std::vector<lnode*> children;
  
  lnode* parent;
  
  lnode(int32_t _id, lnode* _parent) : id(_id), parent(_parent) {}
  
  void add_child(lnode* child) {
    children.push_back(child);
    child->parent = this;
  }

  std::string get_name() {
    if ( id >= 0 ) {
      if ( children.empty() ) {
	char buf[255];
	sprintf(buf,"%d",id);
	return std::string(buf);
      }
      else {
	error("Internal node with non-negative ID %d", id);
	return std::string(".");	
      }
    }
    else if ( children.empty() ) {
      return std::string(".");
    }
    else {
      std::string name("[");
      for(int32_t i=0; i < (int32_t)children.size(); ++i) {
	if (i > 0) name += ",";
	name += children[i]->get_name();
      }
      name += "]";
      return name;
    }
  }

  void print_summary() {
    notice("# of clusters = %d", (int32_t)children.size());
    for(int32_t i=0; i < (int32_t)children.size(); ++i) {
      notice("Cluster %d : Size %d", i+1, children[i]->size());
    }    
  }

  int32_t size() {
    if ( id >= 0 ) { return 1; }
    else {
      int32_t sum = 0; 
      for(int32_t i=0; i < (int32_t)children.size(); ++i) {
	sum += children[i]->size();
      }
      return sum;
    }        
  }

  void print(int32_t level) {
    for(int32_t i=0; i < level; ++i) {
      printf("..");
    }
    if ( id >= 0 ) { printf(" %d\n", id); }
    else { printf(" *\n"); }

    for(int32_t i=0; i < (int32_t)children.size(); ++i) {
      children[i]->print(level+1);
    }
  }
};

class louvain {
public:
  std::map<lnode*, std::map<lnode*, double> > edges;
  lnode* root;
  int32_t nleafs;
  double  sum_edges;
  
  louvain(int32_t _nleafs) : nleafs(_nleafs), sum_edges(0) {
    root = new lnode(-1, NULL);
    for(int32_t i=0; i < nleafs; ++i) {
      lnode* child = new lnode(i, root);
      root->add_child(child);
    }
  }

  void add_edge(int32_t i, int32_t j, double d) {
    edges[root->children[i]][root->children[j]] = d;
    edges[root->children[j]][root->children[i]] = d;
    sum_edges += (d + d);
  }

  void print() {
    root->print(0);
  }

  void print_summary() {
    root->print_summary();
  }

  // make one pass of Louvain algorithm
  // returns true if anything was mofieid
  bool make_cluster_pass(bool shuffle = true) {
    // disconnect the children from the root
    std::vector<lnode*> nodes = root->children;
    //root->children.clear();

    // assign each node to its own cluster, and calculate marginals
    std::vector<double> ks(nodes.size(),0);
    std::vector<double> sum_ins(nodes.size(), 0);
    std::vector<double> sum_tots(nodes.size(), 0);
    std::vector< std::vector<lnode*> > clusts(nodes.size());

    // assign node indices
    std::vector<int32_t> inodes;
    for(int32_t i=0; i < (int32_t)nodes.size(); ++i) inodes.push_back(i);
    
    // shuffle nodes
    if ( shuffle ) {
      for(int32_t i=0; i < (int32_t)inodes.size(); ++i) {
	int32_t j = i + (rand() % ((int32_t)inodes.size() - i));
	int32_t t = inodes[j];
	inodes[j] = inodes[i];
	inodes[i] = t;
      }
    }

    // make clusters as the originals
    for(int32_t i=0; i < (int32_t)nodes.size(); ++i) {
      lnode* nd = nodes[inodes[i]];
      std::map<lnode*,double>& iedges = edges[nd];
      for(std::map<lnode*,double>::iterator it = iedges.begin();
	  it != iedges.end(); ++it) {
	ks[i] += it->second;
	if ( it->first == nd ) sum_ins[i] += it->second;
	else sum_tots[i] += it->second;
      }
      clusts[i].push_back(nd);
    }

    bool modified = false;

    for(int32_t i=0; i < (int32_t)clusts.size(); ++i) {
      // check if it is an isolated node
      if ( clusts[i].size() == 1 ) {
	// try to assign existing community
	double  max_delta_q = 0;
	int32_t best  = -1;
	std::map<lnode*,double>& iedges = edges[clusts[i][0]];	  
	for(int32_t j=0; j < i; ++j) {
	  double k_in = 0;
	  for(int32_t x = 0; x < (int32_t)clusts[j].size(); ++x) {
	    if ( iedges.find(clusts[j][x]) != iedges.end() ) { // node exists
	      k_in += iedges[clusts[j][x]];
	    }
	  }
	  double delta_q = ( (sum_ins[j] + k_in)/sum_edges - ( sum_tots[j] + ks[i] ) * ( sum_tots[j] + ks[i] ) / sum_edges / sum_edges ) - ( sum_ins[j] / sum_edges - sum_tots[j] * sum_tots[j] / sum_edges / sum_edges - ks[i] * ks[i] / sum_edges / sum_edges );
	  if ( delta_q > max_delta_q ) {
	    max_delta_q = delta_q;
	    best  = j; 
	  }
	}
	
	if ( best >= 0 ) { // assign the node to a cluster
	  // calculate the new sum_ins
	  double new_sum_ins = sum_ins[best];
	  std::map<lnode*,double>& iedges = edges[clusts[i][0]];
	  double new_sum_tots = sum_tots[best] + sum_tots[i];	    
	  for(int32_t x=0; x < (int32_t)clusts[best].size(); ++x) {
	    if ( iedges.find(clusts[best][x]) != iedges.end() ) { // node exists
	      if ( clusts[i][0] == clusts[best][x] ) {
		new_sum_ins += iedges[clusts[best][x]];
	      }
	      else {
		new_sum_ins += (2 * iedges[clusts[best][x]]);
		new_sum_tots -= (2 * iedges[clusts[best][x]]);
	      }
	    }	      
	  }
	  
	  //printf("%s --> %d, %lg\t%lg\t%lg\t%lg\t%lg\n", node2str[clusts[i][0]].c_str(), best, max_delta_q, sum_ins[best], new_sum_ins, sum_tots[best], new_sum_tots);
	  
	  sum_ins[best] = new_sum_ins;
	  sum_tots[best] = new_sum_tots;
	  clusts[best].push_back( clusts[i][0] );
	  
	  sum_ins.erase(sum_ins.begin() + i);
	  sum_tots.erase(sum_tots.begin() + i);
	  clusts.erase(clusts.begin() + i);
	  
	  modified = true;
	  --i;
	}
      }
    }

    if ( modified ) {
      // make new levels
      root->children.clear();
      for(int32_t i=0; i < (int32_t)clusts.size(); ++i) {
	lnode* nd = new lnode(-1, NULL);
	for(int32_t j=0; j < (int32_t)clusts[i].size(); ++j) {
	  nd->add_child(clusts[i][j]);
	}
	root->add_child(nd);
      }

      // calculate between edge distances
      std::vector<lnode*>& c = root->children;
      for(int32_t i=0; i < (int32_t)clusts.size(); ++i) {
	edges[c[i]][c[i]] = sum_ins[i];
	for(int32_t j=0; j < i; ++j) {
	  double d = 0;
	  for(int32_t ii=0; ii < (int32_t)clusts[i].size(); ++ii) {
	    std::map<lnode*,double>& iiedges = edges[clusts[i][ii]];
	    for(int32_t jj=0; jj < (int32_t)clusts[j].size(); ++jj) {
	      if ( iiedges.find(clusts[j][jj]) != iiedges.end() ) {
		d += iiedges[clusts[j][jj]];
	      }
	    }
	  }
	  edges[c[i]][c[j]] = d;
	  edges[c[j]][c[i]] = d;	  
	}
      }
      return true;
    }
    else {
      return false;
    }
  }
};

#endif
