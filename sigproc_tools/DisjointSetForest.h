/**
 * \file DisjointSetForest.h
 *
 * \ingroup sigproc_tools
 * 
 * \brief Class def header for a class DisjointSetForest
 *
 * @author koh0207
 */

/** \addtogroup sigproc_tools

    @{*/
#ifndef __SIGPROC_TOOLS_DISJOINTSETFOREST_H__
#define __SIGPROC_TOOLS_DISJOINTSETFOREST_H__

#include <vector>
#include <iostream>
#include <string>

namespace sigproc_tools {

  /**
     \class DisjointSetForest
     User defined class DisjointSetForest ... these comments are used to generate
     doxygen documentation!
  */

  class DisjointSetForest{

    int rootNode = 0;
    
    public:
      
      /// Default constructor
      DisjointSetForest(const unsigned int size, const int rootIndex=0) {
        parent.resize(size+1);
        rank.resize(size+1);
        this->size = size+1;
        rootNode = rootIndex;
      }
      
      /// Default destructor
      ~DisjointSetForest(){}

      std::vector<int> parent;
      std::vector<int> rank;
      size_t size = 0;

      void MakeSet();

      void MakeSet(const std::vector<int>& strongEdges);

      void Union(const int x, const int y);

      int Find(const int x);
  };
}

#endif
/** @} */ // end of doxygen group 

