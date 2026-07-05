
#ifndef Tree_h
#define Tree_h

/// @file Tree.h
/// @brief CRTP-based recursive tree template for hierarchical structures.
///
/// The challenge with recursive tree types in C++ is that a Tree<T> needs to contain
/// vector<Tree<T>>, but Tree<T> is incomplete at that point. TreeGen solves this with
/// a two-template-parameter CRTP pattern: TreeGen<T, TreeT> holds the branches as
/// vector<TreeT>, and the concrete types fix TreeT to themselves.
///
/// - Tree<T> inherits TreeGen<T, Tree<T>> — value-based tree, branches stored by value.
/// - PTree<T> inherits TreeGen<T, PTree<T>*> — pointer-based tree with parent back-links,
///   branches stored as pointers. This allows shared subtrees and upward traversal.
///
/// Used for scene graphs, BVH hierarchies, and nested mesh group structures where
/// you need to traverse both down (branches) and up (parent).

#include <vector>

// this is able to generate properly inheriatable trees
template<typename T,typename TreeT>
class TreeGen{ public:
    //using TreeT = Tree<T>;
    T content;
    std::vector<TreeT> branches;
};

// this can be only specialized
template<typename T>
class Tree : public TreeGen<T,Tree<T>> { public:
    T content;
};

template<typename T>
class PTree : public TreeGen<T,PTree<T>*> { public:
    PTree<T>* parrent;
    T content;
    ~PTree(){ for(PTree<T>* br : this->branches) delete br; this->branches.clear(); }
};


/*
template<typename T,typename TreeT>
class Tree{ public:
    //using TreeT = Tree<T>;
    std::vector<TreeT> branches;
    std::vector<T>     leafs;
};
*/

/*
// labeled tree can be made easily by derived class
template<typename T, typename TreeT, typename L>
class LabeledTree : public Tree<T,TreeT> { public:
    //using TreeT = LabeledTree<T,L>;
    L label;
    //std::vector<TreeT> branches;
    //std::vector<T>     leafs;
};
*/

#endif









