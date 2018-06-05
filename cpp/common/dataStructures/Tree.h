
#ifndef Tree_h
#define Tree_h

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
    //using TreeT = Tree<T>;
    T content;
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









