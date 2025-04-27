#ifndef AVL_HPP
#define AVL_HPP

#include "../base/core.hpp"

struct Node {
    Pair value;
    Node *left;
    Node *right;
    unsigned int height;  // Height of the node
    unsigned int size;    // Number of nodes in the subtree

    Node(Pair val) : value(val), left(nullptr), right(nullptr), height(1), size(1) {}
};

class AVLTree {
   public:
    AVLTree();

    Node *getRoot();

    Pairs getPairs();

    void insert(int index, Pair value);

    Pair lookup(Node *node, int index);

   private:
    Node *root;

    int getBalanceOfNode(Node *node);

    unsigned int getHeightofNode(Node *node);

    unsigned int getSizeOfNode(Node *node);

    void update(Node *node);

    Node *leftRotate(Node *x);

    Node *rightRotate(Node *y);

    Node *balance(Node *node);

    Node *insertByIndex(Node *node, int index, Pair value);

    Pairs inorderTraversal(Node *node);
};

#endif  // AVL_HPP