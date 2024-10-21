#ifndef AVL_HPP
#define AVL_HPP

#include <algorithm>
#include <stack>

#include "core.hpp"

struct Node {
    Pair value;
    Node *left;
    Node *right;
    unsigned int height; // Height of the node
    unsigned int size;   // Number of nodes in the subtree

    Node(Pair val)
        : value(val), left(nullptr), right(nullptr), height(1), size(1) {}
};

class AVLTree {
  public:
    Node *root;

    AVLTree();

    void insert(int index, Pair value);

    Pair lookup(Node *node, int index);

    std::vector<std::array<unsigned int, 2>> getPairs();

    unsigned int getHeight(Node *node);

    unsigned int getSize(Node *node);

    int getBalance(Node *node);

  private:
    void update(Node *node);

    Node *leftRotate(Node *x);

    Node *rightRotate(Node *y);

    Node *balance(Node *node);

    Node *insertByIndex(Node *node, Pair value, int index);

    std::vector<Pair> inorderTraversal(Node *node);
};

#endif // AVL_HPP