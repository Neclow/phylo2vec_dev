#include "avl.hpp"

#include <algorithm>
#include <stack>

AVLTree::AVLTree() { root = nullptr; }

Node *AVLTree::getRoot() { return root; }

Pairs AVLTree::getPairs() { return inorderTraversal(root); }

Pairs AVLTree::inorderTraversal(Node *node) {
    Pairs result;
    std::stack<Node *> stack;
    Node *current = node;

    while (current != nullptr || !stack.empty()) {
        while (current != nullptr) {
            stack.push(current);
            current = current->left;
        }

        current = stack.top();
        stack.pop();
        result.push_back(current->value);

        current = current->right;
    }

    return result;
}

void AVLTree::insert(int index, Pair value) {
    this->root = insertByIndex(this->root, index, value);
}

Node *AVLTree::insertByIndex(Node *node, int index, Pair value) {
    if (!node) {
        return new Node(value);
    }

    int left_size = getSizeOfNode(node->left);

    if (index <= left_size) {  // Insert in the left subtree
        node->left = insertByIndex(node->left, index, value);
    } else {  // Insert in the right subtree
        node->right = insertByIndex(node->right, index - left_size - 1, value);
    }

    update(node);
    return balance(node);
}

void AVLTree::update(Node *node) {
    if (node) {
        node->height = 1 + std::max(getHeightofNode(node->left), getHeightofNode(node->right));
        node->size = 1 + getSizeOfNode(node->left) + getSizeOfNode(node->right);
    }
}

Node *AVLTree::balance(Node *node) {
    int balance = getBalanceOfNode(node);

    // Left Left Case
    if (balance > 1) {
        if (getBalanceOfNode(node->left) >= 0) {
            return rightRotate(node);
        } else {  // Left Right Case
            node->left = leftRotate(node->left);
            return rightRotate(node);
        }
    }

    // Right Right Case
    if (balance < -1) {
        if (getBalanceOfNode(node->right) <= 0) {
            return leftRotate(node);
        } else {  // Right Left Case
            node->right = rightRotate(node->right);
            return leftRotate(node);
        }
    }

    return node;
}

Node *AVLTree::leftRotate(Node *x) {
    Node *y = x->right;
    Node *T2 = y->left;
    y->left = x;
    x->right = T2;
    update(x);
    update(y);
    return y;
}

Node *AVLTree::rightRotate(Node *y) {
    Node *x = y->left;
    Node *T2 = x->right;
    x->right = y;
    y->left = T2;
    update(y);
    update(x);
    return x;
}

Pair AVLTree::lookup(Node *node, int index) {
    if (!node) {
        return {0, 0};  // Index out of bounds
    }

    int left_size = getSizeOfNode(node->left);

    if (index < left_size) {  // Search in left subtree
        return lookup(node->left, index);
    } else if (index == left_size) {  // Found the node
        return node->value;
    } else {  // Search in right subtree
        return lookup(node->right, index - left_size - 1);
    }
}

int AVLTree::getBalanceOfNode(Node *node) {
    return node ? getHeightofNode(node->left) - getHeightofNode(node->right) : 0;
}

unsigned int AVLTree::getHeightofNode(Node *node) { return node ? node->height : 0; }

unsigned int AVLTree::getSizeOfNode(Node *node) { return node ? node->size : 0; }
