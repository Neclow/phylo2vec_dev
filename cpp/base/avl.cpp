#include "avl.hpp"

AVLTree::AVLTree() { root = nullptr; }

unsigned int AVLTree::getHeight(Node *node) { return node ? node->height : 0; }

unsigned int AVLTree::getSize(Node *node) { return node ? node->size : 0; }

void AVLTree::update(Node *node) {
    if (node) {
        node->height =
            1 + std::max(getHeight(node->left), getHeight(node->right));
        node->size = 1 + getSize(node->left) + getSize(node->right);
    }
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

Node *AVLTree::leftRotate(Node *x) {
    Node *y = x->right;
    Node *T2 = y->left;
    y->left = x;
    x->right = T2;
    update(x);
    update(y);
    return y;
}

int AVLTree::getBalance(Node *node) {
    return node ? getHeight(node->left) - getHeight(node->right) : 0;
}

Node *AVLTree::balance(Node *node) {
    int balance = getBalance(node);

    // Left Left Case
    if (balance > 1) {
        if (getBalance(node->left) >= 0) {
            return rightRotate(node);
        } else { // Left Right Case
            node->left = leftRotate(node->left);
            return rightRotate(node);
        }
    }

    // Right Right Case
    if (balance < -1) {
        if (getBalance(node->right) <= 0) {
            return leftRotate(node);
        } else { // Right Left Case
            node->right = rightRotate(node->right);
            return leftRotate(node);
        }
    }

    return node;
}

void AVLTree::insert(int index, Pair value) {
    this->root = insertByIndex(this->root, value, index);
}

Node *AVLTree::insertByIndex(Node *node, Pair value, int index) {
    if (!node) {
        return new Node(value);
    }

    int left_size = getSize(node->left);

    if (index <= left_size) { // Insert in the left subtree
        node->left = insertByIndex(node->left, value, index);
    } else { // Insert in the right subtree
        node->right = insertByIndex(node->right, value, index - left_size - 1);
    }

    update(node);
    return balance(node);
}

Pair AVLTree::lookup(Node *node, int index) {
    if (!node) {
        return {0, 0}; // Index out of bounds
    }

    int left_size = getSize(node->left);

    if (index < left_size) { // Search in left subtree
        return lookup(node->left, index);
    } else if (index == left_size) { // Found the node
        return node->value;
    } else { // Search in right subtree
        return lookup(node->right, index - left_size - 1);
    }
}

std::vector<Pair> AVLTree::inorderTraversal(Node *node) {
    std::vector<Pair> result;
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

std::vector<Pair> AVLTree::getPairs() { return inorderTraversal(root); }
