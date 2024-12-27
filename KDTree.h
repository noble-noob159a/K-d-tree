#ifndef KDTREE_H
#define KDTREE_H
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <queue>
#include <iomanip>
using namespace std;

class DimensionMismatchException : public std::exception {
private:
    string e;
public:
    DimensionMismatchException(int m, int k) {
        this->e = "Data dimension mismatch. Expect: " + to_string(k) + ", got: " + to_string(m);
    }
    const char* what() const throw () {
        stringstream os;
        os << e;
        return os.str().c_str();
    }
};


//T must support operator- for calulating distances (return in type of double or int) and k must be > 0
//No identical points
template<class T, int k>
class KDTree{
public:
    class Node;
protected:
    Node* root;
    int count;
    double imbalanceThreshold;
    int (*comparator)(const T& lhs,const T& rhs);    //custom type comparator return: -1 if less than; 0 if equal; 1 if greater than
public:
    KDTree(Node* root = 0, double imbalanceThreshold = 2.0,int (*comparator)(const T& lhs,const T& rhs) = 0): root(root),count(0),comparator(comparator),imbalanceThreshold(imbalanceThreshold) {}
    KDTree(const vector<vector<T>>& pointList, double imbalanceThreshold = 2.0, int (*comparator)(const T& lhs,const T& rhs) = 0) : root(0), count(0), comparator(comparator), imbalanceThreshold(imbalanceThreshold)
    {
        build(pointList);
    }
    KDTree(const KDTree& tree): root(0),count(0),comparator(0),imbalanceThreshold(2)
    {
        copyFrom(tree);
    }
    KDTree& operator=(const KDTree& tree)
    {
        if (this != &tree) {
            clear();
            copyFrom(tree);
        }
        return *this;
    }
    ~KDTree()
    {
        clear();
    }
    void add(const vector<T>& point, bool b_rebuild = 1) {
        if ((int)point.size() != k) throw DimensionMismatchException(point.size(), k);
        root = addRecursive(root, point, 0);
        ++count;
        if (getImbalanceRatio() > imbalanceThreshold && b_rebuild) rebuild();
    }
    void build(const vector<vector<T>>& pointList)
    {
        clear();
        root = buildRecursive(pointList, 0);
        count = pointList.size();
    }
    vector<vector<T>> rangeSearch(const vector<T>& lowerBound, const vector<T>& upperBound) {
        if ((int)lowerBound.size() != k || (int)upperBound.size() != k) {
            throw DimensionMismatchException(lowerBound.size(), k);
        }
        vector<vector<T>> results;
        rangeSearchRecursive(root, lowerBound, upperBound, 0, results);
        return results;
    }
    vector<T> nearestNeighbor(const vector<T>& target, double* outDist = 0) {
        if ((int)target.size() != k) throw DimensionMismatchException(target.size(), k);

        vector<T> bestPoint;
        double bestDist = numeric_limits<double>::max();

        nearestNeighborRecursive(root, target, 0, bestPoint, bestDist);
        if (outDist) *outDist = sqrt(bestDist);
        return bestPoint;
    }
    vector<vector<T>> kNearestNeighbors(const vector<T>& target, int kNeighbors, vector<double>* outDist = 0) {
        if ((int)target.size() != k) throw DimensionMismatchException(target.size(), k);

        priority_queue<pair<double, vector<T>>> maxHeap; 

        kNearestNeighborsRecursive(root, target, 0, kNeighbors, maxHeap);

        vector<vector<T>> results;
        vector<double> dist;
        while (!maxHeap.empty()) {
            results.push_back(maxHeap.top().second);
            dist.push_back(sqrt(maxHeap.top().first));
            maxHeap.pop();
        }
        reverse(results.begin(), results.end());
        reverse(dist.begin(), dist.end());
        if (outDist) *outDist = dist;
        return results;
        
    }

    int size() { return count; }
    int height() {
        return heightRec(root);
    }
    void clear()
    {
        removeInternalData(root);
        root = nullptr;
        count = 0;
    }
    bool empty() { return count == 0; }
    bool remove(const vector<T>& point, bool b_rebuild = 1) {
        if ((int)point.size() != k) throw DimensionMismatchException(point.size(), k);

        bool removed = false;
        root = removeRecursive(root, point, 0, removed);

        if (removed) {
            --count;
        }
        if (getImbalanceRatio() > imbalanceThreshold && b_rebuild) rebuild();
        return removed;
    }
    bool search(const vector<T>& point)
    {
        return searchRecursive(root, point, 0);
    }
    string toString(int mode = 2, string(*point2str)(T&) = 0) 
    {
        string mark(50, '=');
        stringstream os;
        os << mark << endl;
        os << "Tree Info:" << endl;
        os << left << setw(22) << "Height: " << height() <<endl;
        string rootStr = root ? point2string(root->point, point2str) : "";
        os << left << setw(22) << "Root: " << rootStr << endl;
        vector<vector<T>> points;
        if (mode == 0) { points = bfs(); os << left << setw(22) << "BFS Traversal: ";}
        else if (mode == 1) { points = preOrder(); os << left << setw(22) << "Pre-Order Traversal: "; }
        else if (mode == 2) { points = inOrder(); os << left << setw(22) << "In-Order Traversal: ";  }
        else if (mode == 3) { points = postOrder(); os << left << setw(22) << "Post-Order Traversal: ";}
        string tmp = "";
        for (auto point : points)
            tmp += point2string(point, point2str) + ", ";
        if (tmp.size() > 1)
        {
            tmp.pop_back(); tmp.pop_back();
        }
        os << tmp << endl;
        os << mark << endl;
        return os.str();
    }
    string point2string(const vector<T>& point, string(*point2str)(T&) = 0)
    {     
        string tmp = "(";
        for (auto e : point)
        {
            if (point2str)
                tmp += point2str(e) + ", ";
            else
            {
                stringstream ss;
                ss << e << ", ";
                tmp += ss.str();
            }
        }
        if (tmp.size() > 1)
        {
            tmp.pop_back(); tmp.pop_back();
        }
        tmp += ")";
        return tmp;
    }
    void println(int mode = 2, string(*point2str)(T&) = 0) {
        cout << toString(mode,point2str) << endl;
    }
    vector<vector<T>> bfs() {
        vector<vector<T>> result;
        if (!root) return result;

        queue<Node*> q;
        q.push(root);

        while (!q.empty()) {
            Node* current = q.front();
            q.pop();

            result.push_back(current->point);

            if (current->pLeft) q.push(current->pLeft);
            if (current->pRight) q.push(current->pRight);
        }

        return result;
    }
    vector<vector<T>> preOrder() {
        vector<vector<T>> result;
        preOrderRecursive(root, result);
        return result;
    }
    vector<vector<T>> inOrder() {
        vector<vector<T>> result;
        inOrderRecursive(root, result);
        return result;
    }
    vector<vector<T>> postOrder() {
        vector<vector<T>> result;
        postOrderRecursive(root, result);
        return result;
    }
    void printTreeByLevel(string(*point2str)(T&) = 0) {
        if (!root) {
            cout << "Tree is empty." << endl;
            return;
        }
        string mark(50, '=');
        queue<pair<Node*, int>> q;
        q.push({ root, 0 });
        int currentLevel = 0;
        cout << mark << endl;
        cout << "Tree By Level (left to right):" << endl;
        cout << left << setw(9) <<"Level 0:";

        while (!q.empty()) {
            Node* node = q.front().first;
            int level = q.front().second;
            q.pop();

            if (level > currentLevel) {
                currentLevel = level;
                string tmp = "\nLevel " + to_string(currentLevel) + ": ";
                cout << left << setw(9) << tmp;
            }

            cout << point2string(node->point, point2str);
            cout << ",  ";

            if (node->pLeft) {
                q.push({ node->pLeft, level + 1 });
            }
            if (node->pRight) {
                q.push({ node->pRight, level + 1 });
            }
        }
       cout << endl;
       cout << mark << endl;
    }
    
public:
    class Node
    {
    public:
        vector<T> point;
        Node* pLeft;
        Node* pRight;

        Node(const vector<T>& point, Node* left=0, Node* right=0): pLeft(left), pRight(right)
        {
            int dim = point.size();
            if (dim != k) throw DimensionMismatchException(dim, k);
            this->point = point;
        }
        ~Node() {}
    };

protected:
    double getImbalanceRatio() {
        if (!root) return 0.0; 

        int h = this->height();
        double expectedHeight = floor(log2(count)) + 1;
        double imbalanceRatio = (double)h / expectedHeight;
        return imbalanceRatio;
    }
    void rebuild() {
        vector<vector<T>> allPoints = bfs();
        build(allPoints);
    }
    Node* addRecursive(Node* node, const vector<T>& point, int depth) {
        if (!node) {
            return new Node(point);
        }
        int dim = depth % k;  // Current dimension
        if (compare(point[dim], node->point[dim], comparator) < 0) {
            node->pLeft = addRecursive(node->pLeft, point, depth + 1);
        }
        else {
            node->pRight = addRecursive(node->pRight, point, depth + 1);
        }
        return node;
    }
    Node* buildRecursive(const vector<vector<T>>& points, int depth) {
        if (points.empty()) return nullptr;
        int dim = depth % k;
        vector<vector<T>> sortedPoints = points;
        std::sort(sortedPoints.begin(), sortedPoints.end(), [&](const vector<T>& a, const vector<T>& b) {
            return compare(a[dim], b[dim], this->comparator) < 0;
         });

        int median = (int)sortedPoints.size() / 2;
        Node* node = new Node(sortedPoints[median]);
       for (int i = median - 1; i >= 0; i--)
        {
            if (!compare(sortedPoints[i][dim],node->point[dim],this->comparator)) { node->point = sortedPoints[i]; median = i; }
        }
        vector<vector<T>> leftPoints(sortedPoints.begin(), sortedPoints.begin() + median);
        vector<vector<T>> rightPoints(sortedPoints.begin() + median + 1, sortedPoints.end());

        node->pLeft = buildRecursive(leftPoints, depth + 1);
        node->pRight = buildRecursive(rightPoints, depth + 1);
        return node;
    }
    Node* removeRecursive(Node* current, const vector<T>& point, int depth, bool& removed) {
        if (!current) {
            removed = false;
            return nullptr;
        }

        int axis = depth % k; 
        if (current->point == point) {
            removed = true;

            if (!current->pLeft && !current->pRight) {
                delete current;
                return nullptr;
            }

            if (current->pRight) {
                Node* minNode = findMin(current->pRight, axis, depth + 1);
                current->point = minNode->point;
                current->pRight = removeRecursive(current->pRight, minNode->point, depth + 1, removed);
            }
            else {
                Node* minNode = findMin(current->pLeft, axis, depth + 1);
                current->point = minNode->point;
                current->pRight = removeRecursive(current->pLeft, minNode->point, depth + 1, removed);
                current->pLeft = nullptr;
            }
            return current;
        }

        if (compare(point[axis], current->point[axis], comparator) < 0) {
            current->pLeft = removeRecursive(current->pLeft, point, depth + 1, removed);
        }
        else {
            current->pRight = removeRecursive(current->pRight, point, depth + 1, removed);
        }
        return current;
    }
    Node* findMin(Node* current, int axis, int depth) {
        if (!current) return nullptr;

        int currentAxis = depth % k;

        if (currentAxis == axis) {
            if (current->pLeft) {
                return findMin(current->pLeft, axis, depth + 1);
            }
            return current;
        }

        Node* leftMin = findMin(current->pLeft, axis, depth + 1);
        Node* rightMin = findMin(current->pRight, axis, depth + 1);

        Node* minNode = current;
        if (leftMin && compare(leftMin->point[axis],minNode->point[axis],comparator) < 0) {
            minNode = leftMin;
        }
        if (rightMin && compare(rightMin->point[axis],minNode->point[axis],comparator)<0) {
            minNode = rightMin;
        }

        return minNode;
    }
    void preOrderRecursive(Node* node, vector<vector<T>>& result) {
        if (!node) return;
        result.push_back(node->point);
        preOrderRecursive(node->pLeft, result);
        preOrderRecursive(node->pRight, result);
    }
    void inOrderRecursive(Node* node, vector<vector<T>>& result) {
        if (!node) return;
        inOrderRecursive(node->pLeft, result);
        result.push_back(node->point);
        inOrderRecursive(node->pRight, result);
    }
    void postOrderRecursive(Node* node, vector<vector<T>>& result) {
        if (!node) return;
        postOrderRecursive(node->pLeft, result);
        postOrderRecursive(node->pRight, result);
        result.push_back(node->point);
    }
    int heightRec(Node* node) {
        if (!node) return 0;
        return 1 + std::max(heightRec(node->pLeft), heightRec(node->pRight));
    }
    bool searchRecursive(Node* node, const vector<T>& point, int depth) {
        if (!node) return false;

        int axis = depth % k;
        int cmp = compare(node->point[axis], point[axis], comparator);
        if (cmp == 0) {
            if (node->point == point) return true;   //check all dim
        }
        if (cmp > 0) return searchRecursive(node->pLeft, point, depth + 1);
        return searchRecursive(node->pRight, point, depth + 1);
    }
    void copyFrom(const KDTree& tree) {
        comparator = tree.comparator;
        count = tree.count;
        imbalanceThreshold = tree.imbalanceThreshold;
        root = copyRecursive(tree.root);
    }
    Node* copyRecursive(Node* node) {
        if (!node) return nullptr;
        Node* newNode = new Node(node->point);
        newNode->pLeft = copyRecursive(node->pLeft);
        newNode->pRight = copyRecursive(node->pRight);
        return newNode;
    }
    void removeInternalData(Node* pNode)
    {
        if (!pNode) return;
        removeInternalData(pNode->pLeft);
        removeInternalData(pNode->pRight);
        delete pNode;
    }
    void rangeSearchRecursive(Node* node, const vector<T>& lowerBound, const vector<T>& upperBound, int depth, vector<vector<T>>& results) {
        if (!node) return;
        int axis = depth % k;
        bool inside = true;
        for (int i = 0; i < k; ++i) {
            if (compare(node->point[i],lowerBound[i],comparator) < 0 || compare(node->point[i],upperBound[i],comparator) > 0) {
                inside = false;
                break;
            }
        }
        if (inside) results.push_back(node->point);

        if (node->pLeft && compare(node->point[axis],lowerBound[axis],comparator) > 0) {
            rangeSearchRecursive(node->pLeft, lowerBound, upperBound, depth + 1, results);
        }
        if (node->pRight && compare(node->point[axis],upperBound[axis],comparator) <= 0) {
            rangeSearchRecursive(node->pRight, lowerBound, upperBound, depth + 1, results);
        }
    }
    void nearestNeighborRecursive(Node* node, const vector<T>& target, int depth, vector<T>& bestPoint, double& bestDist) {
        if (!node) return;

        int axis = depth % k;

        double dist = 0;
        for (int i = 0; i < k; ++i) {
            dist += (node->point[i] - target[i]) * (node->point[i] - target[i]);
        }

        if (dist < bestDist) {
            bestDist = dist;
            bestPoint = node->point;
        }

        Node* first = compare(target[axis],node->point[axis],comparator) < 0 ? node->pLeft : node->pRight;
        Node* second = compare(target[axis], node->point[axis], comparator) < 0 ? node->pRight : node->pLeft;

        nearestNeighborRecursive(first, target, depth + 1, bestPoint, bestDist);

        if ((target[axis] - node->point[axis]) * (target[axis] - node->point[axis]) < bestDist) {
            nearestNeighborRecursive(second, target, depth + 1, bestPoint, bestDist);
        }
    }
    void kNearestNeighborsRecursive(Node* node, const vector<T>& target, int depth, int kNeighbors, priority_queue<pair<double, vector<T>>>& maxHeap) {
        if (!node) return;

        int axis = depth % k;

        double dist = 0;
        for (int i = 0; i < k; ++i) {
            dist += (node->point[i] - target[i]) * (node->point[i] - target[i]);
        }

        maxHeap.push({ dist, node->point });
        if ((int)maxHeap.size() > kNeighbors) {
            maxHeap.pop();
        }

        Node* first = compare(target[axis], node->point[axis], comparator) < 0 ? node->pLeft : node->pRight;
        Node* second = compare(target[axis], node->point[axis], comparator) < 0 ? node->pRight : node->pLeft;

        kNearestNeighborsRecursive(first, target, depth + 1, kNeighbors, maxHeap);

        if ((target[axis] - node->point[axis]) * (target[axis] - node->point[axis]) < maxHeap.top().first) {
            kNearestNeighborsRecursive(second, target, depth + 1, kNeighbors, maxHeap);
        }
    }
    static int compare(const T& lhs, const T& rhs, int (*comparator)(const T&, const T&) = 0) {
        if (comparator != 0) return comparator(lhs, rhs);
        else {
            if (lhs < rhs) return -1;
            else if (lhs > rhs) return +1;
            else return 0;
        }
    }
public:
    static inline int BFS = 0;
    static inline int PREORDER = 1;
    static inline int INORDER = 2;
    static inline int POSTORDER = 3;
};

#endif
