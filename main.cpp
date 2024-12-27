#include <iostream>
#include <string>
#include "KDTree.h"
using namespace std;

int main()
{
    //vector<vector<int>> points{{ 1,9 }, { 2,3 }, { 4,1 }, { 3,7 }, { 5,4 }, { 6,8 }, { 7,2 }, { 9,8 }, { 7,9 }, { 9,6 }};
    vector<vector<int>> points{{ 1,2,1 }, { 1,1,3 }, { 2,4,3 }, { 0,0,2 }, { 1,0,0 }, { 0,3,0 }, { 2,2,2 }};
    KDTree<int, 3> tree(points);
	tree.println(KDTree<int, 3>::INORDER);
	tree.printTreeByLevel();
    vector<double> dist;
    auto res = tree.kNearestNeighbors({ 1,1,1 }, 2, &dist);
    cout << "\nK-nearest neighbors of point (1, 1, 1) results (2 points):\n";
    for (int i = 0 ;i<(int) res.size();i++)
        cout << tree.point2string(res[i]) << ": "<< dist[i] << "\n";
	return 0;
}