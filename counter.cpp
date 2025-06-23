#include <iostream>
#include <vector>

using namespace std;
using ll = long long;

int main() {
    int n;
    cin >> n;

    vector<pair<int,int>> vec{{-1,-1},{0,2},{0,3},{0,4}};
    int eta = 2, phi = 3;

    int j = phi, k = 0, v = 0, count = 1, q_start = 1;
    vector<tuple<int,int,int>> list;
    for(int q = 2; q <= (int)vec.size(); q++) {
        if(v != vec[q].first) {
            if(count > 0) {
                list.push_back(make_tuple(v, j, k));
            }
            j = phi, k = 0, v = vec[q].first, count = 1, q_start = q;
        } else if(vec[q].second >= ((k + 1) << j)) {
            if(count > 0) {
                list.push_back(make_tuple(v, j, k));
            }
            while(((k + 2) << j) < vec[q].second) {
                j++, k /= 2;
            }
            k++, count = 1, q_start = q;
        } else {
            count++;
            if(count >= eta) {
                while(j > 0 && ((k + 1) << j) > vec[q].second) {
                    j--, k *= 2;
                }
                if(((k + 1) << j) <= vec[q].second) {
                    if(vec[q_start].second < ((k + 1) << j)) {
                        list.push_back(make_tuple(v, j, k));
                        while(vec[q_start].second < ((k + 1) << j)) {
                            q_start++, count--;
                        }
                    }
                    k++;
                }
            }
        }
    }
    list.push_back(make_tuple(v, j, k));

    for(auto [v, j, k] : list) {
        cout << v << " " << j << " " << k << "\n";
    }


    return 0;
}