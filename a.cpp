#include <functional>
#include <iostream>
#include <vector>

using namespace std;
using ll = long long;

pair<function<ll(ll)>, function<ll(ll)>>
displace(vector<ll> keys, function<ll(ll)> f, function<ll(ll)> g, int r) {
	int n = keys.size();

	vector<pair<ll, ll>> v(n);
    vector<int> cnt(1 << r);
    vector<vector<int>> vals(n + 1);

	for (int i = 0; i < n; i++) {
		v[i] = {f(keys[i]), g(keys[i])};
        cnt[v[i].first]++;
    }

    for(int i = 0; i < (1 << r); i++) {
        vals[cnt[i]].push_back(i);
    }

    for (int i = 1; i < (1 << r); i++) {
        cnt[i] += cnt[i - 1];
    }

    // TODO in place
    vector<pair<ll,ll>> v2(n);
    for(int i = 0; i < n; i++) {
        v2[--cnt[v[i].first]] = v[i];
    }
    v = move(v2);
    // ENDTODO

    vector<int> p(1 << r);
    for(int i = n, p_idx = 0; i >= 0; i--) {
        for(auto val : vals[i]) {
            p[p_idx++] = val;
        }
    }
    
    vector<int> m(1 << r, 0);
    for(int i = 0; i < (1 << r); i++) {
        if(cnt[p[i]] == 0) break;

        

    }


    

}

int main() { return 0; }
