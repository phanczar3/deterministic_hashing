#include <functional>
#include <iostream>
#include <vector>
#include <random>

using namespace std;
using ll = long long;
mt19937 rng(random_device{}());

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

    vector<pair<ll,ll>> v2(n);
    for(int i = 0; i < n; i++) {
        v2[--cnt[v[i].first]] = v[i];
    }
    v = move(v2);

    vector<int> p(1 << r);
    for(int i = n, p_idx = 0; i >= 0; i--) {
        for(auto val : vals[i]) {
            p[p_idx++] = val;
        }
    }

    vector<int> m(1 << r, 0), a(1 << r, 0);
    uniform_int_distribution<int> get_random(0, (1 << r) - 1);
    int cnt_prev = 0;
    for(int i = 0; i < (1 << r); i++) {
        int cur_val = p[i];
        int cnt_val = (cur_val + 1 == 1 << r ? n : cnt[cur_val + 1]) - cnt[cur_val];
        ll max_collisions = (1LL * cnt_val * cnt_prev) >> (r - 1), new_collisions;

        do {
            a[cur_val] = get_random(rng);
            new_collisions = 0;
            for(int j = cnt[cur_val]; j < cnt[cur_val] + cnt_val; j++) {
                new_collisions += m[v[j].second ^ a[cur_val]];
            }
        } while(new_collisions > max_collisions);
        
        for(int j = cnt[cur_val]; j < cnt[cur_val] + cnt_val; j++) {
            m[v[j].second ^ a[cur_val]]++;
        }
        
        cnt_prev += cnt_val;
    }

    function<ll(ll)> h = [a, f, g](ll x) -> ll {
        return g(x) ^ a[f(x)];
    };

    return {h, f};
}

int main() { 
    int n;
    cin >> n;


    vector<ll> keys(n);
    int logn = 31 - __builtin_clz(n), w = 0;
    for(int i = 0; i < n; i++) {
        cin >> keys[i];
        w = max(w, 64 - __builtin_clzll(keys[i]));
    }

    int r = max(w / 2, logn + 4);
    function<ll(ll)> f = [r, w](ll x) -> ll {
        return x >> max(0, w - r);
    };
    function<ll(ll)> g = [r](ll x) -> ll {
        return x & ((1 << r) - 1);
    };

    auto p2 = displace(keys, f, g, r);
    auto p3 = displace(keys, p2.first, p2.second, r);

    for(int i = 0; i < n; i++) {
        cout << keys[i] << " " << p3.first(keys[i]) << "\n";
    }

    return 0; 
}
