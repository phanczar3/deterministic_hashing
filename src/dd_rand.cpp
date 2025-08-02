#include <functional>
#include <iostream>
#include <vector>
#include <random>
#include <chrono>

using namespace std;
using ll = long long;
mt19937 rng(random_device{}());

pair<function<int(ll)>, function<int(ll)>>
displace(vector<ll> keys, function<int(ll)> f, function<int(ll)> g, int r) {
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

    function<int(ll)> h = [a, f, g](ll x) -> ll {
        return g(x) ^ a[f(x)];
    };

    return {h, f};
}

function<int(ll)> dd_rand(vector<ll> keys) {
    int n = keys.size(), ln = 31 - __builtin_clz(n), w = 0;

    for(int i = 0; i < n; i++) {
        w = max(w, 64 - __builtin_clzll(keys[i]));
    }

    int max_len = 2 * ln + 8;
    int delta = ln + 4;
    int times = 1 + max(0, (w - 2 * ln - 8 + delta - 1) / delta);
    int cur_w = w;

    vector<function<int(ll)>> seq(times);
    for(int i = 0; i < times; i++) {
        vector<ll> cur_keys = keys;
        if(i < times - 1) {
            for(int j = 0; j < n; j++) {
                cur_keys[j] >>= cur_w - max_len;
            }
        }

        function<int(ll)> f = [ln](ll x) -> ll {
            return x >> (ln + 4);
        };
        function<int(ll)> g = [ln](ll x) -> ll {
            return x & ((1 << (ln + 4)) - 1);
        };

        auto p2 = displace(cur_keys, f, g, ln + 4);
        auto p3 = displace(cur_keys, p2.first, p2.second, ln + 4);
        seq[i] = p3.first;
        for(int j = 0; j < n && i < times - 1; j++) {
            keys[j] &= (1LL << (cur_w - max_len)) - 1;
            ll hash = (ll)seq[i](cur_keys[j]) << (cur_w - max_len);
            assert((hash & keys[j]) == 0);
            keys[j] |= hash;
        }
        cur_w -= delta;
    }

    function<int(ll)> h = [seq, times, delta, w, max_len](ll x) -> int {
        int cur_w = w;
        for(int i = 0; i < times - 1; i++) {
            ll y = x;
            y >>= (cur_w - max_len);
            ll hash = (ll)seq[i](y) << (cur_w - max_len);
            x &= (1LL << (cur_w - max_len)) - 1;
            assert((hash & x) == 0);
            x |= hash;
            cur_w -= delta;
        }
        return seq[times - 1](x);
    };
    
    return h;
}

int main() { 
    int n;
    cin >> n;

    vector<ll> keys(n);

    for(int i = 0; i < n; i++) {
        cin >> keys[i];
    }

    // for benchmarks
    // chrono::steady_clock::time_point t0 = chrono::steady_clock::now();

    function<int(ll)> f = dd_rand(keys);

    // chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
    // chrono::duration<double> elapsed = chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);

    // cout << elapsed.count() << "\n";

    // for(int i = 0; i < n; i++) {
    //     cout << keys[i] << " " << f(keys[i]) << "\n";
    // }

    return 0; 
}
