#include <functional>
#include <iostream>
#include <vector>

using namespace std;
using ll = long long;

void add_tree(vector<int> &t, int idx, int val) {
    int base = t.size() / 2, cur = base + idx;
    while(cur > 0) {
        t[cur] += val;
        cur /= 2;
    }
}

ll calc_ev(vector<int> &t, int u, int len, vector<pair<int,int>> &v, int st, int en) {
    int r = 31 - __builtin_clz(t.size() / 2);
    ll res = 0;
    for(int i = st; i < en; i++) {
        int seq = (v[i].second >> (r - len)) ^ u;
        res += t[(1 << len) + seq];
    }
    return res;
}

int get_val(vector<int> &t, vector<pair<int,int>> &v, int st, int en) {
    int base = t.size() / 2, cur = 1, u = 0, len = 0;
    while(cur < base) {
        if(calc_ev(t, u << 1, len + 1, v, st, en) < calc_ev(t, (u << 1) + 1, len + 1, v, st, en)) {
            cur = 2 * cur;
            u = u << 1;
        } else {
            cur = 2 * cur + 1;
            u = (u << 1) + 1;
        }
        len++;
    }
    return cur - base;
}

pair<function<ll(ll)>, function<ll(ll)>>
displace(vector<ll> keys, function<ll(ll)> f, function<ll(ll)> g, int r) {
	int n = keys.size();

	vector<pair<int, int>> v(n);
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
    vector<pair<int,int>> v2(n);
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

    vector<int> m(1 << (r + 1), 0), a(1 << r, 0);
    int cnt_prev = 0;
    for(int i = 0; i < (1 << r); i++) {
        int cur_val = p[i];
        int cnt_val = (cur_val + 1 == 1 << r ? n : cnt[cur_val + 1]) - cnt[cur_val];
        ll max_collisions = (1LL * cnt_val * cnt_prev) >> r, new_collisions = 0;

        a[cur_val] = get_val(m, v, cnt[cur_val], cnt[cur_val] + cnt_val);

        for(int j = cnt[cur_val]; j < cnt[cur_val] + cnt_val; j++) {
            new_collisions += m[(1 << r) + (v[j].second ^ a[cur_val])];
        }
        assert(new_collisions <= max_collisions);

        for(int j = cnt[cur_val]; j < cnt[cur_val] + cnt_val; j++) {
            add_tree(m, v[j].second ^ a[cur_val], 1);
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

    auto [f2, g2] = displace(keys, f, g, r);
    auto [f3, g3] = displace(keys, f2, g2, r);

    for(int i = 0; i < n; i++) {
        cout << keys[i] << " " << f3(keys[i]) << "\n";
    }

    return 0; 
}
