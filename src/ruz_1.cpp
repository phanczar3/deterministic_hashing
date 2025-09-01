#include <iostream>
#include <vector>
#include <chrono>
#include <assert.h>

using namespace std;
using ll = long long;

struct mst {
    int k, v, cnt0, cnt1;
};

template<typename T, typename F> 
void counting_sort(vector<T> &v, F f, int psi) {
    vector<int> cnt(1 << psi, 0);
    for(int i = 0; i < (int)v.size(); i++) {
        cnt[f(v[i])]++;
    }

    for(int i = 1; i < (1 << psi); i++) {
        cnt[i] += cnt[i - 1];
    }

    vector<T> v2(v.size());
    for(int i = (int)v.size() - 1; i >= 0; i--) {
        v2[--cnt[f(v[i])]] = v[i];
    }
    v = std::move(v2);
}

const vector<int> phi_preprocessed = {4, 5, 5, 6, 7, 8, 8, 9, 10, 11, 12};

function<int(ll)> find_values(vector<ll> keys, int phi, int psi) {
    
    ll max_val = 0;
    for(int i = 0; i < (int)keys.size(); i++) {
        max_val = max(max_val, keys[i]);
    }
    int w = 64 - __builtin_clzll(max_val);
    assert(phi + psi >= w);


    vector<int> a(1 << phi, 0);

    function<int(ll)> h = [&a, psi](ll x) -> int {
        return (x & ((1LL << psi) - 1)) ^ a[x >> psi];
    };

    for(int pos = 0; pos < psi; pos++) {
        auto get_v = [&](const ll &x) -> int {
            return h(x) >> (psi - pos);
        };

        auto get_phi = [&](const ll &x) -> int {
            return x >> psi;
        };

        auto get_bit = [&](const ll &x) -> int {
            return (x >> (psi - pos - 1)) & 1; 
        };

        counting_sort(keys, get_v, psi);
        counting_sort(keys, get_phi, psi);
        
        // [k, v, cnt0, cnt1]
        vector<mst> list(1);
        list[0] = {get_v(keys[0]), get_phi(keys[0]), 0, 0};
        if(get_bit(keys[0]) == 0) list[0].cnt0++;
        else list[0].cnt1++;
        
        for(int i = 1; i < (int)keys.size(); i++) {
            mst &ref = list[list.size() - 1];
            int pk = ref.k, pv = ref.v;
            int ck = get_phi(keys[i]), cv = get_v(keys[i]), b = get_bit(keys[i]);
            if(pk == ck && pv == cv) {
                if(b == 0) list[list.size() - 1].cnt0++;
                else list[list.size() - 1].cnt1++;
            }
            if(b == 0) list.push_back({ck, cv, 1, 0});
            else list.push_back({ck, cv, 0, 1});
        }

        for(int j = 0; j < phi; j++) {

            vector<mst> new_list;
            int left_start = 0, left_end = 0, right_start = 0, right_end = 0;
            while(left_start < (int)list.size()) {
                if(list[left_start].k % 2 == 0) {
                    while(left_end < (int)list.size() && list[left_start].k == list[left_end].k) left_end++;
                    right_start = right_end = left_end;
                }

                if(right_start < (int)list.size() &&
                ((right_end - right_start > 0 && list[right_start].k + 1 == list[left_start].k) || right_end == right_start)) {
                    while(right_end < (int)list.size() && list[right_start].k == list[right_end].k) right_end++;
                }
                // [left_start, left_end) - left child multisets, [right_start, right_end) - right child multisets

                if(right_start == right_end) {
                    for(int i = left_start; i < left_end; i++) {
                        list[i].k /= 2;
                        new_list.push_back(list[i]);
                    }
                } else if(left_start == left_end) {
                    for(int i = right_start; i < right_end; i++) {
                        list[i].k /= 2;
                        new_list.push_back(list[i]);
                    }
                } else {
                    // count number of collisions between multisets
                    ll coll0 = 0, coll1 = 0;
                    int left_ptr = left_start, right_ptr = right_start;
                    while(left_ptr < left_end && right_ptr < right_end) {
                        mst &ref1 = list[left_ptr], &ref2 = list[right_ptr];
                        int lv = ref1.v, lcnt0 = ref1.cnt0, lcnt1 = ref1.cnt1;
                        int rv = ref2.v, rcnt0 = ref2.cnt0, rcnt1 = ref2.cnt1;

                        if(lv == rv) {
                            coll0 += (ll)lcnt0 * rcnt0 + (ll)lcnt1 * rcnt1;
                            coll1 += (ll)lcnt0 * rcnt1 + (ll)lcnt1 * rcnt0;
                            left_ptr++, right_ptr++;
                        } else if(lv < rv) {
                            left_ptr++;
                        } else {
                            right_ptr++;
                        }
                    }
                    int b = (coll0 <= coll1 ? 0 : 1);
                    left_ptr = left_start, right_ptr = right_start;

                    // new multiset based on two children
                    while(left_ptr < left_end && right_ptr < right_end) {
                        mst &ref1 = list[left_ptr], &ref2 = list[right_ptr];
                        int lk = ref1.k, lv = ref1.v, lcnt0 = ref1.cnt0, lcnt1 = ref1.cnt1;
                        int rv = ref2.v, rcnt0 = ref2.cnt0, rcnt1 = ref2.cnt1;

                        if(lv == rv) {
                            list[left_ptr].k /= 2, list[right_ptr].k /= 2;
                            new_list.push_back({lk / 2, lv, 
                                lcnt0 + (b == 0 ? rcnt0 : rcnt1), 
                                lcnt1 + (b == 0 ? rcnt1 : rcnt0)});
                            left_ptr++, right_ptr++;
                        } else if(lv < rv) {
                            list[left_ptr].k /= 2;
                            new_list.push_back(list[left_ptr++]);
                        } else {
                            list[right_ptr].k /= 2;
                            new_list.push_back(list[right_ptr++]);
                        }
                    }

                    // update a_r accordingly
                    if(b == 1) {
                        int cur_k = list[right_start].k * 2 + 1;
                        for(int i = cur_k << j; i < ((cur_k + 1) << j); i++) {
                            assert(i < (1 << phi));
                            a[i] ^= (1 << (psi - pos - 1));
                        }
                    }
                }
                
                // move pointers to next segments
                left_start = left_end = right_start = right_end = right_end + 1;
            }

            list = std::move(new_list);
        }
    }

    function<int(ll)> h_ret = [a, psi](ll x) -> int {
        return (x & ((1LL << psi) - 1)) ^ a[x >> psi];
    };

    return h_ret;
}

function<int(ll)> ruz_1(vector<ll> keys) {
    assert((int)keys.size() >= 513 && (int)keys.size() <= 1048576);

    int n = keys.size();
    ll max_val = 0;
    
    for(int i = 0; i < n; i++) {
        // only positive keys
        assert(keys[i] >= 0);
        max_val = max(max_val, keys[i]);
    }

    int N =  1 << (31 - __builtin_clz(n - 1) + 1);
    int lN = 31 - __builtin_clz(N);
    int psi = lN, phi = phi_preprocessed[lN - 10];
    int max_len = psi + phi;
    int delta = phi, w = 64 - __builtin_clzll(max_val);
    int times = 1 + max(0, (w - max_len + delta - 1) / delta);
    int cur_w = w;

    vector<function<int(ll)>> seq(times);
    for(int i = 0; i < times - 1; i++) {
        vector<ll> cur_keys = keys;
        for(int j = 0; j < n; j++) {
            cur_keys[j] >>= cur_w - max_len;
        }
        seq[i] = find_values(cur_keys, phi_preprocessed[lN - 10], lN);
        for(int j = 0; j < n; j++) {
            keys[j] &= (1LL << (cur_w - max_len)) - 1;
            ll hash = (ll)seq[i](cur_keys[j]) << (cur_w - max_len);
            assert((hash & keys[j]) == 0);
            keys[j] |= hash;
        }
        cur_w -= delta;
    }
    seq[times - 1] = find_values(keys, phi_preprocessed[lN - 10], lN);
    

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

    auto f = ruz_1(keys);

    // chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
    // chrono::duration<double> elapsed = chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);

    // cout << elapsed.count() << "\n";

    for(int i = 0; i < n; i++) {
        cout << keys[i] << " " << f(keys[i]) << "\n";
    }


    return 0;
}