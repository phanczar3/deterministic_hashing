#include <iostream>
#include <vector>
#include <math.h>

using namespace std;
using ll = long long;

struct label {
    int v, j, k;
    vector<ll> l_keys;
};

struct level {
    int k, v;
    vector<int> w;
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
    v = move(v2);
}

vector<label> labels_of_leafs(vector<ll> keys, const function<int(ll)> &h, const vector<int> &p, int phase, int eta, int phi, int psi) {
    function<int(ll)> get_phi = [&](ll x) -> int { return x >> psi; };
    function<int(ll)> get_v = [&](ll x) -> int { return h(x) >> p[phase - 1]; };
    
    // sort keys by (get_v, get_phi)
    counting_sort(keys, get_phi, psi);
    counting_sort(keys, get_v, psi);

    vector<label> list;
    int j = phi, k = 0, v = get_v(keys[0]), count = 1, i_start = 0;
    for(int i = 1; i < (int)keys.size(); i++) {
        if(v != get_v(keys[i])) {
            list.push_back({v, j, k});
            for(int idx = i_start; idx < i; idx++) {
                list.back().l_keys.push_back(keys[idx]);
            }
            j = phi, k = 0, v = get_v(keys[i]), count = 1, i_start = i;
        } else if(get_phi(keys[i]) >= ((k + 1) << j)) {
            list.push_back({v, j, k});
            for(int idx = i_start; idx < i; idx++) {
                list.back().l_keys.push_back(keys[idx]);
            }
            while( ((k + 2) << j) < get_phi(keys[i])) {
                j++, k /= 2;
            }
            k++, count = 1, i_start = i;
        } else {
            count++;
            while(count >= eta && j > 0) {
                while(j > 0 && ((k + 1) << j) > get_phi(keys[i])) {
                    j--, k *= 2;
                }
                if(((k + 1) << j) <= get_phi(keys[i])) {
                    if(get_phi(keys[i_start]) < ((k + 1) << j)) {
                        list.push_back({v, j, k});
                        while(get_phi(keys[i_start]) < ((k + 1) << j)) {
                            list.back().l_keys.push_back(keys[i_start]);
                            i_start++, count--;
                        }
                    }
                    k++;
                }
            }
        }
    }
    // last label
    assert(i_start < (int)keys.size());
    list.push_back({v, j, k});
    for(int idx = i_start; idx < (int)keys.size(); idx++) {
        list.back().l_keys.push_back(keys[idx]);
    }
    return list;
}

vector<vector<level>> make_leaf_list(vector<label> &labels, int &i_labels, function<int(ll)> &h, int j) {
    int m = labels.size();
    int j_labels = i_labels;
    while(j_labels < m && labels[j_labels].j == j) j_labels++;
    
    for(int idx = i_labels; idx < j_labels; idx++) {
        
    }



    i_labels = j_labels;
}

function<int(ll)> find_values(vector<ll> keys, int phi, int psi) {
    
    assert(phi >= 2);

    ll max_val = 0;
    for(int i = 0; i < (int)keys.size(); i++) {
        max_val = max(max_val, keys[i]);
    }
    int w = 64 - __builtin_clzll(max_val);
    assert(phi + psi >= w);

    int lphi = 31 - __builtin_clz(phi), lpsi = 31 - __builtin_clz(psi);


    vector<int> a(1 << phi, 0);

    function<int(ll)> h = [&a, psi](ll x) -> int {
        return (x & ((1LL << psi) - 1)) ^ a[x >> psi];
    };

    int istar = floor(log2(psi) - log2(log2(phi)) - 1);


    vector<int> p(2 * istar + 3);
    p[0] = 0;
    for(int i = 1; i <= istar; i++) {
        p[i] = ceil(((double)(1LL << i) - 1) * psi / (1LL << i)) - i * lphi;
    }
    for(int i = istar + 1; i <= 2 * istar + 2; i++) {
        p[i] = p[i - 1] + lphi;
    }
    // for(int i = 0; i < 2 * istar + 3; i++) {
    //     cout << p[i] << " ";
    // }

    // check that sequence is correct
    for(int i = 0; i < 2 * istar + 2; i++) {
        assert(p[i] < p[i + 1]);
    }
    assert(p[2 * istar + 2] == psi);

    for(int phase = 1; phase <= 2 * istar + 2; phase++) {

        int eta = (1 << (p[phase] - p[phase - 1] + lphi));

        vector<label> labels = labels_of_leafs(keys, h, p, phase, eta, phi, psi);
        // sort labels by (j, k, v)
        counting_sort(labels, [](label l) { return l.v; }, psi);
        counting_sort(labels, [](label l) { return l.k; }, psi);
        counting_sort(labels, [](label l) { return l.j; }, psi);

        int i_labels = 0;
        for(int j = 0; j < phi; j++) {

            vector<vector<level>> leaf_list = make_leaf_list(labels, i_labels, h, j);


        }
        // for(int i = 0; i < (int)labels.size(); i++) {
        //     cout << labels[i].v << " " << labels[i].j << " " << labels[i].k << "\n";
        // }

        
        // list[0] = {get_v(keys[0]), get_phi(keys[0]), 0, 0};
        // if(get_bit(keys[0]) == 0) list[0].cnt0++;
        // else list[0].cnt1++;
        
        // for(int i = 1; i < (int)keys.size(); i++) {
        //     mst &ref = list[list.size() - 1];
        //     int pk = ref.k, pv = ref.v;
        //     int ck = get_phi(keys[i]), cv = get_v(keys[i]), b = get_bit(keys[i]);
        //     if(pk == ck && pv == cv) {
        //         if(b == 0) list[list.size() - 1].cnt0++;
        //         else list[list.size() - 1].cnt1++;
        //     }
        //     if(b == 0) list.push_back({ck, cv, 1, 0});
        //     else list.push_back({ck, cv, 0, 1});
        // }

        // for(int j = 0; j < phi; j++) {

        //     vector<mst> new_list;
        //     int left_start = 0, left_end = 0, right_start = 0, right_end = 0;
        //     while(left_start < (int)list.size()) {
        //         if(list[left_start].k % 2 == 0) {
        //             while(left_end < (int)list.size() && list[left_start].k == list[left_end].k) left_end++;
        //             right_start = right_end = left_end;
        //         }

        //         if(right_start < (int)list.size() &&
        //         ((right_end - right_start > 0 && list[right_start].k + 1 == list[left_start].k) || right_end == right_start)) {
        //             while(right_end < (int)list.size() && list[right_start].k == list[right_end].k) right_end++;
        //         }
        //         // [left_start, left_end) - left child multisets, [right_start, right_end) - right child multisets

        //         if(right_start == right_end) {
        //             for(int i = left_start; i < left_end; i++) {
        //                 list[i].k /= 2;
        //                 new_list.push_back(list[i]);
        //             }
        //         } else if(left_start == left_end) {
        //             for(int i = right_start; i < right_end; i++) {
        //                 list[i].k /= 2;
        //                 new_list.push_back(list[i]);
        //             }
        //         } else {
        //             // count number of collisions between multisets
        //             ll coll0 = 0, coll1 = 0;
        //             int left_ptr = left_start, right_ptr = right_start;
        //             while(left_ptr < left_end && right_ptr < right_end) {
        //                 mst &ref1 = list[left_ptr], &ref2 = list[right_ptr];
        //                 int lv = ref1.v, lcnt0 = ref1.cnt0, lcnt1 = ref1.cnt1;
        //                 int rv = ref2.v, rcnt0 = ref2.cnt0, rcnt1 = ref2.cnt1;

        //                 if(lv == rv) {
        //                     coll0 += (ll)lcnt0 * rcnt0 + (ll)lcnt1 * rcnt1;
        //                     coll1 += (ll)lcnt0 * rcnt1 + (ll)lcnt1 * rcnt0;
        //                     left_ptr++, right_ptr++;
        //                 } else if(lv < rv) {
        //                     left_ptr++;
        //                 } else {
        //                     right_ptr++;
        //                 }
        //             }
        //             int b = (coll0 <= coll1 ? 0 : 1);
        //             left_ptr = left_start, right_ptr = right_start;

        //             // new multiset based on two children
        //             while(left_ptr < left_end && right_ptr < right_end) {
        //                 mst &ref1 = list[left_ptr], &ref2 = list[right_ptr];
        //                 int lk = ref1.k, lv = ref1.v, lcnt0 = ref1.cnt0, lcnt1 = ref1.cnt1;
        //                 int rv = ref2.v, rcnt0 = ref2.cnt0, rcnt1 = ref2.cnt1;

        //                 if(lv == rv) {
        //                     new_list.push_back({lk / 2, lv, 
        //                         lcnt0 + (b == 0 ? rcnt0 : rcnt1), 
        //                         lcnt1 + (b == 0 ? rcnt1 : rcnt0)});
        //                     left_ptr++, right_ptr++;
        //                 } else if(lv < rv) {
        //                     list[left_ptr].k /= 2;
        //                     new_list.push_back(list[left_ptr++]);
        //                 } else {
        //                     list[right_ptr].k /= 2;
        //                     new_list.push_back(list[right_ptr++]);
        //                 }
        //             }

        //             // update a_r accordingly
        //             if(b == 1) {
        //                 for(int i = list[right_start].k << j; i < ((list[right_start].k + 1) << j); i++) {
        //                     assert(i < (1 << phi));
        //                     a[i] ^= (1 << (psi - pos - 1));
        //                 }
        //             }
        //         }
                
        //         // move pointers to next segments
        //         left_start = left_end = right_start = right_end = right_end + 1;
        //     }

        //     list = move(new_list);
        // }
    }

    function<int(ll)> h_ret = [a, psi](ll x) -> int {
        return (x & ((1LL << psi) - 1)) ^ a[x >> psi];
    };

    return h_ret;
}

function<int(ll)> ruz_1(vector<ll> keys) {
    int n = keys.size();
    ll max_val = 0;
    
    for(int i = 0; i < n; i++) {
        // only positive keys
        assert(keys[i] >= 0);
        max_val = max(max_val, keys[i]);
    }

    int N =  1 << (31 - __builtin_clz(n - 1) + 1);
    int lN = 31 - __builtin_clz(N), llN = 31 - __builtin_clz(lN - 1) + 1;
    int max_len = 2 * lN - 2 * llN;
    int delta = max_len - lN, w = 64 - __builtin_clzll(max_val);
    // assert n big enough
    assert(delta > 0);
    int times = 1 + max(0, (w - max_len + delta - 1) / delta);
    int cur_w = w;

    vector<function<int(ll)>> seq(times);
    for(int i = 0; i < times - 1; i++) {
        vector<ll> cur_keys = keys;
        for(int j = 0; j < n; j++) {
            cur_keys[j] >>= cur_w - max_len;
        }
        seq[i] = find_values(cur_keys, lN - 2 * llN, lN);
        for(int j = 0; j < n; j++) {
            keys[j] &= (1LL << (cur_w - max_len)) - 1;
            ll hash = (ll)seq[i](cur_keys[j]) << (cur_w - max_len);
            assert((hash & keys[j]) == 0);
            keys[j] |= hash;
        }
        cur_w -= delta;
    }
    seq[times - 1] = find_values(keys, lN - 2 * llN, lN);
    

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

    auto f = ruz_1(keys);

    for(int i = 0; i < n; i++) {
        cout << keys[i] << " " << f(keys[i]) << "\n";
    }


    return 0;
}