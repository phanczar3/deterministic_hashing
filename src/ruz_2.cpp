#include <iostream>
#include <vector>
#include <math.h>
#include <map>
#include <chrono>

using namespace std;
using ll = long long;

struct label {
    int v, j, k;
    vector<int> active_bits;
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
    function<int(ll)> get_v = [&](ll x) -> int { return h(x) >> (psi - p[phase - 1]); };
    function<int(ll)> get_bits = [&](ll x) -> int { return (h(x) >> (psi - p[phase])) & ((1 << (p[phase] - p[phase - 1])) - 1); };

    
    // sort keys by (get_v, get_phi)
    counting_sort(keys, get_phi, psi);
    counting_sort(keys, get_v, psi);

    vector<label> list;
    int j = phi, k = 0, v = get_v(keys[0]), count = 1, i_start = 0;
    for(int i = 1; i < (int)keys.size(); i++) {
        if(v != get_v(keys[i])) {
            list.push_back({v, j, k, {}});
            for(int idx = i_start; idx < i; idx++) {
                list.back().active_bits.push_back(get_bits(keys[idx]));
            }
            j = phi, k = 0, v = get_v(keys[i]), count = 1, i_start = i;
        } else if(get_phi(keys[i]) >= ((k + 1) << j)) {
            list.push_back({v, j, k, {}});
            for(int idx = i_start; idx < i; idx++) {
                list.back().active_bits.push_back(get_bits(keys[idx]));
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
                        list.push_back({v, j, k, {}});
                        while(get_phi(keys[i_start]) < ((k + 1) << j)) {
                            list.back().active_bits.push_back(get_bits(keys[i_start]));
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
    list.push_back({v, j, k, {}});
    for(int idx = i_start; idx < (int)keys.size(); idx++) {
        list.back().active_bits.push_back(get_bits(keys[idx]));
    }
    return list;
}

vector<vector<level>> make_leaf_list(vector<label> &labels, int &start_labels, int target_j, int max_group) {
    int end_labels = start_labels;
    while(end_labels < (int)labels.size() && labels[end_labels].j == target_j) end_labels++;
    int no_labels = end_labels - start_labels;
    vector<vector<level>> leaf_list(max_group + 1, vector<level>(no_labels));
    for(int idx_labels = start_labels; idx_labels < end_labels; idx_labels++) {
        level lvl;
        lvl.k = labels[idx_labels].k, lvl.v = labels[idx_labels].v;
        lvl.w.resize(1 << max_group);
        for(int bits : labels[idx_labels].active_bits) {
            assert(bits >= 0 && bits < (1 << max_group));
            lvl.w[bits]++;
        }
        leaf_list[max_group][idx_labels - start_labels] = lvl;
    }

    for(int cur_group = max_group - 1; cur_group >= 1; cur_group--) {
        for(int idx_lvl = 0; idx_lvl < no_labels; idx_lvl++) {
            level lvl, &prev_lvl = leaf_list[cur_group + 1][idx_lvl];
            lvl.k = prev_lvl.k, lvl.v = prev_lvl.v;
            lvl.w.resize(1 << cur_group);
            for(int idx_w = 0; idx_w < (1 << cur_group); idx_w++) {
                lvl.w[idx_w] = prev_lvl.w[2 * idx_w] + prev_lvl.w[2 * idx_w + 1];
            }
            leaf_list[cur_group][idx_lvl] = lvl;
        }
    }
    
    start_labels = end_labels;
    return leaf_list;
}

vector<vector<level>> merge_lists(vector<vector<level>> &list1, vector<vector<level>> &list2) {
    vector<vector<level>> list(list1.size());

    for(int cur_group = 1; cur_group < (int)list1.size(); cur_group++) {
        int left_idx = 0, right_idx = 0;
        while(left_idx < (int)list1[cur_group].size() && right_idx < (int)list2[cur_group].size()) {
            level &left_lvl = list1[cur_group][left_idx], &right_lvl = list2[cur_group][right_idx];

            if(left_lvl.k < right_lvl.k) {
                list[cur_group].push_back(left_lvl);
                left_idx++;
            } else if(left_lvl.k > right_lvl.k) {
                list[cur_group].push_back(right_lvl);
                right_idx++;
            } else {
                if(left_lvl.v < right_lvl.v) {
                    list[cur_group].push_back(left_lvl);
                    left_idx++;
                } else if(left_lvl.v > right_lvl.v) {
                    list[cur_group].push_back(right_lvl);
                    right_idx++;
                } else {
                    assert(false);
                }
            }
        }

        while(left_idx < (int)list1[cur_group].size()) {
            list[cur_group].push_back(list1[cur_group][left_idx]);
            left_idx++;
        }

        while(right_idx < (int)list2[cur_group].size()) {
            list[cur_group].push_back(list2[cur_group][right_idx]);
            right_idx++;
        }
    }

    return list;
}

int get_end(vector<level> lvls, int idx_lvl, int target_k) {
    int idx2_lvl = idx_lvl;
    while(idx2_lvl < (int)lvls.size() && lvls[idx2_lvl].k == target_k) idx2_lvl++;
    return idx2_lvl;
}

void permute_vector(vector<int> &list, vector<int> &perm) {
    assert(list.size() == perm.size());
    // TODO in place (kind of)
    vector<int> new_list(list.size());
    for(int i = 0; i < (int)list.size(); i++) {
        new_list[perm[i]] = list[i];
    }
    list = move(new_list);
}  

pair<ll,ll> count_collisions(vector<level> &lvls, int left_start, int left_end, int right_start, int right_end) {
    ll coll0 = 0, coll1 = 0;
    while(left_start < left_end && right_start < right_end) {
        if(lvls[left_start].v < lvls[right_start].v) left_start++;
        else if(lvls[left_start].v > lvls[right_start].v) right_start++;
        else {
            for(int i = 0; i < (int)lvls[left_start].w.size(); i += 2) {
                coll0 += (ll) lvls[left_start].w[i] * lvls[right_start].w[i] + (ll) lvls[left_start].w[i | 1] * lvls[right_start].w[i | 1];
                coll1 += (ll) lvls[left_start].w[i] * lvls[right_start].w[i | 1] + (ll) lvls[left_start].w[i | 1] * lvls[right_start].w[i];
            }
            left_start++, right_start++;
        }
    }
    return {coll0, coll1};
}

const vector<int> phi_preprocessed = {4, 5, 5, 6, 7, 8, 8, 9, 10, 11, 12};

const map<pair<int,int>,pair<int,vector<int>>> p_preprocessed = {
    {{4,10}, {1, {0, 3, 5, 7, 10}}},
    {{5,11}, {1, {0, 3, 5, 7, 11}}},
    {{5,12}, {1, {0, 4, 6, 8, 12}}},
    {{6,13}, {1, {0, 4, 6, 8, 13}}},
    {{7,14}, {1, {0, 5, 7, 9, 14}}},
    {{8,15}, {1, {0, 4, 7, 10, 15}}},
    {{8,16}, {1, {0, 5, 8, 11, 16}}},
    {{9,17}, {1, {0, 5, 8, 11, 17}}},
    {{10,18}, {1, {0, 6, 9, 12, 18}}},
    {{11,19}, {2, {0, 6, 8, 11, 14, 17, 19}}},
    {{12,20}, {2, {0, 7, 9, 12, 15, 18, 20}}}
};


function<int(ll)> find_values(vector<ll> keys, int phi, int psi) {
    ll max_val = 0;
    for(int i = 0; i < (int)keys.size(); i++) {
        max_val = max(max_val, keys[i]);
    }
    int w = 64 - __builtin_clzll(max_val);
    assert(phi + psi >= w);

    int lphi = 31 - __builtin_clz(phi);

    vector<int> a(1 << phi, 0);

    function<int(ll)> h = [&a, psi](ll x) -> int {
        return (x & ((1LL << psi) - 1)) ^ a[x >> psi];
    };

    int istar;
    vector<int> p;

    auto it = p_preprocessed.find({phi, psi});
    if(it == p_preprocessed.end()) {
        assert(false);
    }
    istar = (*it).second.first, p = (*it).second.second;

    for(int phase = 1; phase <= 2 * istar + 2; phase++) {
        int eta = (1 << (p[phase] - p[phase - 1] + lphi));
        int max_group = p[phase] - p[phase - 1];

        vector<label> labels = labels_of_leafs(keys, h, p, phase, eta, phi, psi);
        // sort labels by (j, k, v)
        counting_sort(labels, [](label l) { return l.v; }, psi);
        counting_sort(labels, [](label l) { return l.k; }, psi);
        counting_sort(labels, [](label l) { return l.j; }, psi);

        int i_labels = 0;
        vector<vector<level>> list(max_group + 1);

        for(int j = 0; j < phi; j++) {
            vector<vector<level>> leaf_list = make_leaf_list(labels, i_labels, j, max_group);
            vector<vector<level>> merged_list = merge_lists(list, leaf_list);
            list = move(merged_list);
            
            int no_trees = list[1].size();
            for(int cur_group = 2; cur_group <= max_group; cur_group++) {
                assert((int)list[cur_group].size() == no_trees);
            }

            vector<int> delta(no_trees, 0);
            vector<vector<int>> perms(no_trees, vector<int>(1, 0));

            vector<vector<level>> new_list(max_group + 1);
            for(int cur_group = 1; cur_group <= max_group; cur_group++) {
                vector<level> &cur_list = list[cur_group]; 

                for(int lvl_start = 0; lvl_start < no_trees; ) {
                    int left_start = lvl_start, left_end = get_end(cur_list, left_start, (cur_list[lvl_start].k / 2) * 2);
                    int right_start = left_end, right_end = get_end(cur_list, right_start, (cur_list[lvl_start].k / 2) * 2 + 1);
                    
                    if(right_start == right_end) {
                        for(int idx_lvl = left_start; idx_lvl < left_end; idx_lvl++) {
                            cur_list[idx_lvl].k /= 2;
                            new_list[cur_group].push_back(cur_list[idx_lvl]);
                        }
                    } else if(left_start == left_end) {
                        for(int idx_lvl = right_start; idx_lvl < right_end; idx_lvl++) {
                            cur_list[idx_lvl].k /= 2;
                            new_list[cur_group].push_back(cur_list[idx_lvl]);
                        }   
                    } else {
                        // permute right children based on delta
                        vector<int> new_perm(perms[right_start].size() * 2);
                        assert(new_perm.size() == (1 << cur_group));
                        for(int i = 0; i < (int)new_perm.size(); i++) {
                            if(i % 2 == 0) new_perm[i] = 2 * perms[right_start][i / 2];
                            else new_perm[i] = 2 * perms[right_start][i / 2] + 1;
                        }

                        perms[right_start] = move(new_perm);
                        for(int idx_lvl = left_start; idx_lvl < left_end; idx_lvl++) 
                            permute_vector(cur_list[idx_lvl].w, perms[right_start]);

                        auto collisions = count_collisions(cur_list, left_start, left_end, right_start, right_end);
                        
                        int b = 0;

                        if(collisions.first > collisions.second) {
                            b = 1;
                            delta[right_start] = delta[right_start] | (1 << (max_group - cur_group));
                            for(int i = 0; i < (int)perms[right_start].size(); i += 2) {
                                swap(perms[right_start][i], perms[right_start][i + 1]);
                            }
                            for(int idx_lvl = right_start; idx_lvl < right_end; idx_lvl++) {
                                for(int i = 0; i < (int)cur_list[idx_lvl].w.size(); i += 2) {
                                    swap(cur_list[idx_lvl].w[i], cur_list[idx_lvl].w[i + 1]);
                                }
                            }
                        }

                        while(left_start < left_end && right_start < right_end) {
                            if(cur_list[left_start].v == cur_list[right_start].v) {
                                for(int i = 0; i < (1 << cur_group); i += 2) { 
                                    if(b == 0) {
                                        cur_list[left_start].w[i] += cur_list[right_start].w[i];
                                        cur_list[left_start].w[i | 1] += cur_list[right_start].w[i | 1];
                                    } else {
                                        cur_list[left_start].w[i] += cur_list[right_start].w[i | 1];
                                        cur_list[left_start].w[i | 1] += cur_list[right_start].w[i];
                                    }
                                }
                                cur_list[left_start].k /= 2;
                                new_list[cur_group].push_back(cur_list[left_start]);
                                left_start++, right_start++;
                            } else if(cur_list[left_start].v < cur_list[right_start].v) {
                                cur_list[left_start].k /= 2;
                                new_list[cur_group].push_back(cur_list[left_start]);
                                left_start++;
                            } else {
                                cur_list[right_start].k /= 2;
                                new_list[cur_group].push_back(cur_list[right_start]);
                                right_start++;
                            }
                        }

                        while(left_start < left_end) {
                            cur_list[left_start].k /= 2;
                            new_list[cur_group].push_back(cur_list[left_start]);
                            left_start++;
                        }

                        while(right_start < right_end) {
                            cur_list[right_start].k /= 2;
                            new_list[cur_group].push_back(cur_list[right_start]);
                            right_start++;
                        }
                    }
                    lvl_start = right_end;

                }
            }

            for(int i_tree = 0; i_tree < no_trees; i_tree++) {
                if(delta[i_tree] != 0) {
                    int cur_k = list[1][i_tree].k;
                    for(int i = cur_k << j; i < ((cur_k + 1) << j); i++) {
                        a[i] ^= (delta[i_tree] << p[phase - 1]);
                    }
                }
            }
            list = move(new_list);
        }
    }

    function<int(ll)> h_ret = [a, psi](ll x) -> int {
        return (x & ((1LL << psi) - 1)) ^ a[x >> psi];
    };

    return h_ret;
}

function<int(ll)> ruz_2(vector<ll> keys) {
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
        int w0 = 0;
        for(ll x : cur_keys) {
            w0 = max(w0, 64 - __builtin_clzll(x));
        }
        for(int j = 0; j < n; j++) {
            cur_keys[j] >>= cur_w - max_len;
        }
        seq[i] = find_values(cur_keys, phi_preprocessed[lN - 10], lN);
        for(int j = 0; j < n; j++) {
            keys[j] &= (1LL << (cur_w - max_len)) - 1;
            ll hash = (ll)seq[i](cur_keys[j]);
            assert(hash < (1LL << (psi)));
            hash <<= (cur_w - max_len);
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

    auto f = ruz_2(keys);

    // chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
    // chrono::duration<double> elapsed = chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    
    // cout << elapsed.count() << "\n";

    for(int i = 0; i < n; i++) {
        cout << keys[i] << " " << f(keys[i]) << "\n";
    }


    return 0;
}