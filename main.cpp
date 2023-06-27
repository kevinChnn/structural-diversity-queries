#include <algorithm>
#include <boost/functional/hash.hpp>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <deque>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;

int32_t g_n = 0;
int32_t g_m = 0;
int32_t g_t = 0;
int32_t g_t_interval = 1;
int32_t g_tau = 0;
vector<pair<int32_t, int32_t>> g_edges;
vector<int32_t> g_offset;

double g_sliding = 0.6; // default sliding window size 60% g_t

struct triplet
{
    int32_t first;
    int32_t second;
    int32_t third;
};

struct my_pair
{
    int32_t first;
    int32_t second;
};

struct ds_item
{
    int32_t parent;
    int32_t size;
};

// disjoint-set
int32_t ds_find(int32_t x, unordered_map<int32_t, ds_item> &ds)
{
    while (ds[x].parent != x)
    {
        ds[x].parent = ds[ds[x].parent].parent;
        x = ds[x].parent;
    }
    return x;
}

void ds_union(int32_t x, int32_t y, unordered_map<int32_t, ds_item> &ds)
{
    x = ds_find(x, ds);
    y = ds_find(y, ds);

    if (x == y)
        return;

    if (ds[x].size < ds[y].size)
        swap(x, y);

    ds[y].parent = x;
    ds[x].size = ds[x].size + ds[y].size;
}

// load graph
void load_graph(char const *file_path)
{
    FILE *f;
    if ((f = fopen(file_path, "r")) == NULL)
    {
        printf("Fail to open file!\n");
        exit(0);
    }

    fscanf(f, "%d %d %d", &g_n, &g_m, &g_t);

    int32_t u, v, t;
    int32_t base_t = -1;
    int32_t prev_t = -1;
    for (int32_t i = 0; i < g_m; ++i)
    {
        fscanf(f, "%d %d %d", &u, &v, &t);

        // select the earliest timestamp as the base time
        if (base_t == -1)
            base_t = t / g_t_interval;
        t = t / g_t_interval - base_t;

        if (t != prev_t)
        {
            g_offset.emplace_back(g_edges.size());
            prev_t = t;
        }

        // store edges
        g_edges.emplace_back(make_pair(u, v));
    }
    g_offset.emplace_back(g_edges.size());
    fclose(f);

    g_t = g_t / g_t_interval - base_t + 1;
}

// Arbitrary Window Query
// Sliding Window Query
// Online
int64_t online_query_one(int32_t t_s, int32_t t_e, int32_t vertex)
{
    auto t1 = chrono::high_resolution_clock::now();

    int32_t sd = 0;

    // handle duplicate edges
    vector<unordered_set<int32_t>> nbr_set(g_n);

    // map old vertex id to new vertex id
    vector<int32_t> new_to_old;
    unordered_map<int32_t, int32_t> old_to_new;

    // adjacent lists
    vector<vector<int32_t>> adj;

    // prevent g_offset[t_e + 1] out of range
    if (t_e >= g_t)
        t_e = g_t - 1;

    for (int32_t i = g_offset[t_s]; i < g_offset[t_e + 1]; ++i)
    {
        int32_t u = g_edges[i].first;
        int32_t v = g_edges[i].second;

        // skip duplicate edges
        if (u > v)
            swap(u, v);
        // check if edge (u,v) is duplicate
        if (nbr_set[u].find(v) != nbr_set[u].end())
            continue;
        // record edge (u,v) for duplicate check
        nbr_set[u].emplace(v);
        // nbr_set[v].emplace(u);

        if (old_to_new.find(u) == old_to_new.end())
        {
            old_to_new.emplace(u, (int32_t)new_to_old.size());
            new_to_old.emplace_back(u);
            adj.emplace_back(vector<int32_t>());
        }
        u = old_to_new[u];

        if (old_to_new.find(v) == old_to_new.end())
        {
            old_to_new.emplace(v, (int32_t)new_to_old.size());
            new_to_old.emplace_back(v);
            adj.emplace_back(vector<int32_t>());
        }
        v = old_to_new[v];

        // add edge (u,v) to adjacent lists
        adj[u].emplace_back(v);
        adj[v].emplace_back(u);
    }

    if (!old_to_new.count(vertex))
        return sd;

    auto t2 = chrono::high_resolution_clock::now();

    // initialize disjoint-set for structural diversity computation
    unordered_map<int32_t, ds_item> ds;
    for (auto &v : adj[old_to_new[vertex]])
    {
        ds.emplace(v, (ds_item){v, 1});
    }

    if (g_tau == 1)
    {
        sd = (int32_t)adj[old_to_new[vertex]].size();
        return sd;
    }

    bool *visited = new bool[new_to_old.size()]();
    // list triangles
    for (int i = 0; i < adj[old_to_new[vertex]].size(); ++i)
    {
        visited[adj[old_to_new[vertex]][i]] = true;
    }

    for (int i = 0; i < adj[old_to_new[vertex]].size(); ++i)
    {
        int v = adj[old_to_new[vertex]][i];
        for (int j = 0; j < adj[v].size(); ++j)
        {
            int w = adj[v][j];
            if (adj[w].size() < adj[v].size() || (adj[w].size() == adj[v].size() && w < v))
            {
                if (visited[w])
                {
                    int32_t sa = ds_find(v, ds);
                    int32_t sb = ds_find(w, ds);
                    if (sa != sb)
                    {
                        if (max(ds[sa].size, ds[sb].size) < g_tau && ds[sa].size + ds[sb].size >= g_tau)
                        {
                            sd += 1;
                        }
                        else if (ds[sa].size >= g_tau && ds[sb].size >= g_tau)
                        {
                            sd -= 1;
                        }
                        ds_union(v, w, ds);
                    }
                }
            }
        }
    }

    delete[] visited;
    // return sd;
    return chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
}

void check_inf(int32_t u, int32_t v, int32_t w, bool *inf, vector<int32_t> &inf_vertices)
{
    if (!inf[u])
    {
        inf[u] = true;
        inf_vertices.emplace_back(u);
    }
    if (!inf[v])
    {
        inf[v] = true;
        inf_vertices.emplace_back(v);
    }
    if (!inf[w])
    {
        inf[w] = true;
        inf_vertices.emplace_back(w);
    }
}

// Arbitrary Window Query
// Baseline
int32_t historical_bl_t_insert(int32_t u, int32_t v, int32_t w, vector<int32_t> &sd, unordered_map<int32_t, ds_item> &ds)
{
    int32_t sa = ds_find(v, ds);
    int32_t sb = ds_find(w, ds);
    if (sa != sb)
    {
        if (max(ds[sa].size, ds[sb].size) < g_tau && ds[sa].size + ds[sb].size >= g_tau)
        {
            sd[u] += 1;
        }
        else if (ds[sa].size >= g_tau && ds[sb].size >= g_tau)
        {
            sd[u] -= 1;
        }
        ds_union(v, w, ds);
    }
    return sd[u];
}

// Arbitrary Window Query
// Baseline
void historical_bl_t_insert_index(int32_t u, int32_t v, int32_t w, vector<int32_t> &sd, vector<vector<triplet>> &idx, unordered_map<int32_t, ds_item> &ds, int32_t t_s, int32_t t_e)
{
    int32_t new_value = historical_bl_t_insert(u, v, w, sd, ds);

    while (idx[u].back().first == t_s && idx[u].back().second == t_e)
    {
        idx[u].pop_back();
    }
    if (new_value == idx[u].back().third)
        return;
    idx[u].emplace_back((triplet){t_s, t_e, new_value});
}

// Arbitrary Window Query
// Baseline
int64_t historical_bl_index(vector<vector<triplet>> &idx)
{
    int64_t init_time = 0;

    auto t3 = chrono::high_resolution_clock::now();
    auto t4 = chrono::high_resolution_clock::now();

    int64_t idx_init_time = 0;

    t3 = chrono::high_resolution_clock::now();
    idx.resize(g_n);

    bool *visited = new bool[g_n]();
    bool *triangle_a = new bool[g_n]();
    bool *triangle_b = new bool[g_n]();
    bool *triangle_c = new bool[g_n]();
    t4 = chrono::high_resolution_clock::now();
    idx_init_time += chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();

    for (int32_t t_e = 0; t_e < g_t; ++t_e)
    {
        auto t1 = chrono::high_resolution_clock::now();

        // remove duplicate edges
        vector<unordered_set<int32_t>> nbr_set(g_n);

        vector<vector<int32_t>> adj(g_n);
        vector<int32_t> adj_offset(g_n, 0);

        vector<unordered_map<int32_t, ds_item>> adj_ds(g_n);
        vector<int32_t> sd(g_n, 0);

        for (int32_t u = 0; u < g_n; u++)
        {
            idx[u].emplace_back((triplet){t_e, t_e, 0});
        }

        for (int32_t t_s = t_e; t_s >= 0; --t_s)
        {
            vector<int32_t> inf_vertices;
            // add edges at t_s
            for (int32_t i = g_offset[t_s]; i < g_offset[t_s + 1]; ++i)
            {
                int32_t u = g_edges[i].first;
                int32_t v = g_edges[i].second;

                // skip duplicate edges
                if (u > v)
                    swap(u, v);
                // check if edge exists
                if (nbr_set[u].find(v) != nbr_set[u].end())
                    continue;
                // add new edges
                nbr_set[u].insert(v);
                // nbr_set[v].insert(u);

                adj[u].emplace_back(v);
                adj[v].emplace_back(u);

                ds_item dsv{v, 1};
                adj_ds[u].insert(make_pair(v, dsv));
                if (g_tau == 1)
                    sd[u]++;

                ds_item dsu{u, 1};
                adj_ds[v].insert(make_pair(u, dsu));
                if (g_tau == 1)
                    sd[v]++;

                if (!visited[u])
                {
                    visited[u] = true;
                    inf_vertices.emplace_back(u);
                }

                if (!visited[v])
                {
                    visited[v] = true;
                    inf_vertices.emplace_back(v);
                }
            }

            // compute new triangles
            for (auto &u : inf_vertices)
            {
                vector<int32_t> ta, tb, tc;
                for (int32_t i = 0; i < adj_offset[u]; ++i)
                {
                    int32_t v = adj[u][i];
                    if (adj[v].size() < adj[u].size() || (adj[v].size() == adj[u].size() && v < u))
                    {
                        triangle_a[v] = true;
                        ta.emplace_back(v);
                    }
                    else
                    {
                        triangle_b[v] = true;
                        tb.emplace_back(v);
                    }
                }
                for (int32_t i = adj_offset[u]; i < adj[u].size(); ++i)
                {
                    int32_t v = adj[u][i];
                    if (adj[v].size() < adj[u].size() || (adj[v].size() == adj[u].size() && v < u))
                    {
                        triangle_c[v] = true;
                        tc.emplace_back(v);
                    }
                }

                for (auto &v : tc)
                {
                    for (int32_t i = 0; i < adj_offset[v]; ++i)
                    {
                        int32_t w = adj[v][i];
                        // CASE 1
                        if (triangle_a[w] || triangle_b[w])
                        {
                            // find triangle (u,v,w)
                            // printf("<%d,%d,%d>\n",u,v,w);
                            historical_bl_t_insert_index(u, v, w, sd, idx, adj_ds[u], t_s, t_e);
                            historical_bl_t_insert_index(v, u, w, sd, idx, adj_ds[v], t_s, t_e);
                            historical_bl_t_insert_index(w, u, v, sd, idx, adj_ds[w], t_s, t_e);
                        }

                        // CASE 2.1
                        if (adj[w].size() < adj[v].size() || (adj[w].size() == adj[v].size() && w < v))
                        {
                            if (triangle_c[w])
                            {
                                // find trinalge (u,v,w)
                                // printf("<%d,%d,%d>\n",u,v,w);
                                historical_bl_t_insert_index(u, v, w, sd, idx, adj_ds[u], t_s, t_e);
                                historical_bl_t_insert_index(v, u, w, sd, idx, adj_ds[v], t_s, t_e);
                                historical_bl_t_insert_index(w, u, v, sd, idx, adj_ds[w], t_s, t_e);
                            }
                        }
                    }

                    for (int32_t i = adj_offset[v]; i < adj[v].size(); ++i)
                    {
                        int32_t w = adj[v][i];
                        // CASE 2.2
                        if (triangle_a[w])
                        {
                            // find triangle (u,v,w)
                            // printf("<%d,%d,%d>\n",u,v,w);
                            historical_bl_t_insert_index(u, v, w, sd, idx, adj_ds[u], t_s, t_e);
                            historical_bl_t_insert_index(v, u, w, sd, idx, adj_ds[v], t_s, t_e);
                            historical_bl_t_insert_index(w, u, v, sd, idx, adj_ds[w], t_s, t_e);
                        }

                        // CASE 3
                        if (adj[w].size() < adj[v].size() || (adj[w].size() == adj[v].size() && w < v))
                        {
                            if (triangle_c[w])
                            {
                                // find trinalge (u,v,w)
                                // printf("<%d,%d,%d>\n",u,v,w);
                                historical_bl_t_insert_index(u, v, w, sd, idx, adj_ds[u], t_s, t_e);
                                historical_bl_t_insert_index(v, u, w, sd, idx, adj_ds[v], t_s, t_e);
                                historical_bl_t_insert_index(w, u, v, sd, idx, adj_ds[w], t_s, t_e);
                            }
                        }
                    }
                }

                for (auto &v : ta)
                    triangle_a[v] = false;
                for (auto &v : tb)
                    triangle_b[v] = false;
                for (auto &v : tc)
                    triangle_c[v] = false;
            }

            for (auto &u : inf_vertices)
                adj_offset[u] = adj[u].size();
            for (auto &u : inf_vertices)
                visited[u] = false;
        }

        auto t2 = chrono::high_resolution_clock::now();
        if (t_e < (int32_t)(g_t * g_sliding))
        {
            init_time += chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        }
    }

    delete[] visited;
    delete[] triangle_a;
    delete[] triangle_b;
    delete[] triangle_c;

    return init_time;
}

// Arbitrary Window Query
// Sliding Window Query
// Our
void our_insert_nbr_triangle(int32_t u, int32_t v, int32_t w, const int32_t &t, vector<map<int32_t, vector<my_pair>, greater<>>> &nbr_triangles)
{
    if (v > w)
        swap(v, w);

    if (nbr_triangles[u].find(t) == nbr_triangles[u].end())
    {
        nbr_triangles[u].emplace(t, vector<my_pair>{(my_pair){v, w}});
    }
    else
    {
        nbr_triangles[u][t].emplace_back((my_pair){v, w});
    }
}

// Arbitrary Window Query
// Our
int64_t historical_our_index(vector<vector<triplet>> &idx_a, vector<vector<triplet>> &idx_b, vector<map<int32_t, vector<my_pair>, greater<>>> &nbr_triangles)
{
    int64_t init_time = 0;

    auto t3 = chrono::high_resolution_clock::now();
    auto t4 = chrono::high_resolution_clock::now();

    int64_t idx_init_time = 0;
    int64_t idx_add_time = 0;
    int64_t idx_update_time = 0;

    t3 = chrono::high_resolution_clock::now();
    idx_a.resize(g_n);
    idx_b.resize(g_n);

    bool *visited = new bool[g_n]();
    bool *inf = new bool[g_n]();

    vector<vector<my_pair>> nbr(g_n);

    int32_t *triangle_a = new int32_t[g_n];
    int32_t *triangle_b = new int32_t[g_n];
    int32_t *triangle_c = new int32_t[g_n];
    memset(triangle_a, -1, sizeof(int32_t) * g_n);
    memset(triangle_b, -1, sizeof(int32_t) * g_n);
    memset(triangle_c, -1, sizeof(int32_t) * g_n);
    bool *triangle_d = new bool[g_n]();

    // hypotenuse
    int32_t *triangle_h = new int32_t[g_n];
    memset(triangle_h, -1, sizeof(int32_t) * g_n);
    bool *triangle_h_new = new bool[g_n]();
    t4 = chrono::high_resolution_clock::now();
    idx_init_time += chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();

    int32_t *ds_parent = new int32_t[g_n];
    int32_t *ds_size = new int32_t[g_n];
    memset(ds_parent, -1, sizeof(int32_t) * g_n);
    memset(ds_size, -1, sizeof(int32_t) * g_n);

    for (int32_t t_e = 0; t_e < g_t; ++t_e)
    {
        auto t1 = chrono::high_resolution_clock::now();

        vector<int32_t> visited_list;
        vector<int32_t> inf_vertices;

        t3 = chrono::high_resolution_clock::now();

        // (unnecessary) remove duplicate edges
        // because there should be no duplicate edge within the same time label after preprocessing the dataset
        unordered_set<pair<int32_t, int32_t>, boost::hash<pair<int32_t, int32_t>>> nbr_set;

        // update neighbors (add new edges)
        for (int32_t t = g_offset[t_e]; t < g_offset[t_e + 1]; ++t)
        {
            int32_t u = g_edges[t].first;
            int32_t v = g_edges[t].second;

            // (unnecessary) remove duplicate edges
            // because there should be no duplicate edge within the same time label after preprocessing the dataset
            if (u > v)
                swap(u, v);
            // check if edge exists
            if (nbr_set.find(make_pair(u, v)) != nbr_set.end())
                continue;
            // add new edges
            nbr_set.emplace(make_pair(u, v));

            if (!visited[u])
            {
                visited[u] = true;
                visited_list.emplace_back(u);
            }
            if (!visited[v])
            {
                visited[v] = true;
                visited_list.emplace_back(v);
            }

            nbr[u].emplace_back((my_pair){v, t_e});
            nbr[v].emplace_back((my_pair){u, t_e});
        }

        // update triangle
        // use bitmap: triangle_a_history, triangle_b_history, triangle_c_history, triangle_d_history;
        for (auto &u : visited_list)
        {
            vector<int32_t> triangle_a_history, triangle_b_history, triangle_c_history, triangle_d_history;

            for (auto item = nbr[u].rbegin(); item != nbr[u].rend(); ++item)
            {
                int32_t v = (*item).first;
                int32_t t = (*item).second;

                // identify new edges
                if (t == t_e)
                {
                    if (nbr[v].size() < nbr[u].size() || (nbr[v].size() == nbr[u].size() && v < u))
                    {
                        triangle_c[v] = t;
                        triangle_c_history.emplace_back(v);
                    }
                    else
                    {
                        triangle_d[v] = true;
                        triangle_d_history.emplace_back(v);
                    }
                    continue;
                }

                // skip B- edges
                if (triangle_d[v])
                    continue;

                // retrieve previous time of new edges
                // skip if already retrieved
                if (triangle_c[v] != -1)
                {
                    if (triangle_c[v] == t_e)
                        triangle_c[v] = t;
                    continue;
                }

                // retrieve latest time of old edges
                // skip if already retrieved
                if (triangle_a[v] != -1)
                    continue;
                if (triangle_b[v] != -1)
                    continue;

                // t is t_max for (u, v) as t is non-increasing in reverse direction of nbr[u]
                if (nbr[v].size() < nbr[u].size() || (nbr[v].size() == nbr[u].size() && v < u))
                {
                    triangle_a[v] = t;
                    triangle_a_history.emplace_back(v);
                }
                else
                {
                    triangle_b[v] = t;
                    triangle_b_history.emplace_back(v);
                }
            }

            // there is no duplicate in triangle_c_history
            for (auto &v : triangle_c_history)
            {
                vector<int32_t> triangle_h_history, triangle_h_new_history;
                for (auto item = nbr[v].rbegin(); item != nbr[v].rend(); ++item)
                {
                    int32_t w = (*item).first;
                    int32_t t = (*item).second;

                    // identify new edges
                    if (t == t_e)
                    {
                        triangle_h_new[w] = true;
                        triangle_h_new_history.emplace_back(w);
                        triangle_h[w] = t;
                        triangle_h_history.emplace_back(w);
                        continue;
                    }

                    // retrieve previous time of new edges
                    // skip if already retrieved
                    if (triangle_h_new[w])
                    {
                        if (triangle_h[w] == t_e)
                            triangle_h[w] = t;
                        continue;
                    }

                    // retrieve latest time of old edges
                    // skip if already retrieved
                    if (triangle_h[w] != -1)
                        continue;

                    triangle_h[w] = t;
                    triangle_h_history.emplace_back(w);
                }

                for (auto &w : triangle_h_history)
                {
                    int32_t t = triangle_h[w];

                    if (triangle_h_new[w])
                    {
                        // CASE 2.2
                        if (triangle_a[w] != -1)
                        {
                            if (triangle_c[v] != t_e && triangle_c[v] >= triangle_a[w] && t != t_e && t >= triangle_a[w])
                            {
                                continue;
                            }
                            // triangle (u,v,w)
                            // (u,v)=t_e, (v,w)=t_e, (u,w)=triangle_a[w]
                            int32_t triangle_time = min(t_e, triangle_a[w]);

                            our_insert_nbr_triangle(u, v, w, triangle_time, nbr_triangles);
                            our_insert_nbr_triangle(v, u, w, triangle_time, nbr_triangles);
                            our_insert_nbr_triangle(w, u, v, triangle_time, nbr_triangles);
                            check_inf(u, v, w, inf, inf_vertices);
                        }

                        // CASE 3
                        if (nbr[w].size() < nbr[v].size() || (nbr[w].size() == nbr[v].size() && w < v))
                        {
                            if (triangle_c[w] != -1)
                            {
                                // find triangle (u,v,w)
                                // (u,v)=t_e, (v,w)=t_e, (u,w)=t_e
                                our_insert_nbr_triangle(u, v, w, t_e, nbr_triangles);
                                our_insert_nbr_triangle(v, u, w, t_e, nbr_triangles);
                                our_insert_nbr_triangle(w, u, v, t_e, nbr_triangles);
                                check_inf(u, v, w, inf, inf_vertices);
                            }
                        }
                        continue;
                    }

                    // CASE 1
                    if (triangle_a[w] != -1 || triangle_b[w] != -1)
                    {
                        if (triangle_c[v] != t_e && (triangle_c[v] >= t || triangle_c[v] >= max(triangle_a[w], triangle_b[w])))
                        {
                            continue;
                        }
                        // find triangle (u,v,w)
                        // (u,v)=t_e, (v,w)=t, (u,w)=max(triangle_a[w],triangle_b[w])
                        // either triangle_a[w] = -1 or triangle_b[w] = -1
                        int32_t triangle_time = min(t, max(triangle_a[w], triangle_b[w]));

                        our_insert_nbr_triangle(u, v, w, triangle_time, nbr_triangles);
                        our_insert_nbr_triangle(v, u, w, triangle_time, nbr_triangles);
                        our_insert_nbr_triangle(w, u, v, triangle_time, nbr_triangles);
                        check_inf(u, v, w, inf, inf_vertices);
                    }

                    // CASE 2.1
                    if (nbr[w].size() < nbr[v].size() || (nbr[w].size() == nbr[v].size() && w < v))
                    {
                        if (triangle_c[w] != -1)
                        {
                            if (triangle_c[v] != t_e && triangle_c[v] >= t && triangle_c[w] != t_e && triangle_c[w] >= t)
                            {
                                continue;
                            }
                            // find triangle (u,v,w)
                            // (u,v)=t_e, (v,w)=t, (u,w)=t_e
                            our_insert_nbr_triangle(u, v, w, t, nbr_triangles);
                            our_insert_nbr_triangle(v, u, w, t, nbr_triangles);
                            our_insert_nbr_triangle(w, u, v, t, nbr_triangles);
                            check_inf(u, v, w, inf, inf_vertices);
                        }
                    }
                }

                for (auto &w : triangle_h_history)
                    triangle_h[w] = -1;
                for (auto &w : triangle_h_new_history)
                    triangle_h_new[w] = false;
            }

            for (auto &v : triangle_a_history)
                triangle_a[v] = -1;
            for (auto &v : triangle_b_history)
                triangle_b[v] = -1;
            for (auto &v : triangle_c_history)
                triangle_c[v] = -1;
            for (auto &v : triangle_d_history)
                triangle_d[v] = false;
        }
        t4 = chrono::high_resolution_clock::now();
        idx_add_time += chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();

        t3 = chrono::high_resolution_clock::now();
        // update idx_a, idx_b
        for (auto &u : inf_vertices)
        {
            vector<int32_t> ds_list;

            int32_t prev_a = 0;
            int32_t prev_b = 0;
            int32_t curr_a = 0;
            int32_t curr_b = 0;

            for (auto const &item : nbr_triangles[u])
            {
                int32_t t = item.first;
                for (auto const t_item : item.second)
                {
                    const int32_t &v = t_item.first;
                    const int32_t &w = t_item.second;

                    if (ds_parent[v] == -1)
                    {
                        ds_parent[v] = v;
                        ds_size[v] = 1;
                        ds_list.emplace_back(v);
                    }
                    if (ds_parent[w] == -1)
                    {
                        ds_parent[w] = w;
                        ds_size[w] = 1;
                        ds_list.emplace_back(w);
                    }

                    int32_t v_root = v;
                    while (ds_parent[v_root] != v_root)
                    {
                        ds_parent[v_root] = ds_parent[ds_parent[v_root]];
                        v_root = ds_parent[v_root];
                    }
                    int32_t w_root = w;
                    while (ds_parent[w_root] != w_root)
                    {
                        ds_parent[w_root] = ds_parent[ds_parent[w_root]];
                        w_root = ds_parent[w_root];
                    }
                    if (v_root == w_root)
                        continue;

                    if (ds_size[v_root] < ds_size[w_root])
                        swap(v_root, w_root);

                    ds_parent[w_root] = v_root;
                    ds_size[v_root] = ds_size[v_root] + ds_size[w_root];

                    ++curr_a;
                }

                if (curr_a != prev_a)
                {
                    prev_a = curr_a;
                    idx_a[u].emplace_back((triplet){t, t_e, curr_a});
                }

                curr_b = curr_a;
                for (auto const &ds_ele : ds_list)
                {
                    if (ds_ele == ds_parent[ds_ele] && ds_size[ds_ele] >= g_tau)
                    {
                        ++curr_b;
                    }
                }

                if (curr_b != prev_b)
                {
                    prev_b = curr_b;
                    idx_b[u].emplace_back((triplet){t, t_e, curr_b});
                }
            }

            for (auto const &ds_ele : ds_list)
            {
                ds_parent[ds_ele] = -1;
                ds_size[ds_ele] = -1;
            }
        }

        for (auto &u : visited_list)
            visited[u] = false;
        for (auto &u : inf_vertices)
            inf[u] = false;
        t4 = chrono::high_resolution_clock::now();
        idx_update_time += chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();

        auto t2 = chrono::high_resolution_clock::now();
        if (t_e < (int32_t)(g_t * g_sliding))
        {
            init_time += chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        }
    }

    delete[] visited;
    delete[] inf;
    delete[] triangle_a;
    delete[] triangle_b;
    delete[] triangle_c;
    delete[] triangle_d;
    delete[] triangle_h;
    delete[] triangle_h_new;
    delete[] ds_parent;
    delete[] ds_size;

    printf("Index(Init) -> %f s\n", idx_init_time / 1000000000.0);
    printf("Index(Add) -> %f s\n", idx_add_time / 1000000000.0);
    printf("Index(Update) -> %f s\n", idx_update_time / 1000000000.0);

    return init_time;
}

// Sliding Window Query
// Baseline
void sliding_baseline_t_insert(int32_t v, int32_t w, int32_t &result, unordered_map<int32_t, ds_item> &ds)
{
    if (ds.find(v) == ds.end())
        ds.emplace(v, (ds_item){v, 1});
    if (ds.find(w) == ds.end())
        ds.emplace(w, (ds_item){w, 1});

    int32_t sa = ds_find(v, ds);
    int32_t sb = ds_find(w, ds);
    if (sa != sb)
    {
        if (max(ds[sa].size, ds[sb].size) < g_tau && ds[sa].size + ds[sb].size >= g_tau)
        {
            result += 1;
        }
        else if (ds[sa].size >= g_tau && ds[sb].size >= g_tau)
        {
            result -= 1;
        }
        ds_union(v, w, ds);
    }
}

// Sliding Window Query
// Baseline
void sliding_baseline_insert_nbr_triangle(int32_t u, int32_t v, int32_t w, vector<set<pair<int32_t, int32_t>>> &nbr_triangles)
{
    if (v > w)
        swap(v, w);

    if (nbr_triangles[u].find(make_pair(v, w)) == nbr_triangles[u].end())
    {
        nbr_triangles[u].emplace(make_pair(v, w));
    }
    else
    {
        printf("ERROR: duplicate triangle (%d, %d, %d)\n", u, v, w);
    }
}

// Sliding Window Query
// Baseline
void sliding_baseline_remove_nbr_triangle(int32_t u, int32_t v, int32_t w, vector<set<pair<int32_t, int32_t>>> &nbr_triangles)
{
    if (v > w)
        swap(v, w);

    if (nbr_triangles[u].find(make_pair(v, w)) != nbr_triangles[u].end())
    {
        nbr_triangles[u].erase(make_pair(v, w));
    }
    else
    {
        printf("ERROR: triangle (%d, %d, %d) not found\n", u, v, w);
    }
}

// Sliding Window Query
// Baseline
int64_t sliding_baseline_index(vector<int32_t> &idx)
{
    int64_t init_time = 0;

    auto t3 = chrono::high_resolution_clock::now();
    auto t4 = chrono::high_resolution_clock::now();

    int64_t idx_init_time = 0;
    int64_t idx_add_time = 0;
    int64_t idx_del_time = 0;
    int64_t idx_update_time = 0;

    t3 = chrono::high_resolution_clock::now();
    idx.resize(g_n);
    for (int32_t i = 0; i < g_n; ++i)
    {
        idx[i] = 0;
    }

    bool *visited = new bool[g_n]();
    bool *inf = new bool[g_n]();

    // record edge count
    vector<map<int32_t, int32_t>> nbr_set(g_n);

    vector<vector<int32_t>> adj(g_n);
    vector<int32_t> adj_offset(g_n, 0);
    vector<set<pair<int32_t, int32_t>>> nbr_triangles(g_n);

    bool *triangle_a = new bool[g_n]();
    bool *triangle_b = new bool[g_n]();
    bool *triangle_c = new bool[g_n]();
    t4 = chrono::high_resolution_clock::now();
    idx_init_time += chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();

    for (int32_t t_e = 0; t_e < g_t; ++t_e)
    {
        auto t1 = chrono::high_resolution_clock::now();

        vector<int32_t> visited_list;
        vector<int32_t> inf_vertices;

        t3 = chrono::high_resolution_clock::now();
        // update neighbors (add new edges)
        for (int32_t t = g_offset[t_e]; t < g_offset[t_e + 1]; ++t)
        {
            int32_t u = g_edges[t].first;
            int32_t v = g_edges[t].second;

            // skip duplicate edges
            if (u > v)
                swap(u, v);
            // check if edge exists
            if (nbr_set[u].find(v) != nbr_set[u].end())
            {
                nbr_set[u][v] += 1;
                continue;
            }
            // add new edges
            nbr_set[u].emplace(v, 1);
            // nbr_set[v].emplace(u, 1);

            if (!visited[u])
            {
                visited[u] = true;
                visited_list.emplace_back(u);
            }

            if (!visited[v])
            {
                visited[v] = true;
                visited_list.emplace_back(v);
            }

            adj[u].emplace_back(v);
            adj[v].emplace_back(u);
        }

        // update triangle
        for (auto &u : visited_list)
        {
            vector<int32_t> ta, tb, tc;
            for (int32_t i = 0; i < adj_offset[u]; ++i)
            {
                int32_t v = adj[u][i];
                if (adj[v].size() < adj[u].size() || (adj[v].size() == adj[u].size() && v < u))
                {
                    triangle_a[v] = true;
                    ta.emplace_back(v);
                }
                else
                {
                    triangle_b[v] = true;
                    tb.emplace_back(v);
                }
            }
            for (int32_t i = adj_offset[u]; i < adj[u].size(); ++i)
            {
                int32_t v = adj[u][i];
                if (adj[v].size() < adj[u].size() || (adj[v].size() == adj[u].size() && v < u))
                {
                    triangle_c[v] = true;
                    tc.emplace_back(v);
                }
            }

            for (auto &v : tc)
            {
                for (int32_t i = 0; i < adj_offset[v]; ++i)
                {
                    int32_t w = adj[v][i];
                    // CASE 1
                    if (triangle_a[w] || triangle_b[w])
                    {
                        // find triangle (u,v,w)
                        // printf("<%d,%d,%d>\n",u,v,w);
                        sliding_baseline_insert_nbr_triangle(u, v, w, nbr_triangles);
                        sliding_baseline_insert_nbr_triangle(v, u, w, nbr_triangles);
                        sliding_baseline_insert_nbr_triangle(w, u, v, nbr_triangles);
                        check_inf(u, v, w, inf, inf_vertices);
                    }

                    // CASE 2.1
                    if (adj[w].size() < adj[v].size() || (adj[w].size() == adj[v].size() && w < v))
                    {
                        if (triangle_c[w])
                        {
                            // find trinalge (u,v,w)
                            // printf("<%d,%d,%d>\n",u,v,w);
                            sliding_baseline_insert_nbr_triangle(u, v, w, nbr_triangles);
                            sliding_baseline_insert_nbr_triangle(v, u, w, nbr_triangles);
                            sliding_baseline_insert_nbr_triangle(w, u, v, nbr_triangles);
                            check_inf(u, v, w, inf, inf_vertices);
                        }
                    }
                }

                for (int32_t i = adj_offset[v]; i < adj[v].size(); ++i)
                {
                    int32_t w = adj[v][i];
                    // CASE 2.2
                    if (triangle_a[w])
                    {
                        // find triangle (u,v,w)
                        // printf("<%d,%d,%d>\n",u,v,w);
                        sliding_baseline_insert_nbr_triangle(u, v, w, nbr_triangles);
                        sliding_baseline_insert_nbr_triangle(v, u, w, nbr_triangles);
                        sliding_baseline_insert_nbr_triangle(w, u, v, nbr_triangles);
                        check_inf(u, v, w, inf, inf_vertices);
                    }

                    // CASE 3
                    if (adj[w].size() < adj[v].size() || (adj[w].size() == adj[v].size() && w < v))
                    {
                        if (triangle_c[w])
                        {
                            // find trinalge (u,v,w)
                            // printf("<%d,%d,%d>\n",u,v,w);
                            sliding_baseline_insert_nbr_triangle(u, v, w, nbr_triangles);
                            sliding_baseline_insert_nbr_triangle(v, u, w, nbr_triangles);
                            sliding_baseline_insert_nbr_triangle(w, u, v, nbr_triangles);
                            check_inf(u, v, w, inf, inf_vertices);
                        }
                    }
                }
            }

            for (auto &v : ta)
                triangle_a[v] = false;
            for (auto &v : tb)
                triangle_b[v] = false;
            for (auto &v : tc)
                triangle_c[v] = false;
        }

        for (auto &u : visited_list)
            adj_offset[u] = adj[u].size();
        for (auto &u : visited_list)
            visited[u] = false;
        visited_list.clear();
        t4 = chrono::high_resolution_clock::now();
        idx_add_time += chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();

        t3 = chrono::high_resolution_clock::now();
        // calculate sliding window
        int32_t sliding = max(0, t_e - (int32_t)(g_t * g_sliding));

        if (sliding != 0)
        {
            vector<vector<int32_t>> del_adj(g_n);

            // update neighbors (remove old edges)
            for (int32_t t = g_offset[sliding - 1]; t < g_offset[sliding]; ++t)
            {
                int32_t u = g_edges[t].first;
                int32_t v = g_edges[t].second;

                // skip duplicate edges
                if (u > v)
                    swap(u, v);
                // check if edge exists
                if (nbr_set[u].find(v) != nbr_set[u].end())
                {
                    nbr_set[u][v] -= 1;
                    if (nbr_set[u][v] == 0)
                    {
                        // remove old edges
                        nbr_set[u].erase(v);
                        // nbr_set[v].erase(u);

                        if (!visited[u])
                        {
                            visited[u] = true;
                            visited_list.emplace_back(u);
                        }

                        if (!visited[v])
                        {
                            visited[v] = true;
                            visited_list.emplace_back(v);
                        }

                        adj[u].erase(remove(adj[u].begin(), adj[u].end(), v), adj[u].end());
                        adj[v].erase(remove(adj[v].begin(), adj[v].end(), u), adj[v].end());

                        del_adj[u].emplace_back(v);
                        del_adj[v].emplace_back(u);
                    }
                }
                else
                {
                    printf("ERROR: edge (%d, %d) does not exist\n", u, v);
                }
            }

            // update triangle
            for (auto &u : visited_list)
            {
                vector<int32_t> ta, tb, tc;
                for (int32_t i = 0; i < adj[u].size(); ++i)
                {
                    int32_t v = adj[u][i];
                    if (adj[v].size() + del_adj[v].size() < adj[u].size() + del_adj[u].size() || (adj[v].size() + del_adj[v].size() == adj[u].size() + del_adj[u].size() && v < u))
                    {
                        triangle_a[v] = true;
                        ta.emplace_back(v);
                    }
                    else
                    {
                        triangle_b[v] = true;
                        tb.emplace_back(v);
                    }
                }
                for (int32_t i = 0; i < del_adj[u].size(); ++i)
                {
                    int32_t v = del_adj[u][i];
                    if (adj[v].size() + del_adj[v].size() < adj[u].size() + del_adj[u].size() || (adj[v].size() + del_adj[v].size() == adj[u].size() + del_adj[u].size() && v < u))
                    {
                        triangle_c[v] = true;
                        tc.emplace_back(v);
                    }
                }

                for (auto &v : tc)
                {
                    for (int32_t i = 0; i < adj[v].size(); ++i)
                    {
                        int32_t w = adj[v][i];
                        // CASE 1
                        if (triangle_a[w] || triangle_b[w])
                        {
                            // find triangle (u,v,w)
                            // printf("<%d,%d,%d>\n",u,v,w);
                            sliding_baseline_remove_nbr_triangle(u, v, w, nbr_triangles);
                            sliding_baseline_remove_nbr_triangle(v, u, w, nbr_triangles);
                            sliding_baseline_remove_nbr_triangle(w, u, v, nbr_triangles);
                            check_inf(u, v, w, inf, inf_vertices);
                        }

                        // CASE 2.1
                        if (adj[w].size() + del_adj[w].size() < adj[v].size() + del_adj[v].size() || (adj[w].size() + del_adj[w].size() == adj[v].size() + del_adj[v].size() && w < v))
                        {
                            if (triangle_c[w])
                            {
                                // find trinalge (u,v,w)
                                // printf("<%d,%d,%d>\n",u,v,w);
                                sliding_baseline_remove_nbr_triangle(u, v, w, nbr_triangles);
                                sliding_baseline_remove_nbr_triangle(v, u, w, nbr_triangles);
                                sliding_baseline_remove_nbr_triangle(w, u, v, nbr_triangles);
                                check_inf(u, v, w, inf, inf_vertices);
                            }
                        }
                    }

                    for (int32_t i = 0; i < del_adj[v].size(); ++i)
                    {
                        int32_t w = del_adj[v][i];
                        // CASE 2.2
                        if (triangle_a[w])
                        {
                            // find triangle (u,v,w)
                            // printf("<%d,%d,%d>\n",u,v,w);
                            sliding_baseline_remove_nbr_triangle(u, v, w, nbr_triangles);
                            sliding_baseline_remove_nbr_triangle(v, u, w, nbr_triangles);
                            sliding_baseline_remove_nbr_triangle(w, u, v, nbr_triangles);
                            check_inf(u, v, w, inf, inf_vertices);
                        }

                        // CASE 3
                        if (adj[w].size() + del_adj[w].size() < adj[v].size() + del_adj[v].size() || (adj[w].size() + del_adj[w].size() == adj[v].size() + del_adj[v].size() && w < v))
                        {
                            if (triangle_c[w])
                            {
                                // find trinalge (u,v,w)
                                // printf("<%d,%d,%d>\n",u,v,w);
                                sliding_baseline_remove_nbr_triangle(u, v, w, nbr_triangles);
                                sliding_baseline_remove_nbr_triangle(v, u, w, nbr_triangles);
                                sliding_baseline_remove_nbr_triangle(w, u, v, nbr_triangles);
                                check_inf(u, v, w, inf, inf_vertices);
                            }
                        }
                    }
                }

                for (auto &v : ta)
                    triangle_a[v] = false;
                for (auto &v : tb)
                    triangle_b[v] = false;
                for (auto &v : tc)
                    triangle_c[v] = false;
            }
        }
        t4 = chrono::high_resolution_clock::now();
        idx_del_time += chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();

        t3 = chrono::high_resolution_clock::now();
        // update idx
        for (auto &u : inf_vertices)
        {
            unordered_map<int32_t, ds_item> ds;
            int32_t sd = 0;

            for (auto const &item : nbr_triangles[u])
            {
                const int32_t &v = item.first;
                const int32_t &w = item.second;
                sliding_baseline_t_insert(v, w, sd, ds);
            }

            idx[u] = sd;
        }
        t4 = chrono::high_resolution_clock::now();
        idx_update_time += chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();

        for (auto &u : visited_list)
            adj_offset[u] = adj[u].size();
        for (auto &u : visited_list)
            visited[u] = false;
        for (auto &u : inf_vertices)
            inf[u] = false;

        auto t2 = chrono::high_resolution_clock::now();
        if (t_e < (int32_t)(g_t * g_sliding))
        {
            init_time += chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        }
    }

    delete[] visited;
    delete[] inf;
    delete[] triangle_a;
    delete[] triangle_b;
    delete[] triangle_c;

    printf("Index(Init) -> %f s\n", idx_init_time / 1000000000.0);
    printf("Index(Add) -> %f s\n", idx_add_time / 1000000000.0);
    printf("Index(Del) -> %f s\n", idx_del_time / 1000000000.0);
    printf("Index(Update) -> %f s\n", idx_update_time / 1000000000.0);

    // compute index size
    double idx_size = g_n * 2;
    for (int32_t i = 0; i < g_n; ++i)
    {
        idx_size += nbr_triangles[i].size() * 2;
    }
    printf("Index(Base) -> %.2f MB\n", idx_size * 4 / 1024 / 1024);

    // compute index size
    idx_size = g_n * 2;
    for (int32_t i = 0; i < g_n; ++i)
    {
        idx_size += nbr_triangles[i].size() * 3;
    }
    printf("Index(Ours) -> %.2f MB\n", idx_size * 4 / 1024 / 1024);

    return init_time;
}

// Sliding Window Query
// Our
int64_t sliding_our_index(vector<deque<my_pair>> &idx_a, vector<deque<my_pair>> &idx_b)
{
    int64_t init_time = 0;

    auto t3 = chrono::high_resolution_clock::now();
    auto t4 = chrono::high_resolution_clock::now();

    int64_t idx_init_time = 0;
    int64_t idx_add_time = 0;
    int64_t idx_del_time = 0;
    int64_t idx_update_time = 0;

    t3 = chrono::high_resolution_clock::now();
    idx_a.resize(g_n);
    idx_b.resize(g_n);

    bool *visited = new bool[g_n]();
    bool *inf = new bool[g_n]();
    bool *del = new bool[g_n]();
    bool *del_inf = new bool[g_n]();

    vector<vector<my_pair>> nbr(g_n);
    vector<map<int32_t, vector<my_pair>, greater<>>> nbr_triangles(g_n);

    int32_t *triangle_a = new int32_t[g_n];
    int32_t *triangle_b = new int32_t[g_n];
    int32_t *triangle_c = new int32_t[g_n];
    memset(triangle_a, -1, sizeof(int32_t) * g_n);
    memset(triangle_b, -1, sizeof(int32_t) * g_n);
    memset(triangle_c, -1, sizeof(int32_t) * g_n);
    bool *triangle_d = new bool[g_n]();

    // hypotenuse
    int32_t *triangle_h = new int32_t[g_n];
    memset(triangle_h, -1, sizeof(int32_t) * g_n);
    bool *triangle_h_new = new bool[g_n]();
    t4 = chrono::high_resolution_clock::now();
    idx_init_time += chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();

    int32_t *ds_parent = new int32_t[g_n];
    int32_t *ds_size = new int32_t[g_n];
    memset(ds_parent, -1, sizeof(int32_t) * g_n);
    memset(ds_size, -1, sizeof(int32_t) * g_n);

    for (int32_t t_e = 0; t_e < g_t; ++t_e)
    {
        auto t1 = chrono::high_resolution_clock::now();

        vector<int32_t> visited_list;
        vector<int32_t> inf_vertices;
        vector<int32_t> del_vertices;
        vector<int32_t> del_inf_vertices;

        t3 = chrono::high_resolution_clock::now();
        // calculate sliding window
        int32_t sliding = max(0, t_e - (int32_t)(g_t * g_sliding));

        if (sliding != 0)
        {
            // record vertices included in removed edges
            for (int32_t t = g_offset[sliding - 1]; t < g_offset[sliding]; ++t)
            {
                int32_t u = g_edges[t].first;
                int32_t v = g_edges[t].second;

                if (!del[u])
                {
                    del[u] = true;
                    del_vertices.emplace_back(u);

                    // update neighbors (remove old edges)
                    for (auto it = nbr[u].begin(); it != nbr[u].end();)
                    {
                        if (it->second < sliding)
                        {
                            it = nbr[u].erase(it);
                        }
                        else
                        {
                            break;
                        }
                    }

                    // record vertices influenced by removed edges
                    if (!del_inf[u])
                    {
                        del_inf[u] = true;
                        del_inf_vertices.emplace_back(u);
                    }
                    for (auto const &item : nbr[u])
                    {
                        if (!del_inf[item.first])
                        {
                            del_inf[item.first] = true;
                            del_inf_vertices.emplace_back(item.first);
                        }
                    }
                }
                if (!del[v])
                {
                    del[v] = true;
                    del_vertices.emplace_back(v);

                    // update neighbors (remove old edges)
                    for (auto it = nbr[v].begin(); it != nbr[v].end();)
                    {
                        if (it->second < sliding)
                        {
                            it = nbr[v].erase(it);
                        }
                        else
                        {
                            break;
                        }
                    }

                    // record vertices influenced by removed edges
                    if (!del_inf[v])
                    {
                        del_inf[v] = true;
                        del_inf_vertices.emplace_back(v);
                    }
                    for (auto const &item : nbr[v])
                    {
                        if (!del_inf[item.first])
                        {
                            del_inf[item.first] = true;
                            del_inf_vertices.emplace_back(item.first);
                        }
                    }
                }
            }
        }

        // remove old triangles from nbr_triangles and remove old indices from idx_a, idx_b
        for (auto const &u : del_inf_vertices)
        {
            // remove old triangles from nbr_triangles
            for (auto it = nbr_triangles[u].begin(); it != nbr_triangles[u].end();)
            {
                if (it->first < sliding)
                {
                    it = nbr_triangles[u].erase(it);
                }
                else
                {
                    break;
                }
            }

            // remove old indices from idx_a, idx_b
            while (!idx_a[u].empty() && idx_a[u].back().first < sliding)
            {
                idx_a[u].pop_back();
            }
            while (!idx_b[u].empty() && idx_b[u].back().first < sliding)
            {
                idx_b[u].pop_back();
            }
        }
        t4 = chrono::high_resolution_clock::now();
        idx_del_time += chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();

        t3 = chrono::high_resolution_clock::now();

        // (unnecessary) remove duplicate edges
        // because there should be no duplicate edge within the same time label after preprocessing the dataset
        unordered_set<pair<int32_t, int32_t>, boost::hash<pair<int32_t, int32_t>>> nbr_set;

        // update neighbors (add new edges)
        for (int32_t t = g_offset[t_e]; t < g_offset[t_e + 1]; ++t)
        {
            int32_t u = g_edges[t].first;
            int32_t v = g_edges[t].second;

            // (unnecessary) remove duplicate edges
            // because there should be no duplicate edge within the same time label after preprocessing the dataset
            if (u > v)
                swap(u, v);
            // check if edge exists
            if (nbr_set.find(make_pair(u, v)) != nbr_set.end())
                continue;
            // add new edges
            nbr_set.emplace(make_pair(u, v));

            if (!visited[u])
            {
                visited[u] = true;
                visited_list.emplace_back(u);
            }
            if (!visited[v])
            {
                visited[v] = true;
                visited_list.emplace_back(v);
            }

            nbr[u].emplace_back((my_pair){v, t_e});
            nbr[v].emplace_back((my_pair){u, t_e});
        }

        // update triangle
        // use bitmap: triangle_a_history, triangle_b_history, triangle_c_history, triangle_d_history;
        for (auto &u : visited_list)
        {
            vector<int32_t> triangle_a_history, triangle_b_history, triangle_c_history, triangle_d_history;

            for (auto item = nbr[u].rbegin(); item != nbr[u].rend(); ++item)
            {
                int32_t v = (*item).first;
                int32_t t = (*item).second;

                // identify new edges
                if (t == t_e)
                {
                    if (nbr[v].size() < nbr[u].size() || (nbr[v].size() == nbr[u].size() && v < u))
                    {
                        triangle_c[v] = t;
                        triangle_c_history.emplace_back(v);
                    }
                    else
                    {
                        triangle_d[v] = true;
                        triangle_d_history.emplace_back(v);
                    }
                    continue;
                }

                // skip B- edges
                if (triangle_d[v])
                    continue;

                // retrieve previous time of new edges
                // skip if already retrieved
                if (triangle_c[v] != -1)
                {
                    if (triangle_c[v] == t_e)
                        triangle_c[v] = t;
                    continue;
                }

                // retrieve latest time of old edges
                // skip if already retrieved
                if (triangle_a[v] != -1)
                    continue;
                if (triangle_b[v] != -1)
                    continue;

                // t is t_max for (u, v) as t is non-increasing in reverse direction of nbr[u]
                if (nbr[v].size() < nbr[u].size() || (nbr[v].size() == nbr[u].size() && v < u))
                {
                    triangle_a[v] = t;
                    triangle_a_history.emplace_back(v);
                }
                else
                {
                    triangle_b[v] = t;
                    triangle_b_history.emplace_back(v);
                }
            }

            // there is no duplicate in triangle_c_history
            for (auto &v : triangle_c_history)
            {
                vector<int32_t> triangle_h_history, triangle_h_new_history;
                for (auto item = nbr[v].rbegin(); item != nbr[v].rend(); ++item)
                {
                    int32_t w = (*item).first;
                    int32_t t = (*item).second;

                    // identify new edges
                    if (t == t_e)
                    {
                        triangle_h_new[w] = true;
                        triangle_h_new_history.emplace_back(w);
                        triangle_h[w] = t;
                        triangle_h_history.emplace_back(w);
                        continue;
                    }

                    // retrieve previous time of new edges
                    // skip if already retrieved
                    if (triangle_h_new[w])
                    {
                        if (triangle_h[w] == t_e)
                            triangle_h[w] = t;
                        continue;
                    }

                    // retrieve latest time of old edges
                    // skip if already retrieved
                    if (triangle_h[w] != -1)
                        continue;

                    triangle_h[w] = t;
                    triangle_h_history.emplace_back(w);
                }

                for (auto &w : triangle_h_history)
                {
                    int32_t t = triangle_h[w];

                    if (triangle_h_new[w])
                    {
                        // CASE 2.2
                        if (triangle_a[w] != -1)
                        {
                            if (triangle_c[v] != t_e && triangle_c[v] >= triangle_a[w] && t != t_e && t >= triangle_a[w])
                            {
                                continue;
                            }
                            // triangle (u,v,w)
                            // (u,v)=t_e, (v,w)=t_e, (u,w)=triangle_a[w]
                            int32_t triangle_time = min(t_e, triangle_a[w]);

                            our_insert_nbr_triangle(u, v, w, triangle_time, nbr_triangles);
                            our_insert_nbr_triangle(v, u, w, triangle_time, nbr_triangles);
                            our_insert_nbr_triangle(w, u, v, triangle_time, nbr_triangles);
                            check_inf(u, v, w, inf, inf_vertices);
                        }

                        // CASE 3
                        if (nbr[w].size() < nbr[v].size() || (nbr[w].size() == nbr[v].size() && w < v))
                        {
                            if (triangle_c[w] != -1)
                            {
                                // find triangle (u,v,w)
                                // (u,v)=t_e, (v,w)=t_e, (u,w)=t_e
                                our_insert_nbr_triangle(u, v, w, t_e, nbr_triangles);
                                our_insert_nbr_triangle(v, u, w, t_e, nbr_triangles);
                                our_insert_nbr_triangle(w, u, v, t_e, nbr_triangles);
                                check_inf(u, v, w, inf, inf_vertices);
                            }
                        }
                        continue;
                    }

                    // CASE 1
                    if (triangle_a[w] != -1 || triangle_b[w] != -1)
                    {
                        if (triangle_c[v] != t_e && (triangle_c[v] >= t || triangle_c[v] >= max(triangle_a[w], triangle_b[w])))
                        {
                            continue;
                        }
                        // find triangle (u,v,w)
                        // (u,v)=t_e, (v,w)=t, (u,w)=max(triangle_a[w],triangle_b[w])
                        // either triangle_a[w] = -1 or triangle_b[w] = -1
                        int32_t triangle_time = min(t, max(triangle_a[w], triangle_b[w]));

                        our_insert_nbr_triangle(u, v, w, triangle_time, nbr_triangles);
                        our_insert_nbr_triangle(v, u, w, triangle_time, nbr_triangles);
                        our_insert_nbr_triangle(w, u, v, triangle_time, nbr_triangles);
                        check_inf(u, v, w, inf, inf_vertices);
                    }

                    // CASE 2.1
                    if (nbr[w].size() < nbr[v].size() || (nbr[w].size() == nbr[v].size() && w < v))
                    {
                        if (triangle_c[w] != -1)
                        {
                            if (triangle_c[v] != t_e && triangle_c[v] >= t && triangle_c[w] != t_e && triangle_c[w] >= t)
                            {
                                continue;
                            }
                            // find triangle (u,v,w)
                            // (u,v)=t_e, (v,w)=t, (u,w)=t_e
                            our_insert_nbr_triangle(u, v, w, t, nbr_triangles);
                            our_insert_nbr_triangle(v, u, w, t, nbr_triangles);
                            our_insert_nbr_triangle(w, u, v, t, nbr_triangles);
                            check_inf(u, v, w, inf, inf_vertices);
                        }
                    }
                }

                for (auto &w : triangle_h_history)
                    triangle_h[w] = -1;
                for (auto &w : triangle_h_new_history)
                    triangle_h_new[w] = false;
            }

            for (auto &v : triangle_a_history)
                triangle_a[v] = -1;
            for (auto &v : triangle_b_history)
                triangle_b[v] = -1;
            for (auto &v : triangle_c_history)
                triangle_c[v] = -1;
            for (auto &v : triangle_d_history)
                triangle_d[v] = false;
        }
        t4 = chrono::high_resolution_clock::now();
        idx_add_time += chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();

        t3 = chrono::high_resolution_clock::now();
        // update idx_a, idx_b
        for (auto &u : inf_vertices)
        {
            idx_a[u].clear();
            idx_b[u].clear();

            vector<int32_t> ds_list;

            int32_t prev_a = 0;
            int32_t prev_b = 0;
            int32_t curr_a = 0;
            int32_t curr_b = 0;

            for (auto const &item : nbr_triangles[u])
            {
                int32_t t = item.first;
                for (auto const t_item : item.second)
                {
                    const int32_t &v = t_item.first;
                    const int32_t &w = t_item.second;

                    if (ds_parent[v] == -1)
                    {
                        ds_parent[v] = v;
                        ds_size[v] = 1;
                        ds_list.emplace_back(v);
                    }
                    if (ds_parent[w] == -1)
                    {
                        ds_parent[w] = w;
                        ds_size[w] = 1;
                        ds_list.emplace_back(w);
                    }

                    int32_t v_root = v;
                    while (ds_parent[v_root] != v_root)
                    {
                        ds_parent[v_root] = ds_parent[ds_parent[v_root]];
                        v_root = ds_parent[v_root];
                    }
                    int32_t w_root = w;
                    while (ds_parent[w_root] != w_root)
                    {
                        ds_parent[w_root] = ds_parent[ds_parent[w_root]];
                        w_root = ds_parent[w_root];
                    }
                    if (v_root == w_root)
                        continue;

                    if (ds_size[v_root] < ds_size[w_root])
                        swap(v_root, w_root);

                    ds_parent[w_root] = v_root;
                    ds_size[v_root] = ds_size[v_root] + ds_size[w_root];

                    ++curr_a;
                }

                if (curr_a != prev_a)
                {
                    prev_a = curr_a;
                    idx_a[u].emplace_back((my_pair){t, curr_a});
                }

                curr_b = curr_a;
                for (auto const &ds_ele : ds_list)
                {
                    if (ds_ele == ds_parent[ds_ele] && ds_size[ds_ele] >= g_tau)
                    {
                        ++curr_b;
                    }
                }

                if (curr_b != prev_b)
                {
                    prev_b = curr_b;
                    idx_b[u].emplace_back((my_pair){t, curr_b});
                }
            }

            for (auto const &ds_ele : ds_list)
            {
                ds_parent[ds_ele] = -1;
                ds_size[ds_ele] = -1;
            }
        }
        t4 = chrono::high_resolution_clock::now();
        idx_update_time += chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();

        for (auto &u : visited_list)
            visited[u] = false;
        for (auto &u : inf_vertices)
            inf[u] = false;
        for (auto &u : del_vertices)
            del[u] = false;
        for (auto &u : del_inf_vertices)
            del_inf[u] = false;

        auto t2 = chrono::high_resolution_clock::now();
        if (t_e < (int32_t)(g_t * g_sliding))
        {
            init_time += chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        }
    }

    delete[] visited;
    delete[] inf;
    delete[] del;
    delete[] del_inf;
    delete[] triangle_a;
    delete[] triangle_b;
    delete[] triangle_c;
    delete[] triangle_d;
    delete[] triangle_h;
    delete[] triangle_h_new;
    delete[] ds_parent;
    delete[] ds_size;

    printf("Index(Init) -> %f s\n", idx_init_time / 1000000000.0);
    printf("Index(Add) -> %f s\n", idx_add_time / 1000000000.0);
    printf("Index(Del) -> %f s\n", idx_del_time / 1000000000.0);
    printf("Index(Update) -> %f s\n", idx_update_time / 1000000000.0);

    return init_time;
}

int historical_our_binary_search_border_ts(vector<triplet> &idx, const int &l, const int &r, const int &n, const int &target)
{
    if (l == r)
    {
        if (idx[l].first >= target)
            return l;
        return n;
    }

    int left = l, right = r, ans = n;
    while (left <= right)
    {
        int mid = ((right - left) >> 1) + left;
        if (idx[mid].first >= target)
        {
            left = mid + 1;
            ans = mid;
        }
        else
        {
            right = mid - 1;
        }
    }
    return ans;
}

int historical_our_binary_search_border_te(vector<triplet> &idx, const int &n, const int &target, const bool &lower)
{
    int left = 0, right = n - 1, ans = n;
    while (left <= right)
    {
        int mid = ((right - left) >> 1) + left;
        if (idx[mid].second > target || (lower && idx[mid].second >= target))
        {
            right = mid - 1;
            ans = mid;
        }
        else
        {
            left = mid + 1;
        }
    }
    return ans;
}

int historical_our_search_index(vector<triplet> &idx, const int &t_s, const int &t_e)
{
    int n = idx.size();

    // find the largest end time <= t_e
    int bound1 = historical_our_binary_search_border_te(idx, n, t_e, true);

    // bound1 as right bound
    if (bound1 == n)
    {
        --bound1;
        // bound2 as left bound
        int bound2 = historical_our_binary_search_border_te(idx, n, idx[bound1].second, true);
        // find the smallest start time >= t_s and corresponding index value
        int pos = historical_our_binary_search_border_ts(idx, bound2, bound1, bound1 + 1, t_s);
        if (pos == bound1 + 1)
            return 0;
        return idx[pos].third;
    }

    // bound1 as left bound
    if (idx[bound1].second == t_e)
    {
        // bound2 as right bound
        int bound2 = historical_our_binary_search_border_te(idx, n, idx[bound1].second, false) - 1;
        // find the smallest start time >= t_s and corresponding index value
        int pos = historical_our_binary_search_border_ts(idx, bound1, bound2, bound2 + 1, t_s);
        if (pos == bound2 + 1)
            return 0;
        return idx[pos].third;
    }

    // none
    if (bound1 == 0)
        return 0;

    // bound1 as right bound
    --bound1;
    // bound2 as left bound
    int bound2 = historical_our_binary_search_border_te(idx, n, idx[bound1].second, true);
    int nn = bound1 - bound2 + 1;
    // find the smallest start time >= t_s and corresponding index value
    int pos = historical_our_binary_search_border_ts(idx, bound2, bound1, bound1 + 1, t_s);
    if (pos == bound1 + 1)
        return 0;
    return idx[pos].third;
}

int32_t historical_bl_query_one(vector<vector<triplet>> &idx, int32_t u, int32_t t_s, int32_t t_e)
{
    int32_t sd = 0;
    if (idx[u].size() != 0)
        sd = historical_our_search_index(idx[u], t_s, t_e);
    return sd;
}

int32_t historical_our_query_one(vector<vector<triplet>> &idx_a, vector<vector<triplet>> &idx_b, int32_t u, int32_t t_s, int32_t t_e)
{
    int32_t sd = 0;
    // each vertex needs 3 binary searches
    // two find left and right borders for the largest end time <= t_e
    // one find within two borders for the smallest start time >= t_s
    if (idx_b[u].size() != 0)
        sd = historical_our_search_index(idx_b[u], t_s, t_e);
    if (idx_a[u].size() != 0)
        sd -= historical_our_search_index(idx_a[u], t_s, t_e);
    return sd;
}

int64_t compute_nbr_triangles_size(vector<map<int32_t, vector<my_pair>, greater<>>> &nbr_triangles)
{
    vector<map<int32_t, vector<my_pair>, greater<>>> final_nbr_triangles(g_n);
    for (int32_t u = 0; u < g_n; ++u)
    {
        unordered_set<pair<int32_t, int32_t>, boost::hash<pair<int32_t, int32_t>>> nbr_set;
        for (auto &it : nbr_triangles[u])
        {
            int32_t t = it.first;
            for (auto &jt : it.second)
            {
                int32_t v = jt.first;
                int32_t w = jt.second;
                if (v > w)
                {
                    swap(v, w);
                }
                if (nbr_set.find(make_pair(v, w)) != nbr_set.end())
                    continue;
                nbr_set.emplace(make_pair(v, w));
                if (final_nbr_triangles[u].find(t) == final_nbr_triangles[u].end())
                {
                    final_nbr_triangles[u].emplace(t, vector<my_pair>{(my_pair){v, w}});
                }
                else
                {
                    final_nbr_triangles[u][t].emplace_back((my_pair){v, w});
                }
            }
        }
    }
    int64_t idx_size = g_n;
    for (int32_t u = 0; u < g_n; ++u)
    {
        idx_size += final_nbr_triangles[u].size();
        for (auto &it : final_nbr_triangles[u])
        {
            idx_size += it.second.size() * 2;
        }
    }
    return idx_size;
}

int32_t sliding_baseline_query_one(vector<int32_t> &idx, int32_t u)
{
    return idx[u];
}

int32_t sliding_our_query_one(vector<deque<my_pair>> &idx_a, vector<deque<my_pair>> &idx_b, int32_t u)
{
    int32_t sd = 0;
    if (idx_b[u].size() != 0)
        sd = idx_b[u].back().second;
    if (idx_a[u].size() != 0)
        sd -= idx_a[u].back().second;
    return sd;
}

int main(int argc, char *argv[])
{
    if (argc < 6 || argc > 7)
    {
        printf("usage: %s [setting] [method] [interval] [graph] [tau] [theta](optional)\n", argv[0]);
        printf("[setting] 1 -> Arbitrary Window Query\n");
        printf("[setting] 2 -> Sliding Window Query\n");
        printf("[method] 1 -> Online\n");
        printf("[method] 2 -> Baseline\n");
        printf("[method] 3 -> Our\n");
        return 0;
    }

    int32_t setting = atoi(argv[1]);
    int32_t method = atoi(argv[2]);
    g_t_interval = atoi(argv[3]);
    load_graph(argv[4]);
    g_tau = atoi(argv[5]);

    // g_t = t_max
    g_t = g_offset.size() - 1;
    printf("t_max = %d ", g_t);

    if (argc == 6)
    {
        printf("%s %s %s %s %s %s\n", argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]);
    }
    else if (argc == 7)
    {
        g_sliding = atof(argv[6]); // set sliding window size to g_sliding * g_t
        printf("%s %s %s %s %s %s %s\n", argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
    }

    auto t1 = chrono::high_resolution_clock::now();
    auto t2 = chrono::high_resolution_clock::now();

    if (method == 1)
    {
        // Online
        int64_t query_time = 0;
        for (int32_t i = 0; i < 1000; ++i)
        {
            // t_e - t_s = g_sliding * g_t
            int32_t t_s = rand() % (int32_t)(g_t * (1.0 - g_sliding));
            int32_t t_e = t_s + (int32_t)(g_t * g_sliding);
            int32_t u = rand() % g_n;

            t1 = chrono::high_resolution_clock::now();
            int64_t scan_time = online_query_one(t_s, t_e, u);
            t2 = chrono::high_resolution_clock::now();
            query_time += chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() - scan_time;
        }
        printf("Query -> %f ns\n", query_time / 1000.0);
        return 0;
    }

    // idx_a -> SNC hierarchy
    // idx_b -> NC hierarchy

    if (setting == 1)
    {
        // Arbitrary Window Query
        if (method == 2)
        {
            // Baseline
            t1 = chrono::high_resolution_clock::now();
            auto t3 = chrono::high_resolution_clock::now();
            vector<vector<triplet>> idx;
            auto t4 = chrono::high_resolution_clock::now();
            int64_t init_time = chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();
            init_time += historical_bl_index(idx);
            t2 = chrono::high_resolution_clock::now();
            // int64_t idx_time = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() - init_time;
            int64_t idx_time = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
            printf("Index -> %f s\n", idx_time / 1000000000.0);

            // compute index size
            double idx_size = g_n;
            for (int32_t i = 0; i < g_n; ++i)
            {
                idx_size += idx[i].size() * 3;
            }
            printf("Index -> %.2f MB\n", idx_size * 4 / 1024 / 1024);

            int64_t query_time = 0;
            for (int32_t i = 0; i < 1000; ++i)
            {
                // t_e - t_s = g_sliding * g_t
                int32_t t_s = rand() % (int32_t)(g_t * (1.0 - g_sliding));
                int32_t t_e = t_s + (int32_t)(g_t * g_sliding);
                int32_t u = rand() % g_n;

                t1 = chrono::high_resolution_clock::now();
                historical_bl_query_one(idx, u, t_s, t_e);
                t2 = chrono::high_resolution_clock::now();
                query_time += chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
            }
            printf("Query -> %f ns\n", query_time / 1000.0);

            // compute t_base
            double t_base = 0;
            for (auto &i : idx)
                t_base += i.size();
            printf("t_base: %.2f\n", t_base / g_n);
        }
        else if (method == 3)
        {
            // Our
            t1 = chrono::high_resolution_clock::now();
            auto t3 = chrono::high_resolution_clock::now();
            vector<vector<triplet>> idx_a;
            vector<vector<triplet>> idx_b;
            vector<map<int32_t, vector<my_pair>, greater<>>> nbr_triangles(g_n);
            auto t4 = chrono::high_resolution_clock::now();
            int64_t init_time = chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();
            init_time += historical_our_index(idx_a, idx_b, nbr_triangles);
            t2 = chrono::high_resolution_clock::now();
            // int64_t idx_time = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() - init_time;
            int64_t idx_time = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
            printf("Index -> %f s\n", idx_time / 1000000000.0);

            // compute index size
            double idx_size = g_n * 2;
            for (int32_t i = 0; i < g_n; ++i)
            {
                idx_size += idx_a[i].size() * 3 + idx_b[i].size() * 3;
            }
            idx_size += compute_nbr_triangles_size(nbr_triangles);
            printf("Index -> %.2f MB\n", idx_size * 4 / 1024 / 1024);

            int64_t query_time = 0;
            for (int32_t i = 0; i < 1000; ++i)
            {
                // t_e - t_s = g_sliding * g_t
                int32_t t_s = rand() % (int32_t)(g_t * (1.0 - g_sliding));
                int32_t t_e = t_s + (int32_t)(g_t * g_sliding);
                int32_t u = rand() % g_n;

                t1 = chrono::high_resolution_clock::now();
                historical_our_query_one(idx_a, idx_b, u, t_s, t_e);
                t2 = chrono::high_resolution_clock::now();
                query_time += chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
            }
            printf("Query -> %f ns\n", query_time / 1000.0);

            // compute t_overline
            unordered_set<int32_t> t_overline;
            double t_overline_a = 0;
            double t_overline_b = 0;
            for (auto &i : idx_a)
            {
                for (auto &tri : i)
                {
                    if (!t_overline.count(tri.second))
                        t_overline.emplace(tri.second);
                }
                t_overline_a += t_overline.size();
                t_overline.clear();
            }
            for (auto &i : idx_b)
            {
                for (auto &tri : i)
                {
                    if (!t_overline.count(tri.second))
                        t_overline.emplace(tri.second);
                }
                t_overline_b += t_overline.size();
                t_overline.clear();
            }
            printf("t_overline: %.2f\n", (t_overline_a + t_overline_b) / (2 * g_n));
        }
    }
    else if (setting == 2)
    {
        // Sliding Window Query
        if (method == 2)
        {
            // Baseline
            t1 = chrono::high_resolution_clock::now();
            vector<int32_t> idx;
            int64_t init_time = sliding_baseline_index(idx);
            t2 = chrono::high_resolution_clock::now();
            int64_t idx_time = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() - init_time;
            printf("Index -> %f s\n", idx_time / 1000000000.0);

            int64_t query_time = 0;
            for (int32_t i = 0; i < 1000; ++i)
            {
                int32_t u = rand() % g_n;

                t1 = chrono::high_resolution_clock::now();
                sliding_baseline_query_one(idx, u);
                t2 = chrono::high_resolution_clock::now();
                query_time += chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
            }
            printf("Query -> %f ns\n", query_time / 1000.0);
        }
        else if (method == 3)
        {
            // Our
            t1 = chrono::high_resolution_clock::now();
            vector<deque<my_pair>> idx_a;
            vector<deque<my_pair>> idx_b;
            int64_t init_time = sliding_our_index(idx_a, idx_b);
            t2 = chrono::high_resolution_clock::now();
            int64_t idx_time = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() - init_time;
            printf("Index -> %f s\n", idx_time / 1000000000.0);

            int64_t query_time = 0;
            for (int32_t i = 0; i < 1000; ++i)
            {
                int32_t u = rand() % g_n;

                t1 = chrono::high_resolution_clock::now();
                sliding_our_query_one(idx_a, idx_b, u);
                t2 = chrono::high_resolution_clock::now();
                query_time += chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
            }
            printf("Query -> %f ns\n", query_time / 1000.0);
        }
    }

    return 0;
}
