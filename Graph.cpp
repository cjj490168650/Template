#include <bits/stdc++.h>
using namespace std;

template<typename T> struct graph
{
    int n,m=0;
    vector<vector<int>> e;
    vector<int> to;
    vector<T> wt;

    graph(const int n): n(n),e(n) {}

    void add_arc(const int u,const int v,const T w=0) {e[u].push_back(m); to.push_back(v); wt.push_back(w); m++;}
    void add_edge(const int u,const int v,const T w=0) {add_arc(u,v,w); add_arc(v,u,w);}

    vector<T> dijkstra(const int st) const
    {
        priority_queue<pair<T,int>,vector<pair<T,int>>,greater<pair<T,int>>> q;
        vector<T> dis(n,-1);
        dis[st]=0;
        q.push({0,st});
        while (!q.empty())
        {
            const auto [d,u]=q.top(); q.pop();
            if (d>dis[u]) continue;
            for (const int i:e[u])
            {
                const int v=to[i];
                const T w=wt[i];
                if (dis[v]==-1 || d+w<dis[v])
                {
                    dis[v]=d+w;
                    q.push(make_pair(dis[v],v));
                }
            }
        }
        return dis;
    }
};