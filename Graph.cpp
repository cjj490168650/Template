#include <bits/stdc++.h>
using namespace std;

template<typename T> struct graph
{
    int n,m=0;
    vector<vector<int>> e;
    vector<int> to;
    vector<T> wt;

    graph(const int n): n(n),e(n) {}

    void add_arc(const int u,const int v,const T w=0) {e[u].push_back(m++); to.push_back(v); wt.push_back(w);}
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

template<typename T> struct tree: graph<T>
{
    using graph<T>::n,graph<T>::e,graph<T>::to;
    using graph<T>::add_edge;

    int maxk;
    vector<vector<int>> fa;
    vector<int> dep;

    tree(const int n,const int maxk=25): graph<T>(n),maxk(maxk),fa(n,vector<int>(maxk+1,-1)),dep(n) {}

    void dfs(const int u,const int f)
    {
        fa[u][0]=f;
        dep[u]=f==-1?0:dep[f]+1;
        for (int i=1;i<=maxk;i++) fa[u][i]=fa[u][i-1]==-1?-1:fa[fa[u][i-1]][i-1];
        for (const int i:e[u])
        {
            const int v=to[i];
            if (v==f) continue;
            dfs(v,u);
        }
    }

    int lca(int u,int v) const
    {
        if (dep[u]<dep[v]) swap(u,v);
        for (int i=maxk;i>=0;i--)
        {
            if (fa[u][i]!=-1 && dep[fa[u][i]]>=dep[v]) u=fa[u][i];
        }
        if (u==v) return u;
        for (int i=maxk;i>=0;i--)
        {
            if (fa[u][i]!=fa[v][i]) u=fa[u][i],v=fa[v][i];
        }
        return fa[u][0];
    }
};

template<typename T> struct tree: graph<T>
{
    using graph<T>::n,graph<T>::e,graph<T>::to,graph<T>::wt;
    using graph<T>::add_edge;

    int clk=0;
    vector<int> fa,dep,siz,son,top,dfn,rnk;

    tree(const int n): graph<T>(n),fa(n),dep(n),siz(n),son(n,-1),top(n),dfn(n),rnk(n) {}

    void dfs1(const int u,const int f)
    {
        fa[u]=f;
        dep[u]=f==-1?0:dep[f]+1;
        siz[u]=1;
        for (const int i:e[u])
        {
            const int v=to[i];
            if (v==f) continue;
            dfs1(v,u);
            siz[u]+=siz[v];
            if (son[u]==-1 || siz[v]>siz[son[u]]) son[u]=v;
        }
    }

    void dfs2(const int u,const int t)
    {
        top[u]=t;
        dfn[u]=clk++;
        rnk[dfn[u]]=u;
        if (son[u]==-1) return;
        dfs2(son[u],t);
        for (const int i:e[u])
        {
            const int v=to[i];
            if (v==fa[u] || v==son[u]) continue;
            dfs2(v,v);
        }
    }

    T query(int u,int v) const
    {
        T ans=0;
        while (top[u]!=top[v])
        {
            if (dep[top[u]]<dep[top[v]]) swap(u,v);
            // query(dfn[top[u]],dfn[u]);
            u=fa[top[u]];
        }
        if (dep[u]>dep[v]) swap(u,v);
        // query(dfn[u],dfn[v]);
        return ans;
    }
};