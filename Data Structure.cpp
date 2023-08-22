#include <bits/stdc++.h>
using namespace std;

struct unionfind
{
    vector<int> f,siz;

    unionfind(const int n): f(n),siz(n,1) {for (int i=0;i<n;i++) f[i]=i;}

    int find(const int u) {return f[u]==u?u:f[u]=find(f[u]);}

    void merge(const int u,const int v)
    {
        int fu=find(u),fv=find(v);
        if (fu==fv) return;
        if (siz[fu]>siz[fv]) swap(fu,fv);
        f[fu]=fv; siz[fv]+=siz[fu];
    }
};

template<typename T,class cmp=less<T>> struct sttable
{
    int n;
    vector<vector<T>> st;
    vector<int> lg2;

    sttable(const vector<T> &a): n(a.size()),st(a.size()),lg2(a.size()+1)
    {
        lg2[0]=-1;
        for (int i=1;i<=n;i++) lg2[i]=lg2[i>>1]+1;
        const int K=lg2.back();
        for (int i=0;i<n;i++) st[i].resize(K+1),st[i][0]=a[i];
        for (int k=1;k<=K;k++)
        {
            for (int i=0;i+(1<<k)-1<n;i++)
            {
                st[i][k]=min(st[i][k-1],st[i+(1<<(k-1))][k-1],cmp());
            }
        }
    }

    T query(const int l,const int r) const
    {
        const int k=lg2[r-l+1];
        return min(st[l][k],st[r-(1<<k)+1][k],cmp());
    }
};

template<typename T> struct fenwick
{
    int n;
    vector<T> val;

    fenwick(const int n): n(n),val(n) {}

    int lowbit(const int x) const {return x&(-x);}
    void update(int p,const T v) {while (p<n) val[p]+=v,p+=lowbit(p+1);}
    T query(int p) const {T re=0; while (p>=0) re+=val[p],p-=lowbit(p+1); return re;}
    T query(const int l,const int r) const {return query(r)-query(l-1);}
};

template<typename T> struct segtree
{
    struct segnode
    {
        int l,r;
        T sum;
        T tmul=1,tadd=0;

        segnode(const int i=0,const T v=0): l(i),r(i),sum(v) {}
        segnode(const int l,const int r,const T x): l(l),r(r),sum(x) {}

        segnode operator+ (const segnode &u) const {return segnode(l,u.r,sum+u.sum);}

        void update(const T mul,const T add)
        {
            sum=sum*mul+add*(r-l+1);
            tadd=tadd*mul+add;
            tmul=tmul*mul;
        }
    };

    vector<segnode> tr;

    segtree(const int n): tr((n+1)<<2) {}

    void push_up(const int id)
    {
        const int lson=(id<<1),rson=(id<<1)|1;
        tr[id]=tr[lson]+tr[rson];
    }

    void push_down(const int id)
    {
        const int lson=(id<<1),rson=(id<<1)|1;
        tr[lson].update(tr[id].tmul,tr[id].tadd);
        tr[rson].update(tr[id].tmul,tr[id].tadd);
        tr[id].tadd=0;
        tr[id].tmul=1;
    }

    void build(const int id,const int l,const int r,const vector<T> &init)
    {
        const int lson=(id<<1),rson=(id<<1)|1,mid=(l+r)>>1;
        if (l==r)
        {
            if (init.empty()) tr[id]=segnode(l);
            else tr[id]=segnode(l,init[l]);
            return;
        }
        build(lson,l,mid,init);
        build(rson,mid+1,r,init);
        push_up(id);
        tr[id].l=l; tr[id].r=r;
    }
    void build(const int l,const int r,const vector<T> &init={}) {build(1,l,r,init);}

    void update(const int id,const int L,const int R,const T mul,const T add)
    {
        const int l=tr[id].l,r=tr[id].r;
        const int lson=(id<<1),rson=(id<<1)|1,mid=(l+r)>>1;
        if (l<r) push_down(id);
        if (L<=l && R>=r)
        {
            tr[id].update(mul,add);
            return;
        }
        if (mid>=R) update(lson,L,R,mul,add);
        else if (mid<L) update(rson,L,R,mul,add);
        else
        {
            update(lson,L,mid,mul,add);
            update(rson,mid+1,R,mul,add);
        }
        push_up(id);
    }
    void update(const int l,const int r,const T mul,const T add) {update(1,l,r,mul,add);}

    segnode query(const int id,const int L,const int R)
    {
        const int l=tr[id].l,r=tr[id].r;
        const int lson=(id<<1),rson=(id<<1)|1,mid=(l+r)>>1;
        if (l<r) push_down(id);
        if (L<=l && R>=r) return tr[id];
        if (mid>=R) return query(lson,L,R);
        if (mid<L) return query(rson,L,R);
        return query(lson,L,mid)+query(rson,mid+1,R);
    }
    segnode query(const int l,const int r) {return query(1,l,r);}
};

template<typename T> struct BST
{
    struct Node
    {
        T val;
        int cnt,siz,lson,rson;
        
        Node(T val): val(val),cnt(1),siz(1),lson(-1),rson(-1) {}
    };

    vector<Node> nodes;
    int rt=-1;

    int siz(int u) const {return u==-1?0:nodes[u].siz;}
    void maintain(int u) {nodes[u].siz=siz(nodes[u].lson)+siz(nodes[u].rson)+nodes[u].cnt;}

    int newnode(T val)
    {
        nodes.push_back(Node(val));
        return nodes.size()-1;
    }

    int find(T val) const
    {
        int p=rt;
        while (p!=-1 && nodes[p].val!=val)
        {
            if (nodes[p].val>val) p=nodes[p].lson;
            else p=nodes[p].rson;
        }
        return p;
    }

    int pre(T val) const
    {
        int p=rt,res=-1;
        while (p!=-1)
        {
            if (nodes[p].val<val) res=p,p=nodes[p].rson;
            else p=nodes[p].lson;
        }
        return res;
    }

    int nxt(T val) const
    {
        int p=rt,res=-1;
        while (p!=-1)
        {
            if (nodes[p].val>val) res=p,p=nodes[p].lson;
            else p=nodes[p].rson;
        }
        return res;
    }

    int rank(T val) const
    {
        int p=rt,k=1;
        while (p!=-1)
        {
            
            if (nodes[p].val==val) return siz(nodes[p].lson)+k;
            if (nodes[p].val>val) p=nodes[p].lson;
            else k+=siz(p)-siz(nodes[p].rson),p=nodes[p].rson;
        }
        return k;
    }

    int kth(int k) const
    {
        if (siz(rt)<k) return -1;
        int p=rt;
        while (p!=-1)
        {
            if (nodes[p].lson!=-1 && k<=siz(nodes[p].lson)) p=nodes[p].lson;
            else if (nodes[p].rson!=-1 && k>siz(p)-siz(nodes[p].rson)) k-=siz(p)-siz(nodes[p].rson),p=nodes[p].rson;
            else break;
        }
        return p;
    }
};

template<typename T> struct Treap: BST<T>
{
    using BST<T>::nodes,BST<T>::rt;
    using BST<T>::siz,BST<T>::maintain,BST<T>::newnode;
    using BST<T>::find,BST<T>::pre,BST<T>::nxt,BST<T>::rank,BST<T>::kth;

    mt19937 rnd;

    bool check(int u,int v) {return int(rnd()%(siz(u)+siz(v)))<siz(u);}

    pair<int,int> split(int u,T val)
    {
        if (u==-1) return {-1,-1};
        if (nodes[u].val>val)
        {
            const auto [r1,r2]=split(nodes[u].lson,val);
            nodes[u].lson=r2; maintain(u);
            return {r1,u};
        }
        else
        {
            const auto [r1,r2]=split(nodes[u].rson,val);
            nodes[u].rson=r1; maintain(u);
            return {u,r2};
        }
    }

    int merge(int u,int v)
    {
        if (u==-1) return v;
        if (v==-1) return u;
        if (check(u,v))
        {
            nodes[u].rson=merge(nodes[u].rson,v); maintain(u);
            return u;
        }
        else
        {
            nodes[v].lson=merge(u,nodes[v].lson); maintain(v);
            return v;
        }
    }

    void insert(T val)
    {
        auto [u1,u2]=split(rt,val-1);
        auto [v1,v2]=split(u2,val);
        if (v1==-1) v1=merge(v1,newnode(val));
        else nodes[v1].cnt++,maintain(v1);
        u2=merge(v1,v2); rt=merge(u1,u2);
    }

    void erase(T val)
    {
        auto [u1,u2]=split(rt,val-1);
        auto [v1,v2]=split(u2,val);
        if (v1==-1) return;
        nodes[v1].cnt--; maintain(v1);
        if (!nodes[v1].cnt) v1=-1;
        u2=merge(v1,v2); rt=merge(u1,u2);
    }
};