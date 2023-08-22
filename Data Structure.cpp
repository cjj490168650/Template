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