#include <bits/stdc++.h>
using namespace std;

template<int P> struct mint
{
    int v;

    mint(const int x=0): v(x) {}

    template<typename T> mint& operator=(const T &a) {*this=mint(a); return (*this);}
    bool operator==(const mint &a) const {return v==a.v;}
    mint operator+(const mint &a) const {const int r=v+a.v; return mint(r<P?r:r-P);}
    mint operator-(const mint &a) const {const int r=v-a.v; return mint(r<0?r+P:r);}
    mint operator-() const {return mint(P-v);}
    mint operator*(const mint &a) const {return mint(1ll*v*a.v%P);}
    mint& operator+=(const mint &a) {(*this)=(*this)+a; return (*this);}
    mint& operator-=(const mint &a) {(*this)=(*this)-a; return (*this);}
    mint& operator*=(const mint &a) {(*this)=(*this)*a; return (*this);}
    mint pow(long long k) const
    {
        mint r=1,a=(*this);
        while (k)
        {
            if (k&1) r*=a;
            a*=a; k>>=1;
        }
        return r;
    }
    mint inv() const {return pow(P-2);}
    mint operator/(const mint &a) const {return (*this)*a.inv();}
    mint& operator/=(const mint &b) {(*this)=(*this)/b; return (*this);}
};

template<int P> struct Matrix
{
    int n;
    vector<vector<mint<P>>> v;

    Matrix(int n, int d=0): n(n),v(n,vector<mint<P>>(n))
    {
        if (d) for (int i=0;i<n;i++) v[i][i]=d;
    }

    Matrix operator* (const Matrix &a) const
    {
        Matrix re(n);
        for (int i=0;i<n;i++)
        {
            for (int j=0;j<n;j++)
            {
                for (int k=0;k<n;k++)
                {
                    re.v[i][j]+=v[i][k]*a.v[k][j];
                }
            }
        }
        return re;
    }
    Matrix operator*= (const Matrix &a) {(*this)=(*this)*a; return (*this);}

    mint<P> elim(Matrix &b)
    {
        mint<P> re=1;
        Matrix a=*this;
        for (int i=0;i<n;i++)
        {
            if (a.v[i][i].v==0)
            {
                for (int j=i+1;j<n;j++)
                {
                    if (a.v[j][i].v)
                    {
                        swap(a.v[i],a.v[j]);
                        swap(b.v[i],b.v[j]);
                        re=-re;
                        break;
                    }
                }
                if (a.v[i][i].v==0) return 0;
            }
            const mint<P> x=a.v[i][i],invx=x.inv();
            re*=x;
            for (int j=0;j<n;j++)
            {
                a.v[i][j]*=invx;
                b.v[i][j]*=invx;
            }
            for (int k=0;k<n;k++)
            {
                if (k==i) continue;
                const mint<P> y=a.v[k][i];
                for (int j=0;j<n;j++)
                {
                    a.v[k][j]-=a.v[i][j]*y;
                    b.v[k][j]-=b.v[i][j]*y;
                }
            }
        }
        return re;
    }

    mint<P> det()
    {
        Matrix e(n);
        return elim(e);
    }

    Matrix inv()
    {
        Matrix e(n,1);
        const mint<P> d=elim(e);
        if (d.v==0) return Matrix(0);
        return e;
    }

    Matrix pow(long long k)
    {
        Matrix r(n,1),a=(*this);
        while (k)
        {
            if (k&1) r*=a;
            a*=a; k>>=1;
        }
        return r;
    }
};