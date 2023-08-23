#include <bits/stdc++.h>
using namespace std;

template<int P> struct mint;

constexpr int P1=233,P2=239;
constexpr int M1=1004535809,M2=2013265921;

template<int P,int M> struct Hash
{
    vector<mint<M>> pow,h;

    Hash(const string &s)
    {
        pow.resize(s.size()+1);
        h.resize(s.size());
        pow[0]=1;
        for (size_t i=0;i<s.size();i++)
        {
            pow[i+1]=pow[i]*P;
            h[i]=i==0?s[i]:h[i-1]*P+s[i];
        }
    }

    mint<M> get(const int l,const int r) const
    {
        if (l==0) return h[r];
        return h[r]-h[l-1]*pow[r-l+1];
    }
};
