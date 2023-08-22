template<const int P> struct mint
{
    int v;

    mint(const int x=0): v((x%P+P)%P) {}
    mint(const long long x): v(int((x%P+P)%P)) {}

    template<typename T> mint& operator=(const T a) {*this=mint(a); return (*this);}
    mint operator+(const mint a) const {const int r=v+a.v; return mint(r<P?r:r-P);}
    mint operator-(const mint a) const {const int r=v-a.v; return mint(r<0?r+P:r);}
    mint operator*(const mint a) const {return mint(1ll*v*a.v);}
    mint& operator+=(const mint a) {(*this)=(*this)+a; return (*this);}
    mint& operator-=(const mint a) {(*this)=(*this)-a; return (*this);}
    mint& operator*=(const mint a) {(*this)=(*this)*a; return (*this);}
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
    mint operator/(const mint a) const {return (*this)*a.inv();}
    mint& operator/=(const mint b) {(*this)=(*this)/b; return (*this);}
};