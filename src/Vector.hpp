template<int n, type Scalar>
class Vector {
    private:
        Scalar *d;
    public:
        Vector()                   : d(new Scalar[n]) {}
        ~Vector()                    {delete[] d;}
        Vector(const Vector& rhs)    : a(new Scalar[n]) {std::copy(rhs.d,rhs.d+n,d);}
        void swap(Vector& rhs)       {std::swap(d,rhs.d);}
        Vector& operator=(const Vector& rhs) 
                                   {Vector tmp(rhs); swap(rhs); return *this;}
};
