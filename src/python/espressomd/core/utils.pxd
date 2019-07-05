cdef extern from "utils/List.hpp" namespace "Utils":
    cppclass List[T]:
        List()
        List(size_t)
        List(size_t, const T &)

        T & operator[](size_t)
        void resize(size_t)
        void push_back(size_t)

        T * data()
        size_t size()

        T * e
        size_t n

cdef extern from "utils/Span.hpp" namespace "Utils":
    cppclass Span[T]:
        Span()
        Span(T * , size_t)

        T & operator[](size_t)

        T * begin()
        T * end()

        T * data()
        size_t size()

    Span[const T] make_const_span[T](T * , size_t)

cdef extern from "utils/Vector.hpp" namespace "Utils":
    cppclass Vector2d:
        pass
    cppclass Vector4d:
        pass

    cppclass Vector3i:
        int & operator[](int i)
        int * data()

    cppclass Vector3d:
        Vector3d()
        Vector3d(const Vector3d &)

        double & operator[](int i)
        double * data()
        Vector3d operator * (double i)
        Vector3d operator / (double i)

    cppclass Vector6d:
        double & operator[](int i)
        double * data()
        Vector6d operator * (double i)
        Vector6d operator / (double i)

    cppclass Vector9d:
        double & operator[](int i)
        double * data()
        Vector9d operator * (double i)
        Vector9d operator / (double i)

    cppclass Vector19d:
        double & operator[](int i)
        double * data()

