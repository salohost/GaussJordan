from libcpp.string cimport string


cdef extern from "<vector>" namespace "std":
	cdef cppclass vector[T]:
		vector()
		void push_back(T&)



cdef extern from "Frac.cpp":
	cdef cppclass Frac:
		Frac()
		Frac(int)
		Frac(int, unsigned)
		int m()
		unsigned n()


cdef extern from "vMatrix.cpp":
	cdef cppclass Matrix[T]:
		Matrix()
		Matrix(size_t, size_t)
		T& operator ()(size_t, size_t)
		size_t rows()
		size_t cols()
		string str()
		void set(size_t, size_t, T)
		string tex_solution(vector[T])


