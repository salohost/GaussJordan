# distutils: language = c++
from sympy import Rational
from IPython.display import display, Math, Latex

cdef class PyMatrix:
	cdef Matrix[Frac] *c_matrix

	def __cinit__(self, L):
		self.c_matrix = new Matrix[Frac](len(L), len(L[0]))

		for i in range(self.c_matrix[0].rows()):
			for j in range(self.c_matrix[0].cols()):
				r = Rational(L[i][j])
				self.c_matrix[0].set(i, j, Frac(r.p, r.q))

	def list(self):
		L = []
		for i in range(self.c_matrix[0].rows()):
			L.append([])
			for j in range(self.c_matrix[0].cols()):
				L[i].append(Rational(
					self.c_matrix[0](i, j).m(),
					self.c_matrix[0](i, j).n()
					))
		return L

	def solve(self, b):
		if(len(b) != self.c_matrix.rows()):
			raise ValueError("Sizes don't match")
		cdef vector[Frac] vec
		for element in b:
			element = Rational(element)
			vec.push_back(Frac(element.p, element.q))
		display(Math(self.c_matrix.tex_solution(vec).decode('UTF-8')))

	def __dealloc__(self):
		del self.c_matrix

	def __str__(self):
		return self.c_matrix.str().decode('UTF-8')

	def _repr_latex_(self):
		return str(self)

	# def __repr__(self):
		# return str(self.c_matrix.print())





