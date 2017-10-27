[TOC]

# Введение

Решение линейных систем алгебраических уравнений одна из основных задач линейной алгебры, возникающая во многих отраслях: в экономике, физике, химии, математическом моделировании. Данная задача эквивалентна поиску прообраза вектора некоторого линейного преобразования. Соответсвенно решение, если оно есть, может быть вектором, в случае если преобразование инъективно или некоторым линейным многообразием, порожденным ядром линейного преобразования, в случае, когда ядро не тривиально. 

Соответсвенно любой СЛАУ:
$$
\begin{cases}
	a_{11}x_1+a_{12}x_2+...+a_{1n}x_n = b_1 \\
	a_{21}x_1+a_{22}x_2+...+a_{2n}x_n = b_2 \\
	... \\
	a_{m1}x_1+a_{m2}x_2+...a_{mn}x_n=b_m
\end{cases}
$$
Мы можем сопоставить матричное уравнение:
$$
Ax=b
$$
Где A - матрица составленная из коэффицентов СЛАУ, а $b$ - вектор составленный из $\{b_i\}$.

Для нахождения $x$ воспользуемся методом Гаусса-Жордана. Т.е. приведем расширенную матрицу $[A|b]$  элементарными преобразованиями над строками к каноническому ступенчатому виду $[A'|b']$, из которого далее легко извлечь ответ: (если система не вырождена)
$$
x \in X = b' + L(Kern(A'))
$$
где $L(Kern(A'))$ - линейная оболочка над ядром.

Для этого реализуем класс на C++ для матрицы над произвольным полем(с помощью шаблонов). И в качестве примера поля возьмем рациональные числа - реализуем класс для дробей(в отличии от float и double, позволяет избежать погрешностей).

В качестве полидрома для испытаний я выбрал Jupter Notebook, т.к. это удачная интерактивная среда, которая умеет отображать tex формулы, что намного интереснее, чем смотреть на матрицу из дробей в обычной консоли. Поэтому так-же далее будет написан класс-обертка на cython для создания Python модуля и дополнительные методы для генерации tex выражений.


# Классы
## Макросы
Определим следующий макрос - синтаксический сахар для for:

~~~c
#define range(i, begin, end) for(size_t i = (begin); i < (end); i++)
~~~
## Матрица
### Приватные поля и методы

~~~c++
template <typename Field>
class Matrix {
    vector< vector<Field>* > _M;
    Matrix& self = *this;
    
    void clear() {
        for(auto& row: _M) delete row;
    }
    
public:
...
~~~

Решение использовать массив из указателей на вектор связаны с частой необходимостью переставлять строки местами. Гораздо проще переставить местами 2 указателя, чем перемещать целые массивы.

Поле **self** используется исключительно для удобства обращения к объекту:

~~~c++
(*this)(i, j); // не очень
self(i, j);    // уже лучше
~~~


И, так как освобождать память нужно будет не только в конструкторе, объявлен **clear()**


### Публичные поля и методы
#### Конструкторы и деструктор
##### По размеру матрицы
~~~c++
    Matrix(size_t rows, size_t cols) {
        _M.resize(rows);
        for(auto& row: _M) row = new vector<Field>(cols);
    }
~~~

##### copy конструктор
Полезно уметь копировать объект.

```Описание assign оператора = будет ниже.```

~~~c++
    Matrix(const Matrix &A) {
        self = A;
    }
~~~

##### move конструктор
Дублировать объект, который скоро будет уничтожен, глупо. Поэтому мы просто приватизируем себе его данные.

~~~c++
    Matrix(Matrix&& A) {
        _M = A._M;
        A._M.resize(0);
    }
~~~
##### деструктор
~~~c++
    ~Matrix() {
        clear();
    }
~~~

#### Assignment operator
##### copy assigment

~~~c++
    Matrix& operator = (const Matrix& A) {
        clear();
        _M.resize(A.rows());
        range(i, 0, A.rows()) _M[i] = new vector<Field>(*(A._M[i]));
        return self;
    }
~~~
##### move assigment

~~~c++
    Matrix& operator = (Matrix&& A) {
        clear();
        _M = A._M;
        A._M.resize(0);
        return self;
    }
~~~

####Основные методы

##### Размер матрицы

Тяжело работать с матрицами, размеры которых не знаешь:

~~~c++
    size_t rows() const {
        return _M.size();
    }
    
    size_t cols() const {
        if(rows() == 0) return 0;
        return _M[0]->size();
    }
~~~

Полезно будет иметь двумерную индексацию:

~~~c++
    Field& operator() (size_t i, size_t j) const {
        return (*_M[i])[j];
    }
~~~

##### Перегруженные операторы:
С арифметикой здесь все очевидно. А вот в поток я вывожу latex формулу для матрицы.

~~~c++
    friend Matrix operator * (const Matrix& A, const Matrix& B) {
        if(A.cols() != B.rows()) {
            throw "Sizes don't match";
        }
        Matrix C(A.rows(), B.cols());
        range(i, 0, A.rows()) {
            range(j, 0, B.cols()) {
                range(k, 0, A.cols()) {
                    C(i, j) += A(i, k) * B(k, j);
                }
            }
        }
        return C;
    }
~~~


~~~c++
    friend Matrix operator + (const Matrix& A, const Matrix& B) {
        if( (A.cols() != B.cols()) && (A.rows() != B.rows()) ) {
            throw "Sizes don't match";
        }
        Matrix C(A.rows(), B.cols());
        range(i, 0, A.rows()) {
            range(j, 0, A.cols()) {
                C(i, j) = A(i, j) + B(i, j);
            }
        }
        return C;
    }
~~~


~~~c++
    friend Matrix operator - (const Matrix& A, const Matrix& B) {
        // return A + (-B);
        if( (A.cols() != B.cols()) || (A.rows() != B.rows()) ) {
            throw "Sizes don't match";
        }
        Matrix C(A.rows(), B.cols());
        range(i, 0, A.rows()) {
            range(j, 0, A.cols()) {
                C(i, j) = A(i, j) - B(i, j);
            }
        }
        return C;
        
    }
~~~

~~~c++
    friend Matrix operator - (const Matrix& A) {
        Matrix B(A);
        range(i, 0, B.rows()) {
            range(j, 0, B.cols()) {
                B(i, j) = -B(i, j);
            }
        }
        return B;
    }
~~~

~~~c++
    friend bool operator == (const Matrix& A, const Matrix& B) {
        if( (A.cols() != B.cols()) || (A.rows() != B.rows()) ) {
            throw "Sizes don't match";
        }
        range(i, 0, A.rows()) {
            range(j, 0, A.cols()) {
                if(A(i, j) != B(i, j)) return false;
            }
        }
        return true;
    }
~~~

~~~c++
    friend bool operator != (const Matrix& A, const Matrix& B) {
        return !(A == B);
    }
~~~

~~~c++
    friend istream& operator >> (istream& in, Matrix &M) {
        range(i, 0, M.rows())
            range(j, 0, M.cols())
                in >> M(i,j);
    }
~~~

~~~c++
    friend ostream& operator << (ostream& out, const Matrix &M) {
        out << "\\begin{bmatrix}\n";
        range(i, 0, M.rows()) {
            out << "    ";
            range(j, 0, M.cols()) {
                out << M(i, j) << (j < M.cols()-1 ? "&" : "");
            }
            out << (i < M.rows()-1 ? "\\\\\n" : "\n");
        }
        out << "\\end{bmatrix}";
        return out;
    }
~~~

~~~c++
    string str() {
        stringstream out;
        out << self;
        return out.str();
    }
~~~

##### Элементарные преобразования со строками:
Самые необходимые методы, для приведения матрицы к каноническому ступенчатому виду.  
Методы для сложения, перестановки, умножения и деления на скаляр строк:

~~~c++
    void add(size_t to, size_t from, Field k) {
        range(j, 0, cols()) self(to, j) += k * self(from, j);
    }
    
    void swap(size_t first, size_t second) {
        std::swap(_M[first], _M[second]);
    }
    
    void mul(size_t row, Field k) {
        range(j, 0, cols()) self(row,j) = self(row, j) * k;
    }
    
    void div(size_t row, Field k) {
        range(j, 0, cols()) self(row, j) = self(row,j) / k;
    }
~~~

##### Приведение матрицы медотом Гаусса-Жордана:
Самый важный метод, ради чего все и задумывалось:

В каждом столбце, начиная с первого(нулевого), ищем, так называемый, очередной опорный элемент с номером *p*. Т.е. находим ненулевой элемент(причем в такой строке, номер которой ниже текущего *p*) в столбце и элементарными преобразованиями обнуляем другие значения в столбце и меняем строку с этим элементом (исключительно для удобства) с p-ой строкой местами(это дает нам право не искать следующий опорный элемент в [0..p-1] строках). Если такой элемент найти невозможно, то поиск продолжается, аналогично, в следующем столбце, а текущий вектор столбец(с поставленной -1 на место соотсветсвущему номеру текущего столбца) далее попадет в nullspace нашей матрицы, т.е. в ядро соответсвующего линейного преобразования, т.е. в пространство решений СЛАУ, если система совместна.

Если в результате приведения оказалось, что напротив нулевой строки матрицы стоит ненулевой элемент столбца расширенной матрицы, то, очевидно, что решения нет, т.к. не существует линейной комбинации из 0 порождающих ненулевое значение.

Иначе решение есть и единственно, в случае если получившаяся матрица - единичная матрица(без учета нулевых строк) или иначе решение - это линейное многообразие, проходящее через вектор-столбец расширенной матрицы и паралелльно пространству решений однородной системы 
$$
\{b'+n|n\in Kern(A)\}
$$
Функция возвращает вердикт о том, есть решение или нет в виде bool и возвращает bool вектор, хранящий информацию о том, какие столбцы войдут а какие нет в nullspace. Значения возвращаются в паре:

~~~c++
    pair< bool, vector<bool> > row_reduce(vector<Field>& b) {
        if(b.size() != rows()) {
            throw "length of b don't match with quantity of rows";
        }
        vector<bool> pivots(cols(), false);
        
        size_t pivot = 0;
        range(k, 0, cols()) {
            range(i, pivot, rows()) {
                if(self(i,k) != Field()) {
                    pivots[k] = true;
                    
                    // swap
                    std::swap(b[pivot], b[i]);
                    swap(pivot, i);
                    
                    // div
                    b[pivot] /= self(pivot, k);
                    div(pivot, self(pivot, k));
                    
                    // обнуляем все элементы в столбце, кроме опорного
                    range(ii, 0, rows()) {
                        if(i != ii) {
                            b[ii] -= b[pivot]*self(ii, k);
                            add(ii, pivot, -self(ii,k));
                        }
                    }
                    pivot++;
                    break;
                }
            }
        }
        bool there_is_a_solution = true;
        range(i, pivot, rows()) {
            if(b[i] != Field()) {
                there_is_a_solution = false;
            }
        }
        
        return make_pair(there_is_a_solution, pivots);
    }
~~~

##### Извлечение ответа
Приведенная к каноническому ступенчатому виду матрица - это конечно здорово. Но хотелось бы видеть ответ как-то так:

$$
X = 

\begin{bmatrix}

    \frac{-1}{3}\\
    \frac{2}{3}\\
    0

\end{bmatrix}

- \begin{bmatrix}
   -1\\
   2\\
   -1
  \end{bmatrix}C_{1}
$$
Воспользуемся tex для генерации ответа:

~~~c++
    string tex_solution(vector<Field>& b) {
        Matrix A(self);
        
        stringstream out;
        
        bool is_there_a_solution;
        vector<bool> pivots;
        
        tie(is_there_a_solution, pivots) = A.row_reduce(b);
        
        if(!is_there_a_solution) {
            out << "There is no soulution\n";
            return out.str();
        }
        
        out << "X = \n";
        out << "\\begin{bmatrix}\n";
        
        range(j, 0, A.cols()) {
            out << "    " << (j < b.size() ? b[j] : Field())
            << ((j != A.cols()-1) ? "\\\\\n" : "\n");
        }
        
        out << "\\end{bmatrix}\n";
        
        size_t numeration = 0;
        range(k, 0, pivots.size()) {
            if(!pivots[k]) {
                out << " + \n" << "\\begin{bmatrix}\n";
                range(i, 0, A.cols()) {
                    out 
                    << "    "
                    << ((i == k) ? -Field(1): ( i < A.rows() ? A(i, k) : Field()))
                    << ((i != A.cols()-1) ? "\\\\\n" : "\n");
                }
                numeration++;
                out << "\\end{bmatrix}C_{" << numeration << "}\n" ;
            }
        }
        return out.str();
    }
~~~

Далее в модуле для Python, именно результат этой функции будет рендериться как решение.


## Frac | Класс для дроби
Реализация для дробей выглядит следующим образом: (Полная реализация в исходниках в Приложении)

~~~c++
class Frac;
bool operator==(const Frac& a, const Frac& b);
bool operator!=(const Frac& a, const Frac& b);
Frac operator+(const Frac& a, const Frac& b);
Frac operator-(const Frac&a);
Frac operator-(const Frac& a, const Frac& b);
Frac operator*(const Frac& a, const Frac& b);
Frac operator/(const Frac& a, const Frac& b);


using namespace std;
ostream& operator << (ostream& out, const Frac& a);
istream& operator >> (istream& in, Frac& a);

class Frac {
private:
    int _m;
    unsigned int _n;
    
    void simplify(int& m, unsigned& n) {...}
    
public:
    Frac(int m, unsigned n): _m(m), _n(n) {...}
    
    Frac(int m) : _m(m), _n(1) {};
         Frac() : _m(0), _n(1) {};
    
         int m() const { return _m;}
    unsigned n() const { return _n;}
    
    Frac operator  =(const Frac& a) {...}
    Frac operator /=(const Frac& a) {...}
    Frac operator +=(const Frac& a) {...}
    Frac operator -=(const Frac& a) {...}
}
~~~

## Обертка для матрицы на cython:

### Определение:

Подключаем написанный на C++ класс:

```cython
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
```

### Обертка:

Далее класс-обертка, который содержит в себе экземпляр нашего класса Matrix<Frac>, который умеет его красиво отображать и рендерить latex формулы:

```cython
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
```

# Софт и источники

1. jupyter http://jupyter.org
2. cython http://cython.org
3. gcc https://gcc.gnu.org
4. python http://python.org
5. typora https://typora.io
6. Статья о Гауссовом преобразовании https://en.wikipedia.org/wiki/Gaussian_elimination







# Приложения

##Тесты

##Исходный код



