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
x \in X = b' + Kern(A)
$$
где $Kern(A)$ - ядро оператора, связанного с нашей матрицей. (Ядро легко находится из приведенной матрицы)

Для этого реализуем класс на C++ для матрицы над произвольным полем(с помощью шаблонов). И в качестве примера поля возьмем рациональные числа - реализуем класс для дробей(в отличии от float и double, позволяет избежать погрешностей).

В качестве полидрома для испытаний я выбрал Jupter Notebook, т.к. это удачная интерактивная среда, которая умеет отображать tex формулы, что намного интереснее, чем смотреть на матрицу из дробей в обычной консоли. Поэтому так-же далее будет написан класс-обертка на cython для создания Python модуля и дополнительные методы для генерации tex выражений.


# Классы
Так-как достаточно часто придется пробегать по элементам матрицы, определим следующий макрос - синтаксический сахар для for:

~~~c
#define range(i, begin, end) for(size_t i = (begin); i < (end); i++)
~~~
## Matrix | Класс для матрицы

Рассмотрим определение класса **Matrix**:

~~~c++
template <typename Field>
class Matrix {
    vector< vector<Field>* > _M;
    Matrix& self = *this;
    void clear();
    
public:
    Matrix(size_t rows, size_t cols);
    Matrix& operator = (const Matrix& A);
    Matrix& operator = (Matrix&& A);
    Matrix(const Matrix &A);
    Matrix(Matrix&& A);
    ~Matrix();

    // Индексация
    Field& operator() (size_t i, size_t j) const;
    
    // Элементарные преобразования
    void add(size_t to, size_t from, Field k);
    void swap(size_t first, size_t second);
    void mul(size_t row, Field k);
    void div(size_t row, Field k);
    
    // Метод Гаусса-Жордана
    pair< bool, vector<bool> > row_reduce(vector<Field>& b);
    
    // Работа с потоками
    friend istream& operator >> (istream& in, Matrix &M);
    friend ostream& operator << (ostream& out, const Matrix &M);
	
	// Операторы
	friend Matrix operator * (const Matrix& A, const Matrix& B);
    friend Matrix operator + (const Matrix& A, const Matrix& B);
    friend Matrix operator - (const Matrix& A);
    friend Matrix operator - (const Matrix& A, const Matrix& B);
    friend bool operator == (const Matrix& A, const Matrix& B);
    friend bool operator != (const Matrix& A, const Matrix& B);

	// Данный метод создан для 
    void set(size_t i, size_t j, const Field& val);

	string tex_solution(vector<Field>& b);
    string str();

};
~~~

Решение использовать массив из указателей на вектор связаны с частой необходимостью переставлять строки местами. Гораздо проще переставить местами 2 указателя, чем перемещать целые массивы.

Поле **self** используется исключительно для удобства обращения к объекту:

~~~c++
(*this)(i, j); // не очень
self(i, j);    // уже лучше
~~~

И, так как освобождать память нужно будет не только в конструкторе, объявлен **clear()**:

~~~c++
	void clear() {
		for(auto& row: _M) delete row;
	}
~~~

Конструктор по размеру матрицы:

~~~c++
    Matrix(size_t rows, size_t cols) {
        _M.resize(rows);
        for(auto& row: _M) row = new vector<Field>(cols);
    }
~~~

Выделяем память под строки, устанавливая размеры массивов, исходя из параметров.

Далее копирующий оператор присваивания и конструктор:

~~~c++
	Matrix& operator = (const Matrix& A) {
        clear();
        _M.resize(A.rows());
        range(i, 0, A.rows()) _M[i] = new vector<Field>(*(A._M[i]));
        return self;
    }

	Matrix(const Matrix &A) {
        self = A;
    }
~~~

Они необходимы, т.к. мы работаем с динамической памятью, а значит нужно явно описать процесс копирования, т.к. копировать значение указателей это далеко не лучшая идея, а компилятор не знает, что мы там напридумывали. Теперь, имеем возможность передавать по значению и возвращать значения из функций.

Перемещающий опертор присваивания и конструктор:

~~~c++
    Matrix& operator = (Matrix&& A) {
        clear();
        _M = A._M;
        A._M.resize(0);
        return self;
    }

    Matrix(Matrix&& A) {
        _M = A._M;
        A._M.resize(0);
    }
~~~

Воспользуемся прелестями c++11 и определим поведение с rvalue. Зачем копировать объект, когда он вот-вот будет уничтожен? Можно просто забрать его ресурсы и сэкономить время на копировании и удалении.

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

Я не буду приводить здесь все операторы, разберу только умножение матрицы, как самое нетривиальное. Пусть:
$$
A =\begin{bmatrix}
a_{11}&...&a_{1n} \\
... & ... &...\\
a_{m1}&...&a_{mn}
\end{bmatrix},B=\begin{bmatrix}
b_{11}&...&b_{1p} \\
... & ... &...\\
b_{n1}&...&b_{np}
\end{bmatrix},
C=A*B=\begin{bmatrix}
c_{11}&...&c_{1p} \\
... & ... &...\\
c_{m1}&...&c_{mp}
\end{bmatrix}
$$
Тогда элемент матрицы $C$ находится как:
$$
c_{ij}=\sum_{k=1}^n a_{ik}b_{kj}
$$
Естественно, нужно учесть согласование размеров:

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

##### Элементарные преобразования со строками:
Самые необходимые методы, для приведения матрицы к каноническому ступенчатому виду.  

Метод для сложения строк. К заданной строке **to** прибавляется строка **from** домноженная на **k**:

~~~c++
    void add(size_t to, size_t from, Field k) {
        range(j, 0, cols()) self(to, j) += k * self(from, j);
    }
~~~

Метод для перестановки строк. Меняем строки с номерами **first** и **second** местами:

~~~c++
    void swap(size_t first, size_t second) {
        std::swap(_M[first], _M[second]);
    }
~~~

Метод для умножения строки с номером **row** на скаляр **k**:

~~~c++
    void mul(size_t row, Field k) {
        range(j, 0, cols()) self(row,j) = self(row, j) * k;
    }
~~~

Метод для деление строки с номером **row** на скаляр **k**

~~~c++
    void div(size_t row, Field k) {
        range(j, 0, cols()) self(row, j) = self(row,j) / k;
    }
~~~

##### Приведение матрицы медотом Гаусса-Жордана:
Самый важный метод, ради чего все и задумывалось:

Изначально у нас есть расширенная матрица:
$$
\left(\begin{array}{@{} r r r | r @{}}
	a_{11}&...&a_{1n} & b_1\\
	... & ... &...& ...\\
	a_{m1}&...&a_{mn}&b_m
  \end{array}\right)
$$
Пытаемся найти ненулевой элемент $a_{i1}$(далее буду называть опорным элементом и буду считать их кол-во счетчиком **p**) и если не находим, то переходим к следющему столбцу, если же нашли, то меняем первую строку и строку номером **i** местами, после чего делим первую строку на найденное значение и обнуляем оставшиеся значения в столбце, а также увеличиваем счетчик **p**.

> Важное замечание. Далее, я буду подразумевать под ${a_{ij}}$ не изначальное значение, а именно значение в данной позиции в матрице на данном этапе. Т.е. каждое $a_{ij}$ меняет свое значение от шагу к шагу.

Получаем следующую матрицу(в случае, если нашелся ненулевой элемент в первом столбце):
$$
\left(\begin{array}{@{} r r r r | r @{}}
	1&a_{12}&...&a_{1n} & b_1\\
	0&a_{22}&...&a_{2n}& b_2\\
	... & ... &...& ...&...\\
	0&a_{m2}&...&a_{mn}&b_m
  \end{array}\right)
$$
Или остается неизменной, если первый столбец был нулевым(при этом **p**=0):
$$
\left(\begin{array}{@{} r r r r | r @{}}
	0&a_{12}&...&a_{1n} & b_1\\
	0&a_{22}&...&a_{2n}& b_2\\
	... & ... &...& ...&...\\
	0&a_{m2}&...&a_{mn}&b_m
  \end{array}\right)
$$
Далее процесс повторяется, за тем исключением, что ненулевой элемент $a_{ij}$ ищется такой, что $p\leq i$. Если же элемент в столбце не нашлеся, то запоминаем номер столбца и переходим к следющему. Запомнить номера подобных столбцов важно, т.к. подобные столбцы далее войдут в ядро оператора, свзяанного с нашей матрицей, т.е. будут решением однородной системы.

В итоге процесса получаем каноническую ступенчатую матрицу. Она может выглядить как-то так:
$$
\left(\begin{array}{@{} r r r r  r| r @{}}
	0&1&0 & -1 & 0 & b_1'\\
	0&0&1& 3 &0 & b_2'\\
	0 & 0 &0& 0&1 & b_3'\\
	0&0&0&0&0&b'_4
  \end{array}\right)
$$
Если приведенная матрица содержит нулевые строки и этим нулевым строкам соотвествует  ненуловое значение в приведенном столбце $b'$, то, очевидно, что система не совместна, т.к. не существует линейной комбинацией из нулей порождающих ненулевое значение.

В ином же случае решение есть и оно выглядит следующим образом:
$$
x\in X=\{b'+n\,|\,n\in Kern(A)\}
$$
В частности, в случае тривиального ядра, получаем одно решение.

Функция возвращает вердикт о том, есть решение или нет в виде bool и возвращает bool вектор, хранящий информацию о том, какие столбцы войдут а какие нет в ядро. Значения возвращаются в паре:

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

Имеем функцию, которая делит $$m,n$$ на $\text{НОД}(m,n)$. Т.к. имеет смысл хранить дробь в сокращенной форме. 

~~~c++
    void simplify(int& m, unsigned& n) {
        unsigned gcd = GCD((unsigned) std::abs(m), n);
        m /= (int)gcd;
        n /= gcd;
    }
~~~

Ну и раз уж был упомянут НОД, то здесь классический алгоритм Евклида:

~~~c++
unsigned GCD(unsigned u, unsigned v) {
    while ( v != 0) {
        unsigned r = u % v;
        u = v;
        v = r;
    }
    return u;
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

Далее класс-обертка, который содержит в себе экземпляр нашего класса Matrix<Frac> и который умеет приводить питовские типы к типам C++ и наоборот

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



