#include <iostream>
#include "vMatrix.hpp"
#include "Frac.hpp"
#include <iomanip>


int main(int argc, const char * argv[]) {
    
    Matrix<Frac> M(2,5);
    M(0, 0) = 2; M(0,1) = 4; M(0,2) = 6; M(0,3) = 0; M(0,4) = 0;
    M(1, 0) = 4; M(1,1) = 10;M(1,2) = 0; M(1,3) = 0; M(1,4) = 1;
    
    vector<Frac> b = {7, 2};
    
//    cout << M << endl;
    
    //   print_solution(M, b, cout);
    
    //   cout << M.tex_solution(b) << endl;
    
    M.row_reduce(b);
    
    cout << M << endl;
    return 0;
    
}
