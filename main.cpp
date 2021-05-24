#include "classes.h"

double f(double x,double y){
  return sin(x)*sin(y)+d0;
}
void solved(matrix &m){ //analytical solution 
  int dim = m.get_dim();
  for(int i=0;i<dim;i++){
        for(int j=0;j<dim;j++){
            double h = m.get_h();
            double x = h*i;
            double y = h*j;
            
            m.input_answer(i,j,f(x,y));
        }
    }
}

int main()
{
    int N = 400; //array size
    double h = PI / (N - 1); //cell size
    double time_ref = 0;
    cout << "2D Array Size:" << N << endl;
    for (int i = 1; i < 9; i++) {
        matrix pot(N, h);

        matrix dens(N, h);
        dens.init_sin_dens();

        matrix ans(N, h);
        solved(ans);

        struct timeval start, end;
        gettimeofday(&start, NULL);



        pot.SOR_smoothing(dens, ans, 1.9, 1e-6, i);

        gettimeofday(&end, NULL);
        double time_use = (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0;
        if (i == 1) {
            cout << "number of threads:" << i << "  ,performance:" << 1.0 << endl;
            time_ref = time_use;
        }
        else cout << "number of threads:" << i << "  ,performance:" << time_ref/time_use << endl;
    }
}
  
