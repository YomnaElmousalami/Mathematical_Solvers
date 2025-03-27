#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
using namespace std;

int main()
{
///load N, A, b from a file
    int N;
    double **A, *b;
    int i, j, k;
    double sum = 0;

    fstream fin;
    fin.open("problem_2_alt.txt",ios::in); ///you can change it to the file you desire and it will work

    ///Initializing A, b, and N
    fin>>N;
    A = new double*[N];
    for(int i=0; i<N; i++)
    {
        A[i]=new double[N];
    }

    b = new double[N];
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
         {
             fin>>A[i][j];
         }
    }

    for(int i=0; i<N; i++)
    {
        fin>>b[i];
    }

    fin.close();

    ///Declare x for: Ax = b
    ///Declare L, U for A = LU
     double *x, **L, **U;
     x = new double[N];
     L = new double *[N];
     U = new double *[N];

    ///Find L and U by factoring A
    for(int i = 0; i<N; i++)
    {
        L[i] = new double[N];
        U[i] = new double[N];
    }

    for (i = 0; i < N; i++) {
        U[i][i] = 1;
    }

    for (j = 0; j < N; j++) {

        for (i = j; i < N; i++) {
            sum = 0;
            for (k = 0; k < j; k++) {
                sum = sum + L[i][k] * U[k][j];
            }
            L[i][j] = A[i][j] - sum;
            //cout<<"L["<<i<<"]["<<j<<"] = "<<L[i][j]<<endl;
        }

        for (i = j; i < N; i++) {
            sum = 0;
            for(k = 0; k < j; k++) {
                sum = sum + L[j][k] * U[k][i];
                }

            U[j][i] = (A[j][i] - sum) / L[j][j];
                //cout<<"U["<<j<<"]["<<i<<"] = "<<U[j][i]<<endl;

        }
    }

    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            if(abs(L[i][j]) < 0.00001)
            {
                L[i][j] = 0;
            }
        }
    }

    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            if(abs(U[i][j]) < 0.00001)
            {
                U[i][j] = 0;
            }
        }
    }

    ///Since Ux = y
    ///Find y for Ly = b

    ///declare y
     double *y;
     y = new double[N];

     ///solve for y in Ly = b
     for(int i = 0; i < N; i++)
     {
         y[i] = b[i];
         for(int j = 0; j < N; j++)
         {
             if(L[i][j] != 0)
             {
                if(j != i)
                {
                    if(L[i][j] > 0)
                    {
                       y[i] -= L[i][j]*y[j];
                    }
                    else
                    {
                        y[i] += abs(L[i][j])*y[j];
                    }

                }
             }
         }
         y[i] = y[i]/L[i][i];
     }

     ///Solve for x if y = Ux
    ///Find all solutions by back-solving
    double div = U[N-1][N-1];
    U[N-1][N-1] = U[N-1][N-1]/div;
    x[N-1] = y[N-1]/div;

    int icounter = 0;
    for(int i = N-2; i > -1; i--)
    {
        x[i] = y[i];
        for(int j = 0; j < N; j++)
        {
            if(U[i][j] == 1 && icounter == 0)
            {
              icounter++;
            }
            else if(icounter == 1)
            {
                if(U[i][j]*x[j] > 0)
                {
                    x[i] = x[i] - U[i][j]*x[j];
                }
                else
                {
                    x[i] = x[i] + abs(U[i][j]*x[j]);
                }
            }
        }
        icounter = 0;
    }

    ///Find two norm
    ///Find r by multipling A*x
    double *r;
    r = new double[N];
    for(int i = 0; i < N; i++)
    {
        r[i] = 0;
        for(int j = 0; j < N; j++)
        {
            r[i] += A[i][j]*x[j]; ///ask if dupA or A
        }
    }

   ///Find the difference of b and r
   double *diff;
   diff = new double[N];
   for(int i = 0; i < N; i++)
   {
       diff[i] = b[i] - r[i];
   }

   double twoNorm = 0;
   for(int i = 0; i < N; i++)
   {
       twoNorm += diff[i]*diff[i];
   }
   twoNorm = sqrt(twoNorm);

    ///display matrix data
   // cout << "Matrix generator" << endl << endl;
   // cout << "Matrix size for A: " << N  << "x" << N << endl;
    //cout << "Matrix A:" << endl << endl;
     fstream f;
   f.open("del.txt",ios::out);
    for(int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            f << "A[ " << i << " ]" << "[ " << j << " ] = " << A[i][j] << " " << endl;
        }
    }
    f.close();

    /* cout << endl << "Matrix size for L and U: " << N  << "x" << N << endl;
     cout << "Matrix L:" << endl << endl;
     for(int i = 0; i < N; i++)
     {
        for (int j = 0; j < N; j++)
        {
            cout << left << setw(8) << setprecision(4) << L[i][j] << setw(8);
        }
        cout << endl;
     }
     cout << endl << "Matrix U:" << endl << endl;
     for(int i = 0; i < N; i++)
     {
        for (int j = 0; j < N; j++)
        {
            cout << left << setw(8) << setprecision(4) << U[i][j] << setw(8);
        }
        cout << endl;
     }
    cout << endl << "Matrix size for y:" << N  << "x" << 1 << endl;
    cout << "Matrix y:" << endl << endl;
    for(int i = 0; i < N; i++)
    {
        cout << y[i] << endl;
    }
    cout << endl << endl << "Matrix size for x: " << N  << "x" << 1 << endl;
    cout << "Matrix x:" << endl << endl;
    for(int i = 0; i < N; i++)
    {
        cout << x[i] << endl;
    }
    cout << endl << endl;
    cout << "Matrix size for r: " << N << "x" << 1 << endl;
    cout << "Matrix r:" << endl << endl;
    for(int i = 0; i < N; i++)
    {
        cout << r[i] << endl;
    }
    cout << endl << endl;
    cout << "Difference: b-r:" << endl << endl;
    for(int i = 0; i < N; i++)
    {
        cout << diff[i] << endl;
    }
    cout << endl << endl;
    cout << "Two Norm: " << twoNorm;*/

        int z, nx = sqrt(N);
     fstream fout;
   fout.open("data.txt",ios::out);
   for(int j=0; j<nx; j++)

        {
        for(int i=0; i<nx; i++)
              {   z = j*nx+i;
                 fout<<x[z]<<" ";
              }
        fout<<endl;
        }
     fout.close();

    ///Memory deallocation
    delete [] r;
    delete [] diff;
    delete [] U;
    delete []L;
    delete [] y;
    delete [] x;
    delete [] A;
    delete [] b;
    return 0;
}
