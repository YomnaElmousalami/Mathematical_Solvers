#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <math.h>
using namespace std;

int main()
{
    ///load N, A, b from a file
    int N;
    double **A, *b;
    int i, j, k;
    double sum = 0;

    fstream fin;
    fin.open("problem_4_alt.txt",ios::in); ///you can change it to the file you desire and it will work

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


    ///Display title of program
    ///Allow user to guess values for x in Ax = b
    double *xold, *xnew;
    xold = new double[N];
    xnew = new double[N];
    cout << "Matrix generator: " << endl;
    //cout << "Guess " << N << " values for x: " << endl;
    for (int i = 0; i < N; i++)
    {
        xold[i] = 0.0;
    }

    for(int i = 0; i < N; i++)
    {
        xnew[i] = xold[i];
    }

    double *bpred;
    bpred = new double[N];
    for(int i = 0; i < N; i++)
        {
            bpred[i] = 0;
        }

    double *diff;
    diff = new double[N];

    ///Solve until the Two-Norm is acceptable
    double twoNorm = 10.0;
    int iteration = 0;
    while(twoNorm > 0.0001)
    {
        //
        iteration++;
        for(int i = 0; i < N; i++)
        {
            xnew[i] = b[i];
            for(int j = 0; j < N; j++)
            {
                if(i!=j)
                {
                    if(A[i][j] > 0)
                    {
                       xnew[i] -= A[i][j]*xnew[j];
                      // cout << "2";
                    }
                    else
                    {
                        //if(A[i][j]*xnew[j] > 0)
                        //{
                         //   xnew[i] += A[i][j]*xnew[j];
                         //   cout << "1";
                       // }
                        //else
                        //{
                             xnew[i] += abs(A[i][j])*xnew[j];
                             //cout << "3";
                        //}

                    }
                }
            }
            xnew[i] = xnew[i]/A[i][i];
        }

    for(int i = 0; i < N; i++)
    {
        bpred[i] = 0;
        for(int j = 0; j < N; j++)
        {
            bpred[i] += A[i][j]*xnew[j];
        }
    }

   for(int i = 0; i < N; i++)
   {   diff[i]=0.0;
       diff[i] = b[i] - bpred[i];
   }
   twoNorm=0.0;
   for(int i = 0; i < N; i++)
   {
       twoNorm += diff[i]*diff[i];
   }
   twoNorm = sqrt(twoNorm);
   }


   cout << endl << endl;
   cout << "Iterations: " << iteration << endl << endl;
 ///display matrix data
   // cout << endl;
    //cout << "Matrix size for A: " << N  << "x" << N << endl;
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

    /*cout << endl << endl << "Matrix size for b: " << N  << "x" << 1 << endl;
    cout << "Matrix b:" << endl << endl;
    for(int i = 0; i < N; i++)
    {
        cout << b[i] << endl;
    }

    cout << endl << endl << "Matrix size for x: " << N  << "x" << 1 << endl;
    cout << "Matrix x:" << endl << endl;
    for(int i = 0; i < N; i++)
    {
        cout << xnew[i] << endl;
    }
    cout << endl << endl;
    cout << endl << endl;
    cout << "Difference: b-r:" << endl << endl;
    for(int i = 0; i < N; i++)
    {
        cout << diff[i] << endl;
    }
    cout << endl << endl;
    cout << "Two Norm: " << twoNorm;*/


        int y, nx = sqrt(N);
     fstream fout;
   fout.open("data.txt",ios::out);
   for(int j=0; j<nx; j++)

        {
        for(int i=0; i<nx; i++)
              {   y = j*nx+i;
                 fout<<xnew[y]<<" ";
              }
        fout<<endl;
        }
     fout.close();


    ///Memory deallocation
    delete [] xnew;
    delete [] xold;
    delete [] A;
    delete [] b;
    delete [] diff;

    return 0;
}

