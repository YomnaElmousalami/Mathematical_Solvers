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
    fin.open("problem_3_alt.txt",ios::in); ///you can change it to the file you desire and it will work

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


    ///The sum of two matrices are D (Diagonal) and C (Everything Else)
    ///Create matrix D
    double**D;
    D = new double*[N];
    for(int i=0; i<N; i++)
    {
        D[i]=new double[N];
    }

    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            if(i == j)
            {
                D[i][j] = A[i][j];
            }
            else
            {
                D[i][j] = 0;
            }
        }
    }

    ///Create matrix C
    double**C;
    C = new double*[N];
    for(int i=0; i<N; i++)
    {
        C[i]=new double[N];
    }

    int columnController = N-1;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            if(j!=i)
            {
                C[i][j] = A[i][j];
            }
            else
            {
                C[i][j] = 0;
            }
        }
        columnController--;
    }

    ///Create D inverse: D^-1
    double *dInv;
    dInv = new double[N];

    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            if(D[i][j] != 0)
            {
                dInv[i] = 1/D[i][j];
            }
        }
    }

    ///Display title of program
    ///Allow user to guess values for x in Ax = b
    double *xold, *xnew;
    xold = new double[N];
    xnew = new double[N];
    cout << "Matrix generator: " << endl;
    cout << "Guess " << N << " values for x: " << endl;
   /* for (int i = 0; i < N; i++)
    {
        cin >> xold[i];
    }*/

    for (int i = 0; i < N; i++)
    {
        xold[i] = 0.0;
    }

///Solve for xnew
///Create a vector for bpred, this will be used for while loop
double *bpred;
bpred = new double[N];
for(int i = 0; i < N; i++)
{
   bpred[i] = 0;
}


int numCorrect = 0;
double twoNorm = 0;
double *r;
r = new double[N];
double *diff;
diff = new double[N];

for(int i = 0; i < N; i++)
{
   if(b[i] == bpred[i])
   {
       numCorrect++;
   }
}

int iteration = 0;
twoNorm =10.0;
while (twoNorm>.00001)
{
    //
    iteration++;
    for(int i = 0; i < N; i++)
    {
        xnew[i] = 0;
        for(int j = 0; j < N; j++)
        {
            xnew[i] += C[i][j]*xold[j];
        }
    }

    for(int i = 0; i < N; i++)
    {
        xnew[i] = b[i] - xnew[i];
    }

    for(int i = 0; i < N; i++)
    {
        xnew[i] = dInv[i]*xnew[i];
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
        {
            xold[i] = 0;
            xold[i] = xnew[i];
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


///test
cout << endl << endl;
cout << "Iteration: " << iteration << endl;
cout << endl << endl;

   ///display matrix data
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
    /*cout << "Matrix size for C: " << N  << "x" << N << endl;
    cout << "Matrix C:" << endl << endl;
    for(int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            cout << left << setw(8) << setprecision(4) << C[i][j] << setw(8);
        }
        cout << endl;
    }
    cout << "Matrix size for D: " << N  << "x" << N << endl;
    cout << "Matrix D:" << endl << endl;
    for(int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            cout << left << setw(8) << setprecision(4) << D[i][j] << setw(8);
        }
        cout << endl;
    }
    cout << "Matrix size for D^-1: " << N  << "x" << 1 << endl;
    cout << "Matrix D^-1:" << endl << endl;
    for(int i = 0; i < N; i++)
    {
        cout << dInv[i] << endl;
    }
    cout << endl << endl << "Matrix size for x: " << N  << "x" << 1 << endl;
    cout << "Matrix x:" << endl << endl;
    for(int i = 0; i < N; i++)
    {
        cout << xnew[i] << endl;
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
    delete [] r;
    delete [] diff;
    delete [] xnew;
    delete [] xold;
    delete [] A;
    delete [] b;
    delete [] C;
    delete [] D;
    delete [] dInv;


return 0;
}
