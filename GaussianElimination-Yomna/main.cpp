#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
using namespace std;

int main()
{
    ///load N, A, b from a file
    ///create duplicates of A and b variables (will be used later)
    int N;
    double **A, *b;

    ///temporaries
    double **dupA, *dupb;

    fstream fin;
    fin.open("problem_1_alt.txt",ios::in); ///you can change it to the file you desire and it will work

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

    ///temporaries
    dupA = new double*[N];
    for(int i = 0; i < N; i++)
    {
        dupA[i] = new double[N];
    }
    dupb = new double[N];
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            dupA[i][j] = A[i][j];
        }
    }
    for(int i = 0; i < N; i++)
    {
        dupb[i] = b[i];
    }

    fin.close();

    ///Declare x for: Ax = b
     double *x;
     x = new double[N];

    ///Solve for x when given N, A, b
    ///Step 1: Find the largest absolute value in column 0 in matrix A
    ///and swap target row T with row 0
    double largestVal = abs(A[0][0]), temp = 0, temp_b = 0;
    int largestRow = 0;

    for(int i = 1; i < N; i++)
    {
        if(largestVal < abs(A[i][0]))
        {
            largestVal = abs(A[i][0]);
            largestRow = i;
        }
    }

    if(largestRow != 0)
    {
        for(int j = 0; j < N; j++)
        {
            temp = A[0][j];
            A[0][j] = A[largestRow][j];
            A[largestRow][j] = temp;
        }

            temp_b = b[0];
            b[0] = b[largestRow];
            b[largestRow] = temp_b;
    }

    ///Step 2: Take the value in row 0, store it, and divide row 0 by that value
    double dividend = A[0][0];
    for(int j = 0; j < N; j++)
    {
        A[0][j] = A[0][j]/dividend;
    }
    if(b[0]==0)
    {
        b[0] = 0;
    }
    else
    {
        b[0] = b[0]/dividend;
    }

    ///Step 3: Use the new value in A[0][0] to eliminate all values bellow
     for(int t =1; t<N; t++)
       {
          double wt = A[t][0];
          for(int c = 0; c < N; c++)
          {
            A[t][c]= A[t][c] - (wt*A[0][c]);
          }
            b[t]=b[t]-(wt*b[0]);
       }


    ///Repeat Steps 2 and 3 for the rest of the rows and columns (need help with this)
    int t = 2;
    while(t != N)
    {
        int temporary = t - 1;
        double div = A[temporary][temporary];
        for(int i = 0; i < N; i++)
        {
          A[temporary][i+1] = A[temporary][i+1]/div;
        }
        b[temporary] = b[temporary]/div;

        for(int target = t; target < N; target++)
        {
            double wte = A[target][temporary];
            for(int c = temporary; c < N; c++)
            {
                A[target][c] = A[target][c] - (wte*A[temporary][c]);
            }
            b[target] = b[target] - (wte*b[temporary]);
        }
        t++;
    }

    ///Step 4: Find all solutions by back-solving
    double div = A[N-1][N-1];
    A[N-1][N-1] = A[N-1][N-1]/div;
    b[N-1] = b[N-1]/div;
    x[N-1] = b[N-1];

    int icounter = 0;
    for(int i = N-2; i > -1; i--)
    {
        x[i] = b[i];
        for(int j = 0; j < N; j++)
        {
            if(A[i][j] == 1 && icounter == 0)
            {
              icounter++;
            }
            else if(icounter == 1)
            {
                if(A[i][j]*x[j] > 0)
                {
                    x[i] = x[i] - A[i][j]*x[j];
                }
                else
                {
                    x[i] = x[i] + abs(A[i][j]*x[j]);
                }
            }
        }
        icounter = 0;
    }

    ///Step 5: Find the two-norm
    ///compute your answer for A*x:
    double *r;
    r = new double[N];
    for(int i = 0; i < N; i++)
    {
        r[i] = 0;
        for(int j = 0; j < N; j++)
        {
            r[i] += dupA[i][j]*x[j]; ///ask if dupA or A
        }
    }
    ///Find the difference: b-r = diff
    double *diff, *squared;
    double squaredSum = 0, twoNorm = 0;
    diff = new double[N];
    for(int i = 0; i < N; i++)
    {
        diff[i] = (dupb[i]) - (r[i]); ///ask if dupB or b
    }

    squared = new double[N];
    for(int i = 0; i < N; i++)
    {
        squared[i] = diff[i]*diff[i];
        squaredSum += squared[i];
    }
    twoNorm = sqrt(squaredSum);

///display output to screen
//cout << "Matrix generator: " << endl;
 //cout << "Matrix size for A: " << N  << "x" << N << endl;
  //  cout << "Matrix A:" << endl << endl;
     fstream out;
   out.open("del.txt",ios::out);
    for(int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            out << "A[ " << i << " ]" << "[ " << j << " ] = " << dupA[i][j] << " " << endl;
        }
        //out << endl;
    }
    out.close();
    cout << N << endl;

    /*for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            cout << "A[ " << i << " ]" << "[ " << j << " ]" << dupA[i][j] << endl;
        }
    }*/

   /* cout << endl << "Matrix size for b:" << N  << "x" << 1 << endl;
    cout << "Matrix b:" << endl << endl;
    for(int i = 0; i < N; i++)
    {
        cout << dupb[i] << endl;
    }
    cout << endl << endl << "Matrix size for x: " << N  << "x" << 1 << endl;
    cout << "Matrix x:" << endl << endl;
    for(int i = 0; i < N; i++)
    {
        if(abs(x[i]) < 0.00001)
        {
            x[i] = 0;
        }
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

       int y, nx = sqrt(N);
     fstream fout;
   fout.open("data.txt",ios::out);
   for(int j=0; j<nx; j++)

        {
        for(int i=0; i<nx; i++)
              {   y = j*nx+i;
                 fout<<x[y]<<" ";
              }
        fout<<endl;
        }
     fout.close();

    ///Memory deallocation
     delete [] A;
     delete [] b;
     delete [] dupA;
     delete [] dupb;
     delete [] x;
     delete []r;
     delete [] diff;
     delete [] squared;
    return 0;
}
