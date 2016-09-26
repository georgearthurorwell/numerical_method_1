#include <iomanip>
#include <iostream>
#include <cmath>
#include <fstream>
#include "myvector.cpp"

using namespace std;

long double scal(myvector a, myvector b, long double h)
{
    long double s = 0;
    if(a.length == b.length)
    {
        for(int i = 1; i < a.length-1; ++i)
        {
            s += a.m[i] * b.m[i]*h;
            //cout << s << endl;
        }
        s += 0.5 * h * (a.m[0] * b.m[0] + a.m[a.length-1] * b.m[a.length-1]);
        return s;
    }
    else
    {
        cout << "error_scal";
        return s;
    }
}

int main()
{
    // DEFINE N, h, y_n
    int N = 11;
    long double l = 1.0;
    long double h = 1.0 / (N - 1);

    long double c[N] = {};
    for(int i = 0; i < N; ++i)
    {
        c[i] = 1.0; //DEFINE c_i 
    }

    myvector * y = new myvector[N];
    long double coords[N] = {};
    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            coords[j] = c[i]*cos(M_PI*i*j/N); // DEFINE y_i
        }
        y[i] = myvector(N, coords);
    }

    //END DEFINE
    cout.setf(ios::showpos);

    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            cout << fixed << setprecision(3) << scal(y[i],y[j]) << " ";
        }
        cout << endl;
    }
    cout << "helloy blya";
}

// int main()
// {
//     int N;
//     long double a;//коэфициент при решение
//     long double * c = new long double[N];
//     myvector * y = new myvector[N];
//     myvector f;
//     long double * v:

//     for(int i =0; i < N; ++i)
//     {
//         v = new long double[N];
//         for(int j = 0; j < N; ++j)
//         {
//             v[j] = ;//тут формула
//         }
//         y[i].set(N, v);
//     }

//     v = new long double[N];
//     for(int i =0; i < N; ++i)
//     {
//         v[i] = ;//тут формула
//     }
//     f.set(N, v);

//     for(int i = 0; i < N; ++i)
//     {
//         c[i] = scal(f, y[i]);
//     }

// }
