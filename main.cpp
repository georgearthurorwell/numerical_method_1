#include <iomanip>
#include <iostream>
#include <cmath>
#include <fstream>
#include "myvector.cpp"

using namespace std;

long double scal(myvector , myvector , long double );

const long double eps = 0.000000000001;

long double scal(myvector a, myvector b, long double h)
{
    long double s = 0;
    if(a.length == b.length)
    {
        for(int i = 1; i < a.length-1; ++i)
        {
            s += a.m[i] * b.m[i]*h;
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
    // DEFINE N, h, y_n, f
    int N = 1001;
    //long double l = 1.0000000000000000;
    long double h = 1.0 / (N - 1);

    long double *c = new long double[N];
    for(int i = 0; i < N; ++i)
    {
        c[i] = pow(2.0,0.5); //DEFINE c_i
    }

    myvector * y = new myvector[N];
    long double *coords = new long double[N];
    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            coords[j] = c[i]*sin(M_PI*i*j*h); // DEFINE y_i
        }
        y[i] = myvector(N, coords);
    }

    myvector f = myvector();
    for(int i = 0; i < N; ++i)
    {
        //coords[i] = cos(2*M_PI*i*h) - 1; // DEFINE f
        coords[i] = cos(2*M_PI*i*h) - 1;
    }
    f = myvector(N, coords);

    //CHECK CORRECT f
    if(abs(f.m[0]) > eps || abs(f.m[N-1]) > eps)
    {
        cout << "ERROR. f does not satisfy border conditions" << f.m[0] << ' ' << f.m[N-1] << endl;
        return -1;
    }

    //END DEFINE

    long double *d = new long double[N];
    for(int i = 0; i < N; ++i)
    {
        d[i] = scal(f,y[i],h);
        //cout << d[i] << ' ';
    }
    
    long double delta = 0;
    int a,b;
    ofstream fout1("delta.txt");
    fout1.setf(ios::showpos);
    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
	    if(i != j){
	        if(abs(delta) < scal(y[i], y[j], h)){
	            delta = scal(y[i],y[j],h);
		    a = i;
		    b = j;
	        }
	    }
	    //delta = scal(y[i], x[i], h);
	    //if(delta ) 
            fout1 << fixed << setprecision(10) << scal(y[i],y[j], h) << " ";
        }
        fout1 << endl;
    }
    fout1.close();
    
    cout << "i=" << a << " j=" << b << " delta=" << delta;    

    for(int i = 0; i < N; ++i)
    {
        coords[i] = 0;
    }
    myvector f_1 = myvector(N, coords);
    for(int i = 0; i < N; ++i)
    {
        y[i] = d[i] * y[i];
        f_1 += y[i];
    }

    ofstream fout("output.txt");
    for(int i = 0; i < N; ++i)
    {
        fout << h*i << ' ' << f_1.m[i] << endl;
    }
    fout.close();
}

// int main()
// {
//     int N;
//     long double a;//���������� ��� �������
//     long double * c = new long double[N];
//     myvector * y = new myvector[N];
//     myvector f;
//     long double * v:

//     for(int i =0; i < N; ++i)
//     {
//         v = new long double[N];
//         for(int j = 0; j < N; ++j)
//         {
//             v[j] = ;//��� �������
//         }
//         y[i].set(N, v);
//     }

//     v = new long double[N];
//     for(int i =0; i < N; ++i)
//     {
//         v[i] = ;//��� �������
//     }
//     f.set(N, v);

//     for(int i = 0; i < N; ++i)
//     {
//         c[i] = scal(f, y[i]);
//     }

// }
