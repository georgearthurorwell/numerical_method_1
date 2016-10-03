#include <iomanip>
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

class myvector{
public:
	myvector()
	{
		length = 0;
		m = 0;
	}
	myvector(int a, long double *b)
	{
		length = a;
		m = b;
	}
	myvector operator=(myvector a)
	{
		if (this != &a)
		{
			delete[] m;
			length = a.length;
			m = new long double[length];
			for (int i = 0; i != length; ++i)
				m[i] = a.m[i];
		}
		return *this;
	}

	void set(int a, long double *b)
	{
		length = a;
		m = b;
	}

	long double scal(myvector a, myvector b)
	{
	    long double s = 0;
	    if(a.length == b.length)
        {
            for(int i = 0; i < a.length; ++i)
            {
                s += a.m[i] * b.m[i];
            }
            return s;
        }
	    else
        {
            cout << "error_scal";
            return s;
        }
	}

	//int nulls()
	//{
	//	int amount = 0;
	//	for (int i = 0; i < this->length - 1; ++i)
	//	{
	//		if (this->m[i] == 0)
	//			++amount;
	//	}
	//	return amount;
	//}
	int length;
	long double * m;
};

myvector operator+(myvector a, myvector b)
{
	if (a.length == b.length)
	{
		long double *m = new long double[a.length];
		myvector c(a.length, m);
		for (int i = 0; i < a.length; ++i)
			c.m[i] = a.m[i] + b.m[i];
		return c;
	}
	else
	{
		cout << "error in +";
		return a;
	}
}

myvector operator-(myvector a, myvector b)
{
	if (a.length == b.length)
	{
		long double *m = new long double[a.length];
		myvector c(a.length, m);
		for (int i = 0; i < a.length; ++i)
			c.m[i] = a.m[i] - b.m[i];
		return c;
	}
	else
	{
		cout << "error in -";
		return a;
	}
}

myvector operator+=(myvector a, myvector b)
{
	if (a.length == b.length)
	{
		for (int i = 0; i < a.length; ++i)
			a.m[i] += b.m[i];
		return a;
	}
	else
	{
		cout << "error in += ";
		return a;
	}
}

myvector operator*(long double d, myvector a)
{
	long double *m = new long double[a.length];
	myvector c(a.length, m);
	for (int i = 0; i < a.length; ++i)
		c.m[i] = a.m[i] * d;
	return c;
}

long double scal(myvector a, myvector b)
{
    long double s = 0;
    if(a.length == b.length)
    {
        for(int i = 0; i < a.length; ++i)
        {
            s += a.m[i] * b.m[i];
            //cout << s << endl;
        }
        return s;
    }
    else
    {
        cout << "error_scal";
        return s;
    }
}

long double max(myvector a)
{
    long double m = a.m[0];
    for(int i = 0; i < a.length; ++i)
    {
        if(a.m[i] > m)
        m = a.m[i];
    }
    return m;
}

int maxcoordinate(myvector a)
{
    int max_ind = 0;
    long double m = a.m[0];
    for(int i = 0; i < a.length; ++i)
    {
        if(a.m[i] > m)
        m = a.m[i];
        max_ind = i;
    }
    return max_ind;
}
