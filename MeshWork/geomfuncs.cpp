#include"externalcommon.h"
#include"internalcommon.h"

namespace meshwork
{
	//3*3矩阵的逆
	void CalcInverse (const Mat<3,3> & m, Mat<3,3> & inv)
	{
		double det = Det (m);
		if (det == 0) 
		{
			inv = 0;
			return;
		}

		double idet = 1.0 / det;
		inv(0,0) =  idet * (m(1,1) * m(2,2) - m(1,2) * m(2,1));
		inv(1,0) = -idet * (m(1,0) * m(2,2) - m(1,2) * m(2,0));
		inv(2,0) =  idet * (m(1,0) * m(2,1) - m(1,1) * m(2,0));

		inv(0,1) = -idet * (m(0,1) * m(2,2) - m(0,2) * m(2,1));
		inv(1,1) =  idet * (m(0,0) * m(2,2) - m(0,2) * m(2,0));
		inv(2,1) = -idet * (m(0,0) * m(2,1) - m(0,1) * m(2,0));

		inv(0,2) =  idet * (m(0,1) * m(1,2) - m(0,2) * m(1,1));
		inv(1,2) = -idet * (m(0,0) * m(1,2) - m(0,2) * m(1,0));
		inv(2,2) =  idet * (m(0,0) * m(1,1) - m(0,1) * m(1,0));
	}

	/*
	// template <>
	void CalcInverse (const Mat<2,3> & m, Mat<3,2> & inv)
	{
	Mat<2,2> a = m * Trans (m);
	Mat<2,2> ainv;
	CalcInverse (a, ainv);
	inv = Trans (m) * ainv;
	}
	*/

	//二阶行列式的值
	double Det (const Mat<2,2> & m) 
	{
		return  m(0,0) * m(1,1) - m(0,1) * m(1,0);
	}
	//三阶行列式的值
	double Det (const Mat<3,3> & m) 
	{
		return 
			m(0,0) * m(1,1) * m(2,2)
			+ m(1,0) * m(2,1) * m(0,2)
			+ m(2,0) * m(0,1) * m(1,2)
			- m(0,0) * m(2,1) * m(1,2)
			- m(1,0) * m(0,1) * m(2,2)
			- m(2,0) * m(1,1) * m(0,2);
	}

	//三阶行列式的特征值
	void EigenValues (const Mat<3,3> & m, Vec<3> & ev)
	{
		const double pi = 3.141592;
		double a, b, c, d;
		double p, q;
		double arg;

		a = -1.;
		b = m(0,0) + m(1,1) + m(2,2);
		c = -( m(0,0)*m(2,2) + m(1,1)*m(2,2) + m(0,0)*m(1,1) - pow(m(0,1),2) - pow(m(0,2),2) - pow(m(1,2),2) );
		d = Det (m);
		p = 3.*a*c - pow(b,2);
		q = 27.*pow(a,2)*d - 9.*a*b*c + 2.*pow(b,2)*b;

		arg = acos((-q/2)/sqrt(-(p*p*p)));

		ev(0) = (2. * sqrt(-p) * cos(arg/3.) - b) / 3.*a;
		ev(1) = (-2. * sqrt(-p) * cos(arg/3.+pi/3) - b) / 3.*a;
		ev(2) = (-2. * sqrt(-p) * cos(arg/3.-pi/3)- b) / 3.*a;
	}


}
