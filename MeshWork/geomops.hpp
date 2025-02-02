#ifndef FILE_GEOMOPS
#define FILE_GEOMOPS

/* *************************************************************************/
/* copy from netgen  and modify it                                         */
/* File:   geomop.hpp                                                 */
/* Function:   overloaded operator of geometry fundamental type            */ 
/* Date:   09. jan. 2016                                                   */
/* Class:  Point - Vector operations     Matrix - Vector operations        */
/* *************************************************************************/

namespace meshwork
{

  /*

  Point - Vector operations

  */
	//向量加法
  template <int D>
  inline Vec<D> operator+ (const Vec<D> & a, const Vec<D> & b)
  {
    Vec<D> res;
    for (int i = 0; i < D; i++)
      res(i) = a(i) + b(i);
    return res;
  }
  //点加向量
  template <int D>
  inline Point<D> operator+ (const Point<D> & a, const Vec<D> & b)
  {
    Point<D> res;
    for (int i = 0; i < D; i++)
      res(i) = a(i) + b(i);
    return res;
  }
  //两点相减得向量
  template <int D>
  inline Vec<D> operator- (const Point<D> & a, const Point<D> & b)
  {
    Vec<D> res;
    for (int i = 0; i < D; i++)
      res(i) = a(i) - b(i);
    return res;
  }
  //点减向量得点
  template <int D>
  inline Point<D> operator- (const Point<D> & a, const Vec<D> & b)
  {
    Point<D> res;
    for (int i = 0; i < D; i++)
      res(i) = a(i) - b(i);
    return res;
  }
  //向量相减
  template <int D>
  inline Vec<D> operator- (const Vec<D> & a, const Vec<D> & b)
  {
    Vec<D> res;
    for (int i = 0; i < D; i++)
      res(i) = a(i) - b(i);
    return res;
  }
  //向量数乘
  template <int D>
  inline Vec<D> operator* (double s, const Vec<D> & b)
  {
    Vec<D> res;
    for (int i = 0; i < D; i++)
      res(i) = s * b(i);
    return res;
  }
  //向量点乘
  template <int D>
  inline double operator* (const Vec<D> & a, const Vec<D> & b)
  {
    double sum = 0;
    for (int i = 0; i < D; i++)
      sum += a(i) * b(i);
    return sum;
  }
  //向量取负
  template <int D>
  inline Vec<D> operator- (const Vec<D> & b)
  {
    Vec<D> res;
    for (int i = 0; i < D; i++)
      res(i) = -b(i);
    return res;
  }
  //点加向量得点
  template <int D>
  inline Point<D> & operator+= (Point<D> & a, const Vec<D> & b)
  {
    for (int i = 0; i < D; i++)
      a(i) += b(i);
    return a;
  }
  //向量相加
  template <int D>
  inline Vec<D> & operator+= (Vec<D> & a, const Vec<D> & b)
  {
    for (int i = 0; i < D; i++)
      a(i) += b(i);
    return a;
  }
  //点减向量
  template <int D>
  inline Point<D> & operator-= (Point<D> & a, const Vec<D> & b)
  {
    for (int i = 0; i < D; i++)
      a(i) -= b(i);
    return a;
  }
  //向量相减
  template <int D>
  inline Vec<D> & operator-= (Vec<D> & a, const Vec<D> & b)
  {
    for (int i = 0; i < D; i++)
      a(i) -= b(i);
    return a;
  }
  //向量数乘
  template <int D>
  inline Vec<D> & operator*= (Vec<D> & a, double s)
  {
    for (int i = 0; i < D; i++)
      a(i) *= s;
    return a;
  }

  //向量数除
  template <int D>
  inline Vec<D> & operator/= (Vec<D> & a, double s)
  {
    for (int i = 0; i < D; i++)
      a(i) /= s;
    return a;
  }

    /*

	Matrix - Vector operations

  */
  //矩阵乘向量
    template <int H, int W>
    inline Vec<H> operator* (const Mat<H,W> & m, const Vec<W> & v)
	{
		Vec<H> res;
		for (int i = 0; i < H; i++)
		{
			res(i) = 0;
			for (int j = 0; j < W; j++)
				res(i) += m(i,j) * v(j);
		}
		return res;
	}
  

  // thanks to VC60 partial template specialization features !!!
	//矩阵乘向量
  inline Vec<2> operator* (const Mat<2,2> & m, const Vec<2> & v)
  {
	  Vec<2> res;
	  for (int i = 0; i < 2; i++)
	  {
		  res(i) = 0;
		  for (int j = 0; j < 2; j++)
			  res(i) += m(i,j) * v(j);
	  }
	  return res;
  }
  //矩阵乘向量
  inline Vec<2> operator* (const Mat<2,3> & m, const Vec<3> & v)
  {
	  Vec<2> res;
	  for (int i = 0; i < 2; i++)
	  {
		  res(i) = 0;
		  for (int j = 0; j < 3; j++)
			  res(i) += m(i,j) * v(j);
	  }
	  return res;
  }
  //矩阵乘向量
  inline Vec<3> operator* (const Mat<3,2> & m, const Vec<2> & v)
  {
	  Vec<3> res;
	  for (int i = 0; i < 3; i++)
	  {
		  res(i) = 0;
		  for (int j = 0; j < 2; j++)
			  res(i) += m(i,j) * v(j);
	  }
	  return res;
  }

  //矩阵乘向量
  inline Vec<3> operator* (const Mat<3,3> & m, const Vec<3> & v)
  {
	  Vec<3> res;
	  for (int i = 0; i < 3; i++)
	  {
		  res(i) = 0;
		  for (int j = 0; j < 3; j++)
			  res(i) += m(i,j) * v(j);
	  }
	  return res;
  }

  /*
    template <int H1, int W1, int H2, int W2>
    inline Mat<H1,W2> operator* (const Mat<H1,W1> & a, const Mat<H2,W2> & b)
    {
    Mat<H1,W2> m;
    for (int i = 0; i < H1; i++)
    for (int j = 0; j < W2; j++)
    {
    double sum = 0;
    for (int k = 0; k < W1; k++)
    sum += a(i,k) * b(k, j);
    m(i,j) = sum; 
    }
    return m;
    }
  */
  //矩阵乘矩阵
  inline Mat<2,2> operator* (const Mat<2,2> & a, const Mat<2,2> & b)
  {
    Mat<2,2> m;
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++)
	{
	  double sum = 0;
	  for (int k = 0; k < 2; k++)
	    sum += a(i,k) * b(k, j);
	  m(i,j) = sum; 
	}
    return m;
  }
  //矩阵乘矩阵
  inline Mat<2,2> operator* (const Mat<2,3> & a, const Mat<3,2> & b)
  {
    Mat<2,2> m;
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++)
	{
	  double sum = 0;
	  for (int k = 0; k < 3; k++)
	    sum += a(i,k) * b(k, j);
	  m(i,j) = sum; 
	}
    return m;
  }
  //矩阵乘矩阵
  inline Mat<3,2> operator* (const Mat<3,2> & a, const Mat<2,2> & b)
  {
    Mat<3,2> m;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 2; j++)
	{
	  double sum = 0;
	  for (int k = 0; k < 2; k++)
	    sum += a(i,k) * b(k, j);
	  m(i,j) = sum; 
	}
    return m;
  }

  //矩阵乘矩阵
  inline Mat<2,3> operator* (const Mat<2,2> & a, const Mat<2,3> & b)
  {
    Mat<2,3> m;
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 3; j++)
	{
	  double sum = 0;
	  for (int k = 0; k < 2; k++)
	    sum += a(i,k) * b(k, j);
	  m(i,j) = sum; 
	}
    return m;
  }

  //矩阵乘矩阵
  inline Mat<3,3> operator* (const Mat<3,3> & a, const Mat<3,3> & b)
  {
    Mat<3,3> m;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	{
	  double sum = 0;
	  for (int k = 0; k < 3; k++)
	    sum += a(i,k) * b(k, j);
	  m(i,j) = sum; 
	}
    return m;
  }

  //矩阵的逆
  template <int H, int W>
  inline Mat<W,H> Trans (const Mat<H,W> & m)
  {
    Mat<W,H> res;
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
	res(j,i) = m(i,j);
    return res;
  }

  template <int D>
  inline ostream & operator<< (ostream & ost, const Vec<D> & a)
  {
    ost << "(";
    for (int i = 0; i < D-1; i++)
      ost << a(i) << ", ";
    ost << a(D-1) << ")";
    return ost;
  }

  template <int D>
  inline ostream & operator<< (ostream & ost, const Point<D> & a)
  {
    ost << "(";
    for (int i = 0; i < D-1; i++)
      ost << a(i) << ", ";
    ost << a(D-1) << ")";
    return ost;
  }

  template <int D>
  inline ostream & operator<< (ostream & ost, const Box<D> & b)
  {
    ost << b.PMin() << " - " << b.PMax();
    return ost;
  }

  template <int H, int W>
  inline ostream & operator<< (ostream & ost, const Mat<H,W> & m)
  {
	  ost << "(";
	  for (int i = 0; i < H; i++)
	  {
		  for (int j = 0; j < W; j++)
			  ost << m(i,j) << "   ";
		  ost << endl;
	  }
	  return ost;
  }


}

#endif
