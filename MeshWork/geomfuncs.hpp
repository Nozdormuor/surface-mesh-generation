#ifndef FILE_GEOMFUNCS
#define FILE_GEOMFUNCS

/* *************************************************************************/
/* copy from netgen  and modify it                                         */
/* File:   geomobjects.hpp                                                 */
/* Function:  functions of  geometry fundamental type                     */ 
/* Date:   09. jan. 2016                                                   */
/* Class: Point Vec Mat Box BoxSphere                                      */
/* *************************************************************************/



namespace meshwork 
{
	//�����ľ���ֵ
  template <int D>
  inline double Abs (const Vec<D> & v)
  {
    double sum = 0;
    for (int i = 0; i < D; i++)
      sum += v(i) * v(i);
    return sqrt (sum);
  }
  //��������ֵ��ƽ��
  template <int D>
  inline double Abs2 (const Vec<D> & v)
  {
    double sum = 0;
    for (int i = 0; i < D; i++)
      sum += v(i) * v(i);
    return sum;
  }
  //����ľ���
  template <int D>
  inline double Dist (const Point<D> & a, const Point<D> & b)
  {
    return Abs (a-b);
  }
  //��������ƽ��
  template <int D>
  inline double Dist2 (const Point<D> & a, const Point<D> & b)
  {
    return Abs2 (a-b);
  }

  //��������ģ��������أ�
  template <int D>
  inline Point<D> Center (const Point<D> & a, const Point<D> & b)
  {
    Point<D> res;
    for (int i = 0; i < D; i++)
      res(i) = 0.5 * (a(i) + b(i));
    return res;
  }
  //��������ģ��������أ�
  template <int D>
  inline Point<D> Center (const Point<D> & a, const Point<D> & b, const Point<D> & c)
  {
    Point<D> res;
    for (int i = 0; i < D; i++)
      res(i) = (1.0/3.0) * (a(i) + b(i) + c(i));
    return res;
  }
  //�ĵ�����ģ��������أ�
  template <int D>
  inline Point<D> Center (const Point<D> & a, const Point<D> & b, const Point<D> & c, const Point<D> & d)
  {
    Point<D> res;
    for (int i = 0; i < D; i++)
      res(i) = (1.0/4.0) * (a(i) + b(i) + c(i) + d(i));
    return res;
  }
  //�����Ĳ��
  inline Vec<3> Cross (Vec<3> v1, Vec<3> v2)
  {
    return Vec<3> 
      ( v1(1) * v2(2) - v1(2) * v2(1),
	v1(2) * v2(0) - v1(0) * v2(2),
	v1(0) * v2(1) - v1(1) * v2(0) );
  }
  //����ʽ��ֵ��ע�����������������
  inline double Determinant (const Vec<3> & col1, const Vec<3> & col2, const Vec<3> & col3)
  {
    return
      col1(0) * ( col2(1) * col3(2) - col2(2) * col3(1)) +
      col1(1) * ( col2(2) * col3(0) - col2(0) * col3(2)) +
      col1(2) * ( col2(0) * col3(1) - col2(1) * col3(0));
  }

  //������׾�����棨�������أ�
  // template <int H, int W>
  inline void CalcInverse (const Mat<2,2> & m, Mat<2,2> & inv)
  {
	  double det = m(0,0) * m(1,1) - m(0,1) * m(1,0);
	  if (det == 0) 
	  {
		  inv = 0;
		  return;
	  }

	  double idet = 1.0 / det;
	  inv(0,0) =  idet * m(1,1);
	  inv(0,1) = -idet * m(0,1);
	  inv(1,0) = -idet * m(1,0);
	  inv(1,1) =  idet * m(0,0);
  }
  //�������׾�����棨�������أ�
  void CalcInverse (const Mat<3,3> & m, Mat<3,3> & inv);

  inline void CalcInverse (const Mat<2,3> & m, Mat<3,2> & inv)
  {
    Mat<2,2> a = m * Trans (m);
    Mat<2,2> ainv;
    CalcInverse (a, ainv);
    inv = Trans (m) * ainv;
  }
  //����3*2�׾�����棨�������أ�
  void CalcInverse (const Mat<3,2> & m, Mat<2,3> & inv);

  inline void CalcInverse (const Mat<3,2> & m, Mat<2,3> & inv)
  {
    Mat<2,2> a = Trans (m) * m;
    Mat<2,2> ainv;
    CalcInverse (a, ainv);
    inv = ainv * Trans (m);
  }

  double Det (const Mat<2,2> & m);
  double Det (const Mat<3,3> & m);

  // eigenvalues of a symmetric matrix
  void EigenValues (const Mat<3,3> & m, Vec<3> & ev);
  void EigenValues (const Mat<2,2> & m, Vec<3> & ev);

}

#endif
