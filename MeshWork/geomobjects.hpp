#ifndef FILE_OBJECTS
#define FILE_OBJECTS

/* *************************************************************************/
/* copy from netgen  and modify it                                                 */
/* File:   geomobjects.hpp                                                 */
/* Function:  describe fundamental type of geomotry                        */ 
/* Date:   08. jan. 2016                                                   */
/* Class: Point Vec Mat Box BoxSphere                                      */
/* *************************************************************************/

namespace meshwork
{

  template <int D> class Vec;
  template <int D> class Point;


  template <int D>
  class Point
  {

  protected:
    double x[D];

  public:
    Point () { ; }
    Point (double ax) 
	{
		for (int i = 0; i < D; i++) x[i] = ax; 
	}
    Point (double ax, double ay) 
    { 
      // static_assert(D==2, "Point<D> constructor with 2 args called");
      x[0] = ax; x[1] = ay; 
    }
    Point (double ax, double ay, double az) 
    {
      // static_assert(D==3, "Point<D> constructor with 3 args called");
      x[0] = ax; x[1] = ay; x[2] = az; 
    }
    Point (double ax, double ay, double az, double au)
    {
		 // static_assert(D==4, "Point<D> constructor with 4 args called");
		x[0] = ax; x[1] = ay; x[2] = az; x[3] = au;
	}

    Point (const Point<D> & p2)
    {
		for (int i = 0; i < D; i++) x[i] = p2.x[i]; 
	}

    explicit Point (const Vec<D> & v)
    {
		for (int i = 0; i < D; i++) x[i] = v(i); 
	}


    Point & operator= (const Point<D> & p2)
    {
      for (int i = 0; i < D; i++) x[i] = p2.x[i]; 
      return *this;
    }

    Point & operator= (double val)
    {
      for (int i = 0; i < D; i++) x[i] = val;
      return *this;
    }
	//两点取小
	const Point & SetToMin (const Point<D> & p2)
	{
		for(int i = 0; i < D; i++)
		if (p2.x[i] < x[i]) x[i] = p2.x[i];
		return *this;
	}
	//两点取大
	const Point & SetToMax (const Point<D> & p2)
	{
		for(int i = 0; i < D; i++)
		if (p2.x[i] > x[i]) x[i] = p2.x[i];
		return *this;
	}

    double & operator() (int i) 
	{ 
		return x[i]; 
	}
    const double & operator() (int i) const
	{
		return x[i]; 
	}
	double & operator[] (int i) 
	{ 
		return x[i]; 
	}
    const double & operator[] (int i) const
	{ 
		return x[i];
	}

    operator const double * () const 
	{ 
		return x;
	}
	//更改第i个值
	void Modify ( int  i , double val )
	{
		x[i] = val; 
	}
	//更改所有的值
	void Modify ( double val )
	 { 
		 for (int i = 0; i < D; i++) x[i] = val;
	}

  };


  template <int D>
  class Vec
  {

  protected:
    double x[D];

  public:
    Vec () { ; } // for (int i = 0; i < D; i++) x[i] = 0; }
    Vec (double ax) 
	{ 
		for (int i = 0; i < D; i++) x[i] = ax;
	}
    Vec (double ax, double ay) 
    { 
      // static_assert(D==2, "Vec<D> constructor with 2 args called");
      x[0] = ax; x[1] = ay; 
    }
    Vec (double ax, double ay, double az)
    { 
      // static_assert(D==3, "Vec<D> constructor with 3 args called");
      x[0] = ax; x[1] = ay; x[2] = az; 
    }
    Vec (double ax, double ay, double az, double au)
    { 
		x[0] = ax; x[1] = ay; x[2] = az; x[3] = au; 
	}

    Vec (const Vec<D> & p2)
    { 
		for (int i = 0; i < D; i++) x[i] = p2.x[i]; 
	}

    explicit Vec (const Point<D> & p)
    { 
		for (int i = 0; i < D; i++) x[i] = p(i); 
	}

    Vec (const Vec<D> & p1, const Vec<D> & p2)
    { 
		for(int i=0; i<D; i++) x[i] = p2(i)-p1(1);
	}
  
    Vec & operator= (const Vec<D> & p2)
    {
      for (int i = 0; i < D; i++) x[i] = p2.x[i]; 
      return *this;
    }

    Vec & operator= (double s)
    {
      for (int i = 0; i < D; i++) x[i] = s;
      return *this;
    }

    double & operator() (int i) 
	{
		return x[i]; 
	}
    const double & operator() (int i) const 
	{ 
		return x[i]; 
	}

    operator const double* () const 
	{ 
		return x; 
	}

	void Modify ( int  i , double val )
	{ 
		x[i] = val; 
	}
	void Modify ( double val )
	{ 
		for (int i = 0; i < D; i++) 
			x[i] = val; 
	}
	//向量的长度
    double Length () const
    {
      double l = 0;
      for (int i = 0; i < D; i++)
	      l += x[i] * x[i];
      return sqrt (l);
    }
	//向量长度的平方
    double Length2 () const
	{
		double l = 0;
		for (int i = 0; i < D; i++)
			l += x[i] * x[i];
		return l;
	}
	//向量归一化
    void Normalize ()
	{
		double l = Length();
		if (l != 0)
			for (int i = 0; i < D; i++)
				x[i] /= l;
	}
	//获取向量方向
    Vec<D>& GetNormal () const
	{
		Vec<D> vec; 
		double l = Length();
		if (l != 0)
			for (int i = 0; i < D; i++)
				vec.x[i] = x[i] / l;
		return vec;
	}
	
  };




  template <int H, int W=H>
  class Mat
  {
  protected:
    double x[H*W];

  public:
    Mat () { ; }
    Mat (const Mat & b)
    { 
		for (int i = 0; i < H*W; i++) x[i] = b.x[i]; 
	}
  
    Mat & operator= (double s)
    {
      for (int i = 0; i < H*W; i++) x[i] = s;
      return *this;
    }

    Mat & operator= (const Mat & b)
    {
      for (int i = 0; i < H*W; i++) x[i] = b.x[i]; 
      return *this;
    }

    double & operator() (int i, int j) 
	{ 
		return x[i*W+j]; 
	}
    const double & operator() (int i, int j) const 
	{ 
		return x[i*W+j]; 
	}
    double & operator() (int i) 
	{ 
		return x[i]; 
	}
    const double & operator() (int i) const 
	{ 
		return x[i]; 
	}
	//获取第i列
    Vec<H> Col (int i) const
	{
		Vec<H> hv; 
		for (int j = 0; j < H; j++)
			hv(j) = x[j*W+i];
		return hv; 
	}
	//获取第i行
    Vec<W> Row (int i) const
	{
		Vec<W> hv; 
		for (int j = 0; j < W; j++)
			hv(j) = x[i*W+j];
		return hv; 
	}

/*
	解矩阵方程还没有实现
    void Solve (const Vec<H> & rhs, Vec<W> & sol) const
    {
      Mat<W,H> inv;
      CalcInverse (*this, inv);
      sol = inv * rhs;
    }
*/
	void Modify ( int  i , int  j, double val ) 
	{ 
		x[i*W+j] = val; 
	}
	void Modify ( int  i , double val )
	{ 
		x[i] = val; 
	}
	void Modify ( double val )
	{ 
		for (int i = 0; i < H*W; i++) x[i] = val; 
	}
  };

  template <int D>
  class Box
  {
  protected:
    Point<D> pmin, pmax;
  public:
    Box () { ; }

    Box ( const Point<D> & p1)
	{
		for (int i = 0; i < D; i++)
			pmin(i) = pmax(i) = p1(i);
	}

    Box ( const Point<D> & p1, const Point<D> & p2)
	{
		for (int i = 0; i < D; i++)
		{
			pmin(i) = min(p1(i), p2(i));
			pmax(i) = max(p1(i), p2(i));
		}
	}

    const Point<D> & PMin () const 
	{ 
		return pmin; 
	}
    const Point<D> & PMax () const 
	{ 
		return pmax; 
	}
  
    void Set (const Point<D> & p)
    { 
		pmin = pmax = p; 
	}

    void Add (const Point<D> & p)
	{ 
		for (int i = 0; i < D; i++)
		{
			if (p(i) < pmin(i)) pmin(i) = p(i);
			else if (p(i) > pmax(i)) pmax(i) = p(i);
		}
	}

	double & operator[] (int i) 
	{
		if(i<=D) 
			return pmin[i];
		else 
			return pmax[i-D];
	}
	const	double & operator[] (int i) const
	{
		if(i<=D) 
			return pmin[i];
		else 
			return pmax[i-D];
	}
/*
    template <typename T1, typename T2>
    void Set (const IndirectArray<T1, T2> & points)
    {
      Set (points[points.Begin()]);
      for (int i = points.Begin()+1; i < points.End(); i++)
        Add (points[i]);
    }

    template <typename T1, typename T2>
    void Add (const IndirectArray<T1, T2> & points)
    {
      for (int i = points.Begin(); i < points.End(); i++)
        Add (points[i]);
    }
*/
	//盒子的中心
    Point<D> Center () const 
	{ 
		Point<D> c;
		for (int i = 0; i < D; i++)
			c(i) = 0.5 * (pmin(i)+pmax(i)); 
		return c;
	}
	//对角线的长度
    double Diam () const 
	{ 
		return Abs (pmax-pmin); 
	}
	//获取第nr个点
    Point<D> GetPointNr (int nr) const
	{
		Point<D> p;
		for (int i = 0; i < D; i++)
		{
			p(i) = (nr & 1) ? pmax(i) : pmin(i);
			nr >>= 1;
		}
		return p;
	}

    bool Intersect (const Box<D> & box2) const
	{
		for (int i = 0; i < D; i++)
			if (pmin(i) > box2.pmax(i) ||
				pmax(i) < box2.pmin(i)) return 0;
		return 1;
	}
	//判断点是否在盒子内
    bool IsIn (const Point<D> & p) const
	{
		for (int i = 0; i < D; i++)
			if (p(i) < pmin(i) || p(i) > pmax(i)) return 0;
		return 1;
	}
	//增加盒子
    void Increase (double dist)
	{
		for (int i = 0; i < D; i++)
		{
			pmin(i) -= dist;
			pmax(i) += dist;
		}
	}
  };


  template <int D>
  class BoxSphere : public Box<D>
  {
  protected:
	  ///
	  Point<D> c;
	  ///
	  double diam;
	  ///
	  double inner;
  public:
	  ///
	  BoxSphere () { ; }
	  ///
	  BoxSphere (const Box<D> & box) : Box<D> (box) 
	  { 
		  CalcDiamCenter();
	  }

	  ///
	  BoxSphere ( Point<D> apmin, Point<D> apmax ): Box<D> (apmin, apmax)
	  {
		  CalcDiamCenter();
	  }

	  ///
	  const Point<D> & Center () const 
	  { 
		  return c; 
	  }
	  ///
	  double Diam () const 
	  { 
		  return diam; 
	  }
	  ///
	  double Inner () const 
	  { 
		  return inner; 
	  }


	  ///
	  void GetSubBox (int nr, BoxSphere & sbox) const
	  {
		  for (int i = 0; i < D; i++)
		  {
			  if (nr & 1)
			  {
				  sbox.pmin(i) = c(i);
				  sbox.pmax(i) = this->pmax(i);
			  }
			  else
			  {
				  sbox.pmin(i) = this->pmin(i);
				  sbox.pmax(i) = c(i);
			  }
			  sbox.c(i) = 0.5 * (sbox.pmin(i) + sbox.pmax(i));
			  nr >>= 1;
		  }
		  sbox.diam = 0.5 * diam;
		  sbox.inner = 0.5 * inner;
	  }


	  ///
	  void CalcDiamCenter ()
	  {
		  c = Box<D>::Center ();
		  diam = Dist (this->pmin, this->pmax);

		  inner = this->pmax(0) - this->pmin(0);
		  for (int i = 1; i < D; i++)
			  if (this->pmax(i) - this->pmin(i) < inner)
				  inner = this->pmax(i) - this->pmin(i);
	  }

  };


}


#endif
