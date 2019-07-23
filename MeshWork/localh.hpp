#ifndef LOCALH
#define LOCALH

/**************************************************************************/
/* File:   localh.hh                                                      */
/* Function: Control of 3D mesh grading and restrict sphere radius        */
/* Date:   24. jan. 2016                                                   */
/**************************************************************************/


namespace meshwork
{

  /// box for grading
  class GradingBox
  {
    /// xmid
    double xmid[3];
    /// half edgelength
    double h2;
    ///
    GradingBox * childs[8];
    ///
    GradingBox * father;
    ///
    double hopt;
    ///
  public:

    //struct 
    //{
    //  unsigned int cutboundary:1;
    //  unsigned int isinner:1;
    //  unsigned int oldcell:1;
    //  unsigned int pinner:1;
    //} flags;

    ///
    GradingBox (const double * ax1, const double * ax2);
    ///
    void DeleteChilds();
    ///

    Point<3> PMid() const 
	{ 
		return Point<3> (xmid[0], xmid[1], xmid[2]); 
	}
    double H2() const 
	{ 
		return h2; 
	}

    friend class LocalH;

  };




  /**
     Control of 3D mesh grading
  */
  class LocalH 
  {
    ///
    GradingBox * root;
    ///
    double grading;
    ///
    vector<GradingBox*> boxes;
    ///
    Box<3> boundingbox;
  public:
    ///
    LocalH (const Point<3> & pmin, const Point<3> & pmax, double grading);
    ///
    LocalH (const Box<3> & box, double grading);
    ///
    ~LocalH();
    ///
    void Delete();
    ///
    void SetGrading (double agrading)
	{ 
		grading = agrading;
	}
    ///
    void SetH (const Point<3> & x, double h);
    ///
    double GetH (const Point<3> & x) const;
    /// minimal h in box (pmin, pmax)
    double GetMinH (const Point<3> & pmin, const Point<3> & pmax) const;

    /// mark boxes intersecting with boundary-box
    // void CutBoundary (const Point3d & pmin, const Point3d & pmax)
    // { CutBoundaryRec (pmin, pmax, root); }
    void CutBoundary (const Box<3> & box)
    { 
		CutBoundaryRec (box.PMin(), box.PMax(), root); 
	}
  

    /// widen refinement zone
    void WidenRefinement ();

    int GetNBoxes () 
	{
		return boxes.size();
	} 
    const Box<3> & GetBoundingBox () const
    { 
		return boundingbox; 
	}
    ///
    void PrintMemInfo (ostream & ost) const;
  private:
    /// 
    double GetMinHRec (const Point<3> & pmin, const Point<3> & pmax,const GradingBox * box) const;
    ///
    void CutBoundaryRec (const Point<3> & pmin, const Point<3> & pmax,GradingBox * box);

    friend ostream & operator<< (ostream & ost, const LocalH & loch);
  };




  inline ostream & operator<< (ostream & ost, const GradingBox & box)
  {
    ost << "gradbox, pmid = " << box.PMid() << ", h2 = " << box.H2() 
	<< endl;
    return ost;
  }

  inline ostream & operator<< (ostream & ost, const LocalH & loch)
  {
    for (int i = 0; i < loch.boxes.size(); i++)
      ost << "box[" << i << "] = " << *(loch.boxes[i]);
    return ost;
  }

}

#endif
