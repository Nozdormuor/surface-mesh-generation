#ifndef FILE_TABLE
#define FILE_TABLE

/**************************************************************************/
/* File:   table.hpp                                                      */
/* Date:   15. April. 2016                                                    */
/**************************************************************************/

namespace meshwork 
{
	///类的可见性设计还存在问题
	/// generic class TABLE.
	template <class T>
	class Table
	{
	public:

		///
		class linestruct
		{
		public:
			///
			int size;
			/// 
			int maxsize;
			///
			T * col;
			T operator[] (int i) const
			{ 
				return col[i];
			}
		};

		///
		vector<linestruct> data;

	public:
		///
		Table (int size)
		{
			data.resize(size);
			for (int i = 0; i < size; i++)
			{
				data[i].maxsize = 0;
				data[i].size = 0;
				data[i].col = NULL;
			}
		}
		///
		Table ()
		{
			;
		}
		///
		~Table ()
		{
			for (int i = 0; i < data.size(); i++)
			{
				if(data[i].col != NULL)
					delete [] (T*)data[i].col;
			}	
		}
		///
		void SetSize (int size)
		{
			for (int i = 0; i < data.size(); i++)
			{
				if(data[i].col != NULL)
					delete [] (T*)data[i].col;
			}	

			data.resize(size);
			for (int i = 0; i < size; i++)
			{
				data[i].maxsize = 0;
				data[i].size = 0;
				data[i].col = NULL;
			} 
		}
		///
		void ChangeSize (int size)
		{
			int oldsize = data.size();
			if (size == oldsize) 
				return;

			if (size < oldsize)
				for (int i = size; i < oldsize; i++)
				{
					if(data[i].col != NULL)
						delete [] (T*)data[i].col;
				}	

				data.SetSize(size);

				for (int i = oldsize; i < size; i++)
				{
					data[i].maxsize = 0;
					data[i].size = 0;
					data[i].col = NULL;
				}    
		}
		///
		void Add (int i, const T & acont)
		{
			int oldsize = data.size();
			if(data.size()<=i)
			{
				data.resize(i+1);
				for (int j = oldsize; j <= i; j++)
				{
					data[j].maxsize = 0;
					data[j].size = 0;
					data[j].col = NULL;
				} 
			}
			linestruct & line = data[i];
			if (line.size == line.maxsize)
			{
				T * p = new T [(line.maxsize+5)];

				memcpy (p, line.col, line.maxsize*sizeof(T));
				delete [] (T*)line.col;

				line.col = p;
				line.maxsize += 5;
			}

			line.size++;

			((T*)data[i].col)[data[i].size-1] = acont;
		}
		///共几行
		int Size () const
		{
			return data.size();
		}
		///row行共几列
		int RowSize (int row) const
		{
			return data[row].size;
		}
		///没有检查是否越界，请慎用
		vector<T> operator[] (int i) const
		{ 
			return vector<T> (data[i].col,data[i].col+data[i].size*sizeof(T));
		}
		///没有检查是否越界，请慎用
		void Set (int i, int nr, const T & acont)
		{ 
			data[i].col[j] = acont; 
		}
		///没有检查是否越界，请慎用
		T Get(int i,int j) const
		{
			return data[i][j];
		}
	};

	template <class T>
	ostream & operator<< (ostream & ost, const Table<T> & table)
	{
		for (int i = 0; i < table.Size(); i++)
		{
			ost << i << ": ";
			ost << "(" << table.data[i].size << ") ";
			for (int j = 0; j < table.data[i].size; j++)
				ost << table.data[i].col[j] << " ";
			ost << endl;
		}
		return ost;
	}




}
#endif