//
// Created by ccfy on 18-7-28.
//
#include <iostream>
#include <iterator>
#include <algorithm>
using namespace std;

/***********************************************************************
int abc(int a, int b,int c)
{
	if(a <=0 || b <=0 || c <=0)
	{
		throw "All para shoule be > 0";
	}
	return a + b * c;
}

int main(int argc, char** argv)
{
	try {
		cout<< abc(2,0,1);
	}
	catch (char* e)
	{
		cout<<"The Parameters are 2, 0 ,1 for abc function"<<endl;
		cout<<"An exception has been throw"<<endl;
		cout<<e<<endl;
		return 1;
	}
	return 0;
}

template <class T>
void make2dArray(T** &x, int row, int col)
{
	x = new T* [row];
	for(int i=0;i<row;i++)
	{
		x[i] = new T[col];
	}
}

template <class T>
void delete2dArray(T** &x , int row)
{
	for(int i=0;i<row;i++)
	{
		delete [] x[i];
	}
	delete []x;
}


template <class T>
void permucation(T list[], int k, int m)
{
	if(k == m)
	{
		copy(list,list+m+1,ostream_iterator<T>(cout,""));
		cout<<endl;
	}
	else
	{
		for(int i=k;i<=m;i++)
		{
			swap(list[k], list[i]);
			permucation(list,k+1,m);
			swap(list[k],list[i]);
		}
	}
}

template <class T>
void permutation(T list[], int k, int m)
{
	do{
		copy(list,list+m+1,ostream_iterator<T>(cout,""));
		cout<<endl;
	}while(next_permutation(list+k,list+m+1));
}

template <class T>
void ranking(T a[],int n, T r[])
{
	for(int i=0;i<n;i++)
	{
		r[i] = 0;
	}
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			if(j == i)
				continue;
			if(a[j]<a[i])
				r[i] = r[i] +1;
			if(a[j] == a[i] && j < i)
				r[i] = r[i] +1;
		}
	}
}

template <class T>
void countSort(T a[],int n,T r[])
{
	T * temp = new T[n];
	for(int i=0;i<n;i++)
	{
		temp[r[i]] = a[i];
	}
	for(int i=0;i<n;i++)
	{
		a[i] = temp[i];
	}
	delete []temp;
}

template <class T>
void countSort_inplace(T a[],int n,T r[])
{
	for(int i=0;i<n;i++)
	{
		if(r[i] == i)
			continue;
		else
		{
			swap(a[i],a[r[i]]);
			swap(r[i],r[r[i]]);
			countSort_inplace(a,n,r);
		}
	}
}

template <class T>
void countSort_inplace(T a[],int n,T r[])
{
	for(int i=0;i<n;i++)
	{
		while (r[i]!=i)
		{
			swap(a[i],a[r[i]]);
			swap(r[i],r[r[i]]);
		}
	}
}

template <class T>
void sparseMatrix<T>::transpose(sparseMatrix<T>& b)
{
	b.cols = rows;
	b.rows = cols;
	b.terms.reSet(terms.size());
	int *colSize = new int[cols+1];
	int *rowNext = new int[rows+1];
	for(int i =1;i<=cols;i++)
	{
		colSize[i] = 0;
	}
	// get element number of origin matrix col
	for(auto i = terms.begin();i != terms.end();i++)
	{
		colSize[(*i).col]++;
	}
	rowNext[1] = 0;
	for(int i = 2;i<=cols;i++)
	{
		rowNext[i] = rowNext[i - 1] + colSize[i-1];
	}
	// rowNext[i] is the index of first element of ith row of transpose matrix
	MatrixTerm<T> mTerm;
	for(arrayList<MatrixTerm<T>>::iterator it= terms.begin();it!=terms.end();it++)
	{
		mTerm.row = (*it).col;
		mTerm.col = (*it).row;
		mTerm.value = (*it).value;
		int j = rowNext[(*it).col];
		//先遇到的同一列的元素，在原有的线性表中行靠前
		rowNext[(*it).col]++;
		b.terms.set(j,mTerm);
	}
}

template <class T>
void sparseMatrix<T>::add(sparseMatrix<T>& b, sparseMatrix<T>& c)
{
	if(rows != b.rows || cols !=b.cols)
		throw MatrixSizeMisMatch();
	c.rows = rows;
	c.cols = cols;
	c.terms.clear();
	int csize =0;
	arrayList<MatrixTerm<T>>::iterator it = terms.begin();
	arrayList<MatrixTerm<T>>::iterator ib = b.terms.begin();
	arrayList<MatrixTerm<T>>::iterator itend = terms.end();
	arrayList<MatrixTerm<T>>::iterator ibend = terms.end();
	while(it!=itend && ib!=ibend)
	{
		int tindex = (*it).row * cols + (*it).col;
		int bindex = (*ib).row * cols + (*ib).col;
		if(tindex < bindex)
		{
			c.terms.insert(csize++,(*it));
			it ++;
		}
		else if(tindex > bindex)
		{
			c.terms.insert(csize++,(*ib));
			ib ++;
		}
		else
		{
			if(((*it).value + (*ib).value)!=0)
			{
				MatrixTerm<T> mTerm;
				m.row = (*it).row;
				m.col = (*it).col;
				m.value = (*it).value + (*ib).value;
				c.terms.insert(csize++,mTerm);
				it++;
				ib++;
			}
		}
	}
	for(;it!=itend;it++)
	{
		c.insert(csize++, (*it));
	}
	for(;ib!=ib.end();ib++)
	{
		c.insert(csize++,*ib);
	}
}



void printMatchedPairs(string expr)
{
	arrayStack<int>s;
	int length = expr.size();
	for(int i=0;i<length;i++)
	{
		if(expr.at(i) == '(')
			s.push(i);
		else
		{
			if(expr.at(i) == ')')
			{
				try
				{
					cout<<s.top()<<' '<<i<<endl; //left right match
					s.pop();
				}
				catch (stackEmpty)
				{
					cout<<"no match for right parenthesis "<<"at "<<i<<endl;
				}
			}
		}
		while(!s.empty())
		{
			cout<<"no match for left parenthesis at "<<i<<endl;
			s.pop();
		}
	}
}

template <class T>
void changeLength1D(T* &a,int oldLength,int newLength)
{
	T* temp = new T[newLength];
	int num = min(oldLength,newLength);
	copy(a,a+num,temp);
	delete []a;
	a = temp;
}


void towerOfHanoi(int n,int x,int y,int z)
{
	//move n hanoi from x to y vai z 数学归纳,递归推导
	if(n>0)
	{
		towerOfHanoi(n-1,x,z,y);
		cout<<"move top disk from "<<x<<" to "<<y<<endl;
		towerOfHanoi(n-1,z,y,x);
	}
}


//use stack as the hanoi's data structure
arrayStack<int>tower[4];
void moveAndShow(int n,int x,int y, int z)
{
	if(n>0)
	{
		moveAndShow(n-1,x,z,y);
		int d = tower[x].top();
		tower[x].pop();
		tower[y].push(d);
		showstate();
		moveAndShow(n-1,z,y,x);
}

void towerOfHanoi(int n)
{
   //push disk into stack one from big to small
	for(int i=n;i>0;i--)
	{
		tower[1].push(i);
	}
	moveAndShow(n,1,2,3);
}

*****************************************************************/

//列车重排列 入轨 缓冲轨道 出轨
arrayStack<int> *track;
int numberOfCars;
int numberOfTracks;
// 缓冲轨道中编号最小的车厢以及轨道
int smallestCar;
int itsTrack;

void outputFromHoldingTrack()
{
    track[itsTrack].pop();
	cout<<"move car "<<smallestCar<<" from holding track to output track"<<endl;
	smallestCar = numberOfCars + 1;
	for(int i=1;i<=numberOfTracks;i++)
	{
		if(!track[i].empty() && track[i].top()<smallestCar)
		{
			smallestCar = track[i].top();
			itsTrack = i;
		}
	}
}

bool putintoHoldingTrack(int c)
{
	//find best track to insert c
	int bestTrack = 0;
	int bestTop = numberOfCars + 1; // current min number of track id of holding track
	for(int i = 1;i<=numberOfTracks;i++)
	{
		if(!track[i].empty())
		{
			int topcar = track[i].top();
			if(c<topcar && topcar < bestTop)
			{
				bestTrack =i;
				bestTop = topcar;
			}
		}
		else
		{
			if(bestTrack == 0)
				bestTrack = i;
		}
	}
	if(bestTrack == 0)
		return false;
	track[bestTrack].push(c);
	cout<<"move car"<<c<<" from input track to holding track"<<endl;
	if(c<smallestCar)
	{
		smallestCar = c;
		itsTrack = bestTrack;
	}
	return true;
}

bool railroad(int inputOrder[],int thenumberOfTrack,int thenumberOfCar)
{
	numberOfCars = thenumberOfCar;
	numberOfTracks = thenumberOfTrack;
	track  = new arrayStack<int>[numberOfTracks];
	int nextCarToOutput = 1;
	smallestCar = numberOfCars + 1;
	for(int i=0;i<numberOfCars;i++)
	{
		if(inputOrder[i] == nextCarToOutput)
		{
			cout<<"move car "<<inputOrder[i]<<" from input track to output track"<<endl;
			nextCarToOutput++;
			while(smallestCar == nextCarToOutput)
			{
				outputFromHoldingTrack();
				nextCarToOutput++;
			}
		}
		else
		{
			if(!putintoHoldingTrack(inputOrder[i]))
				return false;
		}
	}
	return true;
}

int main()
{

}