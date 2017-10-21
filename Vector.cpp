#include <iostream>
#include <fstream>
#include <string>
#include "math.h"
#include "general.h"
#include "Vector.h"
using namespace std;

/*
	int count_size(const string & filename) 
	{
	ifstream fin(filename.c_str());
	int fsize = std::distance(std::istream_iterator<double>(fin),
		std::istream_iterator<double>());
	fin.close();
	return fsize;
	}

Implement the following functions:
	bool Read(const string& filename, Vector& A)	//read the data from a file and create a Vector A;
													//if A has data already, delete them first.
													//return false if error occurs
	{
	cout<<"Reading data from: "<<filename<<endl;	
	int fsize = count_size(filename);
	//Write the vector to file
	fstream myfile(filename, ios_base::in);    		// read file
	int val=0;
		for (int i = 0; i < fsize; i++)
		{
			myfile >> val;
			A.buf[i] = val;
		}
		myfile.close();	
		A.size = fsize; 					 	//do NOT forget to close the file
		return true;
	}
	
	
	bool Write(const string& filename, const Vector& A)	//write the Vector data to a file;
														//return false if error occurs
														//refer to Read()
	{
	//Write the vector to file
	fstream fout(filename, ios_base::out);
	for (int i = 0; i < A.size; i++)
	{
		fout << A.buf[i] <<"\t";
	}
	fout.close();
	return true;
	}
	*/


	//Implement the following functions:

	bool Read(const string& filename, Vector& A)
		//read the data from a file and create a Vector A;
		//if A has data already, delete them first.
		//return false if error occurs
	{
		string line;
		cout << "Reading data from: " << filename << endl;
		ifstream ins(filename.c_str());
		BypassComment(ins);//skip comments

						   //count data in input file
		int ct = 0;
		int x;
		while (ins >> x) {
			ct++;
			//cout<<ct<<" data counted: "<<x<<endl;
		}
		cout << "Total data in file: " << ct << endl;
		ins.close();

		//load values into vector		
		A.size = ct;
		A.buf = new double[ct];

		//debug: print previous contents of buffer		
		for (int indata = 0; indata<A.size; indata++) {
			cout << "Prevous value #" << indata << ": " << A.buf[indata] << endl;
		}

		int i = 0;
		int y;
		ifstream ins2(filename.c_str());	//reopen file
		BypassComment(ins2);			//skip comments
										//Load contents into new value
		while (ins2 >> y) {
			A.buf[i] = y;
			cout << "New value #" << i << ": " << A.buf[i] << endl;
			i++;
		}

		//close the file
		ins2.close();


		//print loaded matrix
		cout << "New vector: " << endl;
		for (int j = 0; j<A.size; j++) {
			cout << ".buf[" << j << "] = " << A.buf[j] << endl;
		}
		cout << "[end]" << endl;


		return true;
	}


	bool Assign (Vector& A,const Vector& B) //assign B to A, if A has data, delete them first
											//return false if error occurs
	{
		Clear(A);
		A.size = B.size;
		A.buf = new double[A.size];
		for (int i = 0; i < B.size; i++)
		{
			A.buf[i] = B.buf[i];
		}                        
		return true;
	}


	bool Equal(const Vector& A,const  Vector& B) //compare the two, return ture if they are equal to each other, otherwise return false;
	{
		//Please note that if two floating point sizebers are compared, you may get wrong result!
		if (A.size != B.size)
			return false;
		else
		for (int i = 0; i < B.size; i++)
		{
			if (A.buf[i] != B.buf[i])
				return false;
			else
				return true;
		}
	}									

	bool Add(const Vector&  A, const Vector& B, Vector& C)// C=A+B
	{
		Clear(C);
		C.size = A.size;
		if(B.size < C.size)
			C.size = B.size;
		C.buf = new double[C.size];
		for(int i = 0; i < C.size; i++)
		{
			C.buf[i] = A.buf[i] + B.buf[i];
		}
		return true;
	}


	bool Subtract(const Vector& A,const Vector& B, Vector& C)//C=A-B
	{									                
		Clear(C);
		C.size = A.size;
		if(B.size < C.size)
			C.size = B.size;
		C.buf = new double[C.size];
		for(int i = 0; i < C.size; i++)
		{
			C.buf[i] = A.buf[i] - B.buf[i];
		}
		return true;
	}


	bool Multiply(double a,const Vector& B,Vector& C) //C=a*B;
	{
		Clear(C);
		C.size = B.size;
		C.buf = new double[C.size];
		for(int i = 0; i < C.size; i++)
		{
			C.buf[i] = a * B.buf[i];
		}
		return true;
	}


	bool Multiply(const Vector& A, const Vector& B, Vector C)  //dot product
	{
		Clear(C);
		C.size = A.size;
		if(B.size < C.size)
			C.size = B.size;
		C.buf = new double[C.size];
		for(int i = 0; i < C.size; i++)
		{
			C.buf[i] = A.buf[i] * B.buf[i];
		}
		return true;
	}


	void Initialize(Vector& A)
	{
		A.size=0;			//
		A.buf=NULL;			//pointing to nothing
	}


	void Clear(Vector& A)  //remove all the data from the Vector.
	{
		A.size=0;//
		if(A.buf) delete[] A.buf;		
		A.buf=NULL;//pointing to nothing
	}


	void Reverse(Vector& A)
	{
		int size = A.size;
		double* buf = new double[size];
		for (int i = size - 1; i >= 0; i--)
			buf[size - 1 - i] = A.buf[i];
		for (int i = 0; i < size; i--)
			A.buf[i] = buf[i];
		delete[] buf;
	}


	bool Remove(Vector& A, double b)          // b 这个数是吧？
	{
		int bCount = 0;
		for(int i = 0; i < A.size; i++)
		{
			if(A.buf[i] == b)
				bCount++;
		}
		if(bCount <= 0)
			return true;
		double* buf = new double[A.size-bCount];
		int j = 0;
		for (int i = 0; i < A.size; i++)
		{
			if(A.buf[i] != b)
				buf[j++] = A.buf[i];
		}
		delete[] A.buf;
		A.buf = buf;
	}

	bool Remove(Vector& A, int location)    //remove location of vector.
	{
		if(location >= A.size)
			return false;
		double* buf = new double[A.size-1];
		int j = 0;
		for (int i = 0; i < A.size; i++)
		{
			if(i != location)
				buf[j++] = A.buf[i];
		}
		delete[] A.buf;
		A.buf = buf;             
	}	
	
	bool Insert(Vector& A, double b, int location)
	{
		if(location >= A.size)
			return false;
		double* buf = new double[A.size+1];
		int j = 0;
		for (int i = 0; i < A.size; i++)
		{
			if(i == location)
				buf[j++] = b;
			buf[j++] = A.buf[i];
		}
		delete[] A.buf;
		A.buf = buf;
		return true;
	}


	void Sort(Vector& A, bool Ascending)
	{
		double temp = 0;
		for(int i = 0; i < A.size; i++)
		{
			for(int j = 0; j < A.size - i; j++)
			{
				if(A.buf[j] > A.buf[j+1])
				{
					temp = A.buf[j];
					A.buf[j] = A.buf[j+1];
					A.buf[j+1] = temp;
				}
			}
		}
	}


	void Unit(Vector& A)    //Make Vector A unit vector of A. A->|A|=1
	{
		double sum(0.0);
		double mod;
		for (int i = 0; i < A.size; i++)
		{
			sum += A.buf[i] * A.buf[i];
		}
		mod = sqrt(sum);
		for (int i = 0; i < A.size; i++)
		{
			A.buf[i] /= mod;
		};
	}


	bool IsUnit(const Vector& A)		//check if Vector is Unit Vector |A|=1?
	{
		double sum(0.0);
		double mod;
		for (int i = 0; i < A.size; i++)
		{
			sum += A.buf[i] * A.buf[i];
		}
		mod = sqrt(sum);
		if (mod == 1)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
