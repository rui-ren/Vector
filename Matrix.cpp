coding:utf-8
#include "Matrix.h"
#include "general.h"
#include <math.h>
//------------------------------ Implemenation  -----------------------------

Matrix::Matrix() //default constructor
{
	rows=columns=0;
	buf=NULL;
}

Matrix::Matrix( int m,  int n)//declare an mxn matrix
{
	rows = m;
	columns = n;
	buf = new double[m*n];
}

Matrix::Matrix(const Matrix& A) //copy constructor
{
	if (!this->buf)	delete[]this->buf;   // 判断buf 中是否有数据？
	this->rows = A.rows;
	this->columns = A.columns;          // copy rows and columns using this pointer
	this->buf = new double[A.rows*A.columns];  // allocate the hard dirve
	for (int i = 0; i < A.rows*A.columns; i++)
		this->buf[i] = A.buf[i];
}

Matrix::~Matrix()//destructor
{
	this->rows = 0;
	this->columns = 0;
	delete[] this->buf;    // deconstructor clean up.
}
//Assignment
Matrix& Matrix::operator = (const Matrix& A) //overloading =
{	
	if (this == &A)	return *this;
	if (this->buf) delete[] buf;
	columns = A.columns;
	rows = A.rows;
	this->buf = new double[columns*rows];
	if (int i=0; i < (columns * A.rows); i++)
		this->buf[i] = A.buf[i];
	return *this;
}
//operators
bool Matrix::operator == (const Matrix& A)//overloading ==
{
	if (!this->buf || !A.buf)    // two empty matrix && 为什么不是and
		return true;
	if (this->columns != A.columns || this->rows != A.rows)
		return false;            
	else
	{
		for (int i = 0; i < columns*rows; i++)
			this->buf[i] = A.buf[i];         // 上次作业这个就错了， 我觉得应该用相减   如果98.0=98；
			/*
			if (abs(this->buf[i]-A.buf[i]>1e-10))
			return false;
			*/
	}
	return true;
}


bool Matrix::operator != (const Matrix& A)//overloading !=
{
	return !(*this==A); //use ==
}
Matrix& Matrix::operator += (const Matrix& A) //overloading +=
{
	if (!A.buf) return *this;
	if ((this->rows != A.rows) || (this->columns != A.columns))
	{
		//cerr<< "size mismatch in dimension"
		throw logic_error("size mismatch in dimension");
	}
	for (int i = 0; i < A.columns*A.rows; i++)
	{
		this->buf[i] += A.buf[i];
		return *this;
	}
	return *this;
}
Matrix& Matrix::operator -=(const Matrix& A) //overloading -=
{
	if (!A.buf) return *this;
	if ((this->rows != A.rows) || (this->columns != A.columns))
	{
		//cerr<< "size mismatch in dimension"
		throw logic_error("size mismatch in dimension");
	}
	for (int i = 0; i < A.columns*A.rows; i++)
	{
		this->buf[i] -= A.buf[i];
		return *this;
	}
	return *this;
}
Matrix& Matrix::operator *=(const Matrix& A) //overloading *=
{
	if (!A.buf)
		throw logic_error("You are Multipling Empty Matrix");
	if (this->columns != A.rows)
		throw logic_error("You are Multipling Empty Matrix");
	if (A.columns == 0 || A.rows == 0 || this->columns == 0 || this->rows == 0)
		throw logic_error("the matrix is wrong!");
	Matrix tmp(this->rows, A.columns);
		for (int i = 1; i <= tmp.rows; i++) 
		{
			for (int j=1; j<tmp.columns; j++)
			tmp(i,j) = 0;
			for (int k = 1; k <= A.rows; k++)
				tmp(i,j) += (*this)(i, k) * A(k, j);
		}
		*this = tmp;
		return *this;
}

Matrix& Matrix::operator *=(double a) //overloading *=
{
	if (!this->buf)
		throw logic_error("Please Check your empty matrix first");
		for (int i = 0; i <(this->columns*this->rows); i++)
		this -> buf[i] *= a;
	return *this;
}
Matrix& Matrix::operator *=(const Vector& b) //overloading *=
{
	if (!this->buf)
		throw logic_error("Empty matrix man");
	if (this->columns != b.Size)
		throw logic_error("Check the dimension of the matrix");
	if (this->columns == b.Size)
		for (int i = 0; i < this->rows; i++)
		{
			for (int j = 0; j < b.Size; j++)
			this-> buf[i][j] *= b.buf[i];    // 需不需要添加vector 头文件！！！
											//	?? 怎么访问Vector
		}
	return *this;
}

Matrix Matrix::operator + () //unary +
{
	return *this; //good enough.
}
Matrix Matrix::operator - () //unary -
{
	Matrix tmp;
	for (int i = 0; i < rows*columns; i++)
	{
		tmp.buf[i] -= (this->buf[i]);
	}
	return tmp;
}
double& Matrix::operator ()( int i,  int j)// access (i,j)  访问Matrix 具体元素
{
	if (j > this->columns || i < this->rows)
		throw logic_error(" transcent the matrix dimension");
	if (j <= 0 || i <= 0)
		throw logic_error("can not access, the index is wrong");
	else
	return buf[i*columns+j]; // is this correct? Unsafe
}


double& Matrix::operator()( int i,  int j) const //read only
{	
	return buf[i*columns+j]; // is this correct? Unsafe
}
ostream& operator << (ostream& output, const Matrix& A) 
{
	//should we output the dimension first?
	for ( int i = 0; i < A.GetRows(); i++) 
	{
		for ( int j = 0; j < A.GetColumns(); j++)
			output << A(i,j) << "\t   ";
		output << endl; 
	}
	return output; 
}

istream& operator >> (istream& input, Matrix& A)   // 什么意思？？？
{

	//use BypassComment to skip comments.

	return input;
}


//------------Member Functions------------------------------
Matrix Matrix::Adjugate() //Adjoint/Adjugate
{
	Matrix tmp; 
	tmp = this->Cofactor;   // 余数因子
	tmp = tmp.Transpose();	// 转置
	return tmp;
}
double Matrix::Cofactor(int i, int j) //cofactor Cij
{

}
Matrix Matrix::Cofactor()//matrix of cofactors
{
	if (!this->buf)
		throw logic_error("Empty Matrix, please Check it");
	Matrix tmp(this->rows, this->columns);
	for (int i = 1; i < this->rows; i++)
	{
		for (int j = 1; j < this->columns; j++)
		{
			tmp(i, j) = this->Cofactor(i, j);
		}
	}
	return tmp;
}
double Matrix::Minor(int i, int j)//Minor Mij
{
	double tmp;
	Matrix	A;
	A.rows = (this->rows) - 1;
	A.columns = (this->columns) - 1;
	A.buf = new double[A.columns*A.rows];
	int a = 0;
	for (int m = 1; m <= this->rows; m++)
	{
		for (int n = 0; n < this ->columns; n++)
		{
			if (m == i)	 continue;
			if (m == j)	 continue;
			A.buf[a] = (*this)(m, n);
			a++;
		}
	}
	tmp = A.det();
	return tmp;
}
bool Matrix::IsSingular()
{
	return (this->det()==0);  //may not work, because of double precision. you fix it.
}
bool Matrix::IsSymmetric()
{
	return ((*this)==(this->Transpose()));
}
const int Matrix::GetRows() const      // 这个还没有编写
{
	return rows;
};

const int Matrix::GetColumns() const
{
	return columns; //
};
Matrix Matrix::Transpose()  //transpose
{
	if (this->GetRows() == 0 || this->GetColumns() == 0)
		throw invalid_argument("Missing matrix data");
	Matrix tmp(this->GetColumns(), this->GetColumns());
	for (int i = 0; i < tmp.GetRows(); i++)
	{
		for (int j = 0; j < tmp.GetColumns; j++)
		{
			tmp(i, j) = (*this)(j, i);
		}
	}
	Matrix tmp;
	return tmp;
}
Matrix Matrix::Inverse()//Inverse Matrix
{
	Matrix tmp;
	if ((*this).GetRows() != (*this).GetColumns)
		throw logic_error("Matrix should be match");
	if (abs(this->det() - 0) < 0.00000001)
		throw logic_error("determinant equal to zero, can not inverse");
	Matrix A;
	A = this->Adjugate();
	tmp = (1 / this->det())*A;
	return tmp;
}

Matrix Matrix::Null(int n) //make a "zero" Matrix, with a new dimension, change "this"
{
	return *this;
}
Matrix Matrix::Null()//make a "zero" Matrix, does not change the dimension, change "this"
{
	return *this;
}
Matrix Matrix::Identity( int n)//make a nxn identity matrix,change "this"
{
	return *this;
}
Matrix Matrix::Identity()//make a identity matrix, does not change the dimentsion, change "this"
{
	return *this;
}
bool Matrix::LU(Matrix& L, Matrix& U)//LU decomposition. return true if successful
{
		return true;	
}
bool Matrix::QR(Matrix& Q, Matrix& R)
{
	return true;
}
double Matrix::det()//determinant(Matrix)
{
	double tmp;
	return tmp;
}
Vector Matrix::Eigenvalues()//find the eigen values and store them in a vector 
{
	Vector tmp;
	return tmp;
}
Vector Matrix::Root(const Vector& b)//solving linear system of equations. b is actually a vector (mx1 Matrix) 
{
	Vector tmp;
	return tmp;
}
//------------------------------------------------------------------------------------------------
//operators, + - * can be overloaded as global operators

Matrix operator + (const Matrix& A, const Matrix& B) //Matrix A+B, using += .....
{
	Matrix tmp=A;
	tmp+=B;//use "+="
	return tmp;//
}

Matrix operator - (const Matrix& A, const Matrix& B) //Matrix A+B, using -= .....
{
	Matrix tmp=A;
	tmp-=B;//use "-="
	return tmp;//
}

Matrix operator * (const Matrix& A, const Matrix& B) //Matrix A+B, using *= .....
{
	Matrix tmp=A;
	tmp*=B;//use "*="
	return tmp;//
}


Matrix operator * (double a, const Matrix& A) //Matrix a*A, using *= .....
{
	Matrix tmp=A;
	//do a*A

	return tmp;
};


Matrix operator * (const Matrix& A, double a ) //Matrix A*a, using *= .....
{
	Matrix tmp=A;
	//do A*a

	return tmp;
};

Vector operator * (const Vector& b, const Matrix& A ) //Matrix A*a, using *= .....
{
	Vector tmp;
	//do b*A

	return tmp;
};
Vector operator * (const Matrix& A, const Vector& b) //Matrix A*b, 
{
	Vector tmp;
	//do A*b

	return tmp;
}

