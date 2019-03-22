// Version 20180208-02: Constant multiplication now works.

#include <math.h>
#include <stdio.h>

template <class num>
class Complex {
	public:
	
	// Variables

	num x;
	num y;

	//Functions

	// Addition
	
	__host__ __device__ Complex operator+(Complex rhs) const {
		return Complex{x + rhs.x, y + rhs.y};
	}

	__host__ __device__ Complex operator+(num rhs) const {
		return Complex{x + rhs, y}; 
	}

	__host__ __device__ friend Complex operator+(num lhs, Complex rhs) {
		return Complex{lhs + rhs.x, rhs.y};
	}
 
	// Subtraction

	__host__ __device__ Complex operator-(Complex rhs) const {
		return Complex{x - rhs.x, y - rhs.y}; 
	}

	__host__ __device__ Complex operator-(num rhs) const {
		return Complex{x - rhs, y}; 
	}
 
	__host__ __device__ friend Complex operator-(num lhs, Complex rhs) {
		return Complex{lhs - rhs.x, -rhs.y}; 
	}

	// Multiplication

	__host__ __device__ Complex operator*(Complex rhs) const {
		return Complex{x * rhs.x - y * rhs.y, x * rhs.y + y * rhs.x}; 
	}

	__host__ __device__ Complex operator*(num rhs) const {
		return Complex{x * rhs, y * rhs}; 
	}

	__host__ __device__ friend Complex operator*(num lhs, Complex rhs) {
		return Complex{lhs * rhs.x, lhs * rhs.y};
	}
 
	// Division

	__host__ __device__ Complex operator/(Complex rhs) const { 
		return Complex{(x * rhs.x + y * rhs.y)/(rhs.x * rhs.x + rhs.y * rhs.y), (- x * rhs.y + y * rhs.x)/(rhs.x * rhs.x + rhs.y * rhs.y)}; 
	}

	__host__ __device__ Complex operator/(num rhs) const { 
		return Complex{x / rhs, y / rhs}; 
	}

	__host__ __device__ friend Complex operator/(num lhs, Complex rhs) {
		return Complex{lhs * rhs.x/(rhs.x * rhs.x + rhs.y * rhs.y), - lhs * rhs.y/(rhs.x * rhs.x + rhs.y * rhs.y)};
	}
};


//Functions

template <class num>
__host__ __device__ Complex<num> conj(Complex<num> a) {
	return Complex<num>{a.x,-a.y};
}

template <class num>
__host__ __device__ num modsq(Complex<num> a) {
	return a.x*a.x + a.y*a.y;
}
	
template <class num>
__host__ __device__ num abs(Complex<num> a) {
	return sqrtf( a.x*a.x + a.y*a.y );
}
	
template <class num>
__host__ __device__ num arg(Complex<num> a) {

	if (a.x > 0) {
		return atan(a.y/a.x);
	}
	if (a.x < 0 && a.y >= 0) {
		return M_PI + atan(a.y/a.x); 
	} 
	if (a.x < 0 && a.y < 0) {
		return -M_PI + atan(a.y/a.x); 
	}
	if (a.y > 0) {
		return M_PI/2;
	}
	if (a.y < 0) {
		return -M_PI/2;
	}
	
	printf("\nWARNING: arg(0)\n");
	return 0;
	
}
