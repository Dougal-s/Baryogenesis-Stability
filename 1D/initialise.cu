// Version 20180514-01: Fixed bug in flaton initial conditions.

#include <fstream>

#include <curand.h>
#include <curand_kernel.h>
#include <cufft.h>
#include <cufftXt.h>

#include <common/select_GPU.cu>
#include <common/complex.cu>
#include <1D/parameters.hpp>

using namespace std;

////////////////////////////
//Initialization Functions//
////////////////////////////

__global__ void PreFFTInitPhi(float2 d_PhiC[X], int f, double temperature, double msq) {

	const uint index = {blockIdx.x * blockDim.x + threadIdx.x};
	const uint stride = {blockDim.x * gridDim.x};
	
	curandState_t state;
	curand_init(index*(f+1) + 2*X*seed, 0, 0, &state);
	
	int kx;
	double qsq,r,theta;
	
	for (int x=index;x<X;x+= stride) {
	
		kx = (x+X/2-1)%X - X/2 + 1;
		
		qsq = 2*M_PI * 2*M_PI * (kx/(X*dx)) * (kx/(X*dx));
		
		if (thermal_spectrum) {	
			r = sqrt(1 / (X*dx*(exp(sqrt(qsq + msq) / temperature) - 1))); // exactly giving BE distribution
			// When given msq=0, division by zero for kx=0, eventually failing FFT
		} else {
			r = (1/sqrt(X*dx)) * ((curand(&state)%1000+0.5)/1000) / sqrt(qsq+1); // arbitrary non-thermal spectrum
		}
		theta = 2*M_PI * (curand(&state)%1000+0.5)/1000; // random phase
		d_PhiC[x].x = r*cos(theta); // Real part
		d_PhiC[x].y = r*sin(theta); // Imaginary part
		
	}
}

__global__ void PreFFTInitPhiDot(float2 d_PhiDotC[X], int f) {

	const uint index = {blockIdx.x * blockDim.x + threadIdx.x};
	const uint stride = {blockDim.x * gridDim.x};
	
	curandState_t state;
	curand_init(index*(f+1) + X + 2*X*seed, 0, 0, &state);
	
	int kx;
	double qsq,r,theta;
	
	for (int x=index;x<X;x+= stride) {
	
		kx = (x+X/2-1)%X - X/2 + 1;
		
		qsq = 2*M_PI * 2*M_PI * (kx/(X*dx)) * (kx/(X*dx));
			
		if (thermal_spectrum) {	
			r = (1/sqrt(X*dx)) * ((curand(&state)%1000+0.5)/1000) / sqrt( exp( sqrt(qsq+1) ) - 1 ); // crude thermal spectrum
		} else {
			r = (1/sqrt(X*dx)) * ((curand(&state)%1000+0.5)/1000) / sqrt(qsq+1); // arbitrary non-thermal spectrum
		}
		theta = 2*M_PI * (curand(&state)%1000+0.5)/1000; // random phase
		d_PhiDotC[x].x = r*cos(theta); // Real part
		d_PhiDotC[x].y = r*sin(theta); // Imaginary part

	}
}

__global__ void cpy(float2 d_PhiC[X], float2 d_PhiDotC[X], Complex<double> d_Phi[F][X], Complex<double> d_PhiDot[F][X], int f) 
{
	const uint index = {blockIdx.x * blockDim.x + threadIdx.x};
	const uint stride = {blockDim.x * gridDim.x};
	
	for (int x=index; x<X; x+=stride) {
		d_Phi[f][x].x = d_PhiC[x].x;
		d_Phi[f][x].y = d_PhiC[x].y;
		d_PhiDot[f][x].x = d_PhiDotC[x].x;
		d_PhiDot[f][x].y = d_PhiDotC[x].y;
	}
}

__global__ void PostFFTInit(Complex<double> d_Phi[F][X], Complex<double> d_PhiDot[F][X]) 
{
	const uint index = {blockIdx.x * blockDim.x + threadIdx.x};
	const uint stride = {blockDim.x * gridDim.x};
	
	// Initial conditions
	
	for (int x=index; x<X; x+=stride) {
			
		d_Phi[0][x] = (m0/phi0) * d_Phi[0][x];
		d_PhiDot[0][x] = (m0/phi0) * d_PhiDot[0][x];

     	for (short int f=1; f<5; f++) {
			 d_Phi[f][x] = (m0/l0) * d_Phi[f][x];
			 d_PhiDot[f][x] = (m0/l0) * d_PhiDot[f][x];
		}
			
		d_Phi[3][x].x += cos((M_PI-argAnu)/2.0); // real l
		d_Phi[3][x].y += sin((M_PI-argAnu)/2.0); // imag l
	
		d_Phi[1][x].x = sqrt( modsq(d_Phi[2][x]) + modsq(d_Phi[3][x]) - 0.5*modsq(d_Phi[4][x]) + epsilon*epsilon ); // real h_u
		d_Phi[1][x].y = 0.0; // imag h_u
		
	}
}


void print(Complex<double> Phi[F][X], Complex<double> PhiDot[F][X], float ndt, ofstream &PhiFile, ofstream &PhiDotFile) {

	Complex<float> tmp[F];

	for (uint32_t x=0; x<X; x++) {
		PhiFile.write( reinterpret_cast<const char*>(&ndt), sizeof(ndt) );
		PhiFile.write( reinterpret_cast<const char*>(&x), sizeof(x) );
		for (int f=0; f<F; f++) {
			tmp[f].x = Phi[f][x].x;
			tmp[f].y = Phi[f][x].y;
		}
		PhiFile.write( reinterpret_cast<const char*>(&tmp), sizeof(tmp) );
	}

	for (uint32_t x=0; x<X; x++) {
		PhiDotFile.write( reinterpret_cast<const char*>(&ndt), sizeof(ndt) );
		PhiDotFile.write( reinterpret_cast<const char*>(&x), sizeof(x) );
		for (int f=0; f<F; f++) {
			tmp[f].x = PhiDot[f][x].x;
			tmp[f].y = PhiDot[f][x].y;
		}
		PhiDotFile.write( reinterpret_cast<const char*>(&tmp), sizeof(tmp) );
	}
}

Complex<double> Phi[F][X];
Complex<double> PhiDot[F][X];

int main() {
	cout << "\n////////////////////////////////////////////////////////////////////////////////\n\n";
		
		SelectDevice();
		
		Complex<double> (*d_Phi)[X];
		Complex<double> (*d_PhiDot)[X];
		
		cudaMalloc(&d_Phi, (F*X)*sizeof(Complex<double>));
		cudaMalloc(&d_PhiDot, (F*X)*sizeof(Complex<double>));
		
		float2 (*d_PhiC);
		float2 (*d_PhiDotC);
		cudaMalloc(&d_PhiC, (X)*sizeof(float2));
		cudaMalloc(&d_PhiDotC, (X)*sizeof(float2));
		
		dim3 block = {1024};
		dim3 grid = {(X+block.x-1)/block.x};
		
		cudaError_t __err;
		
		cout << "\n=============================================\n";
		cout << "||                                         ||";

		cout << "\n|| Initializing                            ||" << endl;
		cout << "||                                         ||"; 
		
		__err = cudaGetLastError();
		
		if (__err != cudaSuccess) {
			cout << "||                                         ||\n";
			cout << "|| Failed to create variables              ||\n";
			cout << "|| "  << setw(40) << left << cudaGetErrorString(__err) << "||\n";
			cout << "=============================================\n";
			return -1;
		}
		
		//////////////////
		//Initialization//
		//////////////////
		
		// FFT initialitation
		
		cout << "\n|| Running inverse fast fourier transform  ||\n";
		cout << "||                                         ||\n";
		cufftHandle plan;
		cufftPlan1d( &plan, X, CUFFT_C2C, 1);
		
		for (short int f=0; f<5; f++) {
		
			// Pre FFT initialization
			
			PreFFTInitPhi<<<grid, block>>>(d_PhiC, f, 1, 1); // TODO: make T & msq as external parameters
			PreFFTInitPhiDot<<<grid, block>>>(d_PhiDotC, f);
			cudaDeviceSynchronize();
			
			cufftExecC2C( plan, (cufftComplex *) d_PhiC, (cufftComplex *) d_PhiC, CUFFT_INVERSE );
			cufftExecC2C( plan, (cufftComplex *) d_PhiDotC, (cufftComplex *) d_PhiDotC, CUFFT_INVERSE );
			
			cpy<<<grid, block>>>(d_PhiC, d_PhiDotC, d_Phi, d_PhiDot, f);
			cudaDeviceSynchronize();
		}
	
		for (short int f=5; f<F; f++) {
		
			// Pre FFT initialization
			
			PreFFTInitPhi<<<grid, block>>>(d_PhiC, f, 1, 0); // msq=0 is causing problem
			PreFFTInitPhiDot<<<grid, block>>>(d_PhiDotC, f);
			cudaDeviceSynchronize();
			
			cufftExecC2C( plan, (cufftComplex *) d_PhiC, (cufftComplex *) d_PhiC, CUFFT_INVERSE );
			cufftExecC2C( plan, (cufftComplex *) d_PhiDotC, (cufftComplex *) d_PhiDotC, CUFFT_INVERSE );
			
			cpy<<<grid, block>>>(d_PhiC, d_PhiDotC, d_Phi, d_PhiDot, f);
			cudaDeviceSynchronize();
		}
		
		cufftDestroy(plan);
		cudaFree(d_PhiC);
		cudaFree(d_PhiDotC);
		
		// Post FFT initialization
		
		PostFFTInit<<<grid, block>>>( d_Phi, d_PhiDot);
		cudaDeviceSynchronize();
		
		cudaMemcpy(Phi, d_Phi, (F*X)*sizeof(Complex<double>), cudaMemcpyDeviceToHost);
		cudaMemcpy(PhiDot, d_PhiDot, (F*X)*sizeof(Complex<double>), cudaMemcpyDeviceToHost);
		
		// Printing data
		
		ofstream PhiFile("initialPhi.bin", ios::binary);
		ofstream PhiDotFile("initialPhiDot.bin", ios::binary);
		
		print(Phi, PhiDot, 0, PhiFile, PhiDotFile);
		
		PhiFile.close();
		PhiDotFile.close();
		
		cudaFree(d_Phi);
		cudaFree(d_PhiDot);
		
		__err = cudaGetLastError();
		
		if (__err != cudaSuccess) {
			cout << "||                                         ||\n";
			cout << "|| Failed to create variables              ||\n";
			cout << "|| "  << setw(40) << left << cudaGetErrorString(__err) << "||\n";
			cout << "=============================================\n";
			return -1;
		}
		
		cout << "|| Completed Initializing                  ||\n";
		cout << "=============================================\n";
		
		cout << "\n\n////////////////////////////////////////////////////////////////////////////////\n" << endl;
		
		return 0;
}
