// Version 20180514-01: Eliminated flaton and AD damping. Eliminated option to turn off back-reaction.

// Libraries
#include <fstream>
#include <string>
#include <iostream>
#include <chrono>

#include <common/select_GPU.cu>
#include <common/complex.cu>
#include <1D/parameters.hpp>

using namespace std;

// field macros

#define phi d_Phi[0][x]
#define hu d_Phi[1][x]
#define hd d_Phi[2][x]
#define l d_Phi[3][x]
#define d d_Phi[4][x]

////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions
////////////////////////////////////////////////////////////////////////////////////////////////////

// read initial data from input files

void Read(Complex<double> Phi[F][X], Complex<double> PhiDot[F][X], string inputPhiStr, string inputPhiDotStr, double ndt) {

	// Variables

	float n;
	uint32_t x;
	Complex<float> tmp[F];

	// read Phi

	ifstream inputPhi(inputPhiStr, std::ios::binary);

	while (inputPhi.read( reinterpret_cast<char*>(&n), sizeof(n))) {
		if (n != ndt) {
			inputPhi.seekg((sizeof(x) + sizeof(tmp)), std::ios::cur);
			continue;
		}
		inputPhi.read( reinterpret_cast<char*>(&x), sizeof(x));
		inputPhi.read( reinterpret_cast<char*>(&tmp), sizeof(tmp));
		for (int f = 0; f < F; f++)
			Phi[f][x] = Complex<double>{tmp[f].x, tmp[f].y};
	}

	inputPhi.close();

	// read PhiDot

	ifstream inputPhiDot(inputPhiDotStr, std::ios::binary);

	while (inputPhiDot.read( reinterpret_cast<char*>(&n), sizeof(n))) {
		if (n != ndt) {
			inputPhiDot.seekg( (sizeof(x) + sizeof(tmp)), std::ios::beg);
			continue;
		}
		inputPhiDot.read( reinterpret_cast<char*>(&x), sizeof(x));
		inputPhiDot.read( reinterpret_cast<char*>(&tmp), sizeof(tmp));
		for (int f = 0; f < F; f++)
			PhiDot[f][x] = Complex<double>{tmp[f].x, tmp[f].y};
	}

	inputPhiDot.close();
}


// Calculate a single timestep of Phi

__global__ void StepPhi(Complex<double> d_Phi[F][X], Complex<double> d_PhiDot[F][X], double dt)
{
	const uint2 index = {blockIdx.x * blockDim.x + threadIdx.x, blockIdx.y * blockDim.y + threadIdx.y};
	const uint2 stride = {blockDim.x * gridDim.x, blockDim.y * gridDim.y};

	for (int x = index.x; x<X; x += stride.x) {
		for (int f = index.y; f<F; f += stride.y) {
			d_Phi[f][x] = d_Phi[f][x] + dt * d_PhiDot[f][x];
		}
	}
}


// Calculate a single timestep of PhiDot

__global__ void StepPhiDot(Complex<double> d_Phi[F][X], Complex<double> d_PhiDot[F][X], int n, double dt)
{
	const uint index = {blockIdx.x * blockDim.x + threadIdx.x};
	const uint stride = {blockDim.x * gridDim.x};

	// charges

	const double Q[13] = {0,1,-1,-1,0.5,0,0,0,0,0,0,0,0}; // phi,h_u,h_d,l,d,8xg

	// dependent parameters

    const Complex<double> Amu = Complex<double>{ modAmu*cosf(argAmu), modAmu*sinf(argAmu)};
    const Complex<double> Anu = Complex<double>{ modAnu*cosf(argAnu), modAnu*sinf(argAnu)};
    const Complex<double> Ad = Complex<double>{ modAd*cosf(argAd), modAd*sinf(argAd)};

	const double kmu = mu; // phi_0 = 1
	const double knu = (( sqrt(modsq(Anu)) + sqrt( modsq(Anu) - 6.0*(mLsq+mHusq) ) ) / 6.0); // l_0 = 1

	// evolution variables

	const double Tphi = Tphi_i * exp(-GammaTphi*n*time_step);

	Complex<double> d2Phidx2[F];
	double Tsq;
	Complex<double> dVdPhi[F];
	double PhiQQPhi;
	double PhiQQQPhi;
	Complex<double> PhiQPhidot;
	Complex<double> PhiQd2Phidx2;
	Complex<double> PhiQdVdPhi;
	Complex<double> pi[F];
	double piQpi;
	Complex<double> PhiQQpi;
	double Pi;

	for (int x=index; x<X; x+=stride) {

		// evaluate Laplacian

		for (short int f=0; f<13; f++) {
			d2Phidx2[f] = ( d_Phi[f][(x+1)%X] - 2*d_Phi[f][x] + d_Phi[f][(x+X-1)%X] ) / (dx*dx);
		}

		/* potential

		V = V0 - Tphi*Tphi*Tphi*Tphi * (m0*m0/(phi0*phi0)) * exp(-conj(phi)*phi*phi0*phi0/(Tphi*Tphi*m0*m0)) + 0.5*alpha*mphisq * ( log(conj(phi)*phi+mssq) - 1 ) * conj(phi)*phi + mHusq*conj(hu)*hu + mHdsq*conj(hd)*hd + mLsq*conj(l)*l + mdsq*conj(d)*d - Amu*kmu*phi*phi*hu*hd + 0.5*Anu*knu*l*l*hu*hu - 0.5*Ad*ld*hd*d*d - conj(Amu*kmu*phi*phi*hu*hd) + conj(0.5*Anu*knu*l*l*hu*hu) - conj(0.5*Ad*ld*hd*d*d) + conj(kmu*phi*phi*hu+0.5*ld*d*d)*(kmu*phi*phi*hu+0.5*ld*d*d) + conj(kmu*phi*phi*hd-knu*l*l*hu)*(kmu*phi*phi*hd-knu*l*l*hu) + conj(knu*l*hu*hu)*(knu*l*hu*hu) + conj(ld*hd*d)*(ld*hd*d) + conj(ldg*g*d)*(ldg*g*d) + conj(g*g)*(g*g);

		*/

		// evaluate potential derivatives

		Tsq = 0;
		for (short int f=5; f<13; f++) {
			Tsq += modsq(d_Phi[f][x]);
		}

		dVdPhi[0] = Tphi * Tphi * phi * exp(-modsq(phi)*phi0*phi0/(Tphi*Tphi*m0*m0)) + 0.5*alpha*mphisq * ( log(modsq(phi)+mssq) - mssq/(modsq(phi)+mssq) ) * phi + 2.0*(l0/phi0)*(l0/phi0) * ( 0.0 - conj(Amu*kmu*phi*hu*hd) + conj(kmu*phi*hu)*(kmu*phi*phi*hu+0.5*ld*d*d) + conj(kmu*phi*hd)*(kmu*phi*phi*hd-knu*l*l*hu) );

		dVdPhi[1] = mHusq*hu - conj(Amu*kmu*phi*phi*hd) + conj(Anu*knu*l*l*hu) + conj(kmu*phi*phi)*(kmu*phi*phi*hu+0.5*ld*d*d) - conj(knu*l*l)*(kmu*phi*phi*hd-knu*l*l*hu) + 2.0*conj(knu*l*hu)*(knu*l*hu*hu);

		dVdPhi[2] = mHdsq*hd - conj(Amu*kmu*phi*phi*hu) - 0.5*conj(Ad*ld*d*d) + conj(kmu*phi*phi)*(kmu*phi*phi*hd-knu*l*l*hu) + conj(ld*d)*ld*hd*d;

		dVdPhi[3] = mLsq*l + conj(Anu*knu*l*hu*hu) - 2.0*conj(knu*l*hu)*(kmu*phi*phi*hd-knu*l*l*hu) + conj(knu*hu*hu)*knu*l*hu*hu;

		dVdPhi[4] = mdsq*d - conj(Ad*ld*hd*d) + conj(ld*d)*(kmu*phi*phi*hu+0.5*ld*d*d) + conj(ld*hd)*ld*hd*d + ldg*ldg*Tsq*d;

		for (short int f=5;f<13;f++) {
			dVdPhi[f] = conj(ldg*(l0*d/m0))*ldg*(l0*d/m0)*d_Phi[f][x] + conj(d_Phi[f][x])*d_Phi[f][x]*d_Phi[f][x];
		}

		// evaluate evolution variables

		PhiQQPhi = 0;
		for (short int f=1; f<5; f++) {
			PhiQQPhi += Q[f] * Q[f] * modsq(d_Phi[f][x]);
		}

		PhiQQQPhi = 0;
		for (short int f=1; f<5; f++) {
			PhiQQQPhi += Q[f] * Q[f] * Q[f] * modsq(d_Phi[f][x]);
		}

		PhiQPhidot = Complex<double>{0,0};
		for (short int f=1; f<5; f++) {
			PhiQPhidot = PhiQPhidot + conj(d_Phi[f][x]) * Q[f] * d_PhiDot[f][x];
		}

		PhiQd2Phidx2 = Complex<double>{0,0};
		for (short int f=1; f<5; f++) {
			PhiQd2Phidx2 = PhiQd2Phidx2 + conj(d_Phi[f][x]) * Q[f] * d2Phidx2[f];
		}

		PhiQdVdPhi = Complex<double>{0,0};
		for (short int f=1; f<5; f++) {
			PhiQdVdPhi = PhiQdVdPhi + conj(d_Phi[f][x]) * Q[f] * dVdPhi[f];
		}

		for (short int f=1; f<5; f++) {
			pi[f] = d_PhiDot[f][x] + dt * d2Phidx2[f] -  dt * dVdPhi[f] - ( ( PhiQPhidot + dt * PhiQd2Phidx2 - dt * PhiQdVdPhi ) / PhiQQPhi ) * Q[f] * d_Phi[f][x];
		}

		piQpi = 0;
		for (short int f=1; f<5; f++) {
			piQpi += Q[f] * modsq(pi[f]);
		}

		PhiQQpi = Complex<double>{0,0};
		for (short int f=1; f<5; f++) {
			PhiQQpi = PhiQQpi + conj(d_Phi[f][x]) * Q[f] * Q[f] * pi[f];
		}

		Pi = ( 2.0 * piQpi ) / ( PhiQQPhi + dt * PhiQQpi.x + sqrt( ( PhiQQPhi + dt * PhiQQpi.x ) * ( PhiQQPhi + dt * PhiQQpi.x ) - dt * dt * piQpi * PhiQQQPhi ) );

		// step

		d_PhiDot[0][x] = d_PhiDot[0][x] + dt * d2Phidx2[0] - dt * dVdPhi[0];

		for (short int f=1; f<5; f++) {
			d_PhiDot[f][x] = pi[f] - 0.5 * dt * Pi * Q[f] * d_Phi[f][x];
		}

		for (short int f=5; f<13; f++) {
			d_PhiDot[f][x] = d_PhiDot[f][x] + dt * d2Phidx2[f] - dt * dVdPhi[f];
		}

	}
}


// write continuous data to output file

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


////////////////////////////////////////////////////////////////////////////////////////////////////
//Main Loop
////////////////////////////////////////////////////////////////////////////////////////////////////

Complex<double> Phi[F][X];
Complex<double> PhiDot[F][X];

int main () {

	/////////////
	//Variables//
	/////////////

	chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();

	cout << "\n////////////////////////////////////////////////////////////////////////////////\n\n";

	SelectDevice();

	string inputPhi, inputPhiDot, fileName;

	fileName = "Phi.bin";
	ofstream PhiFile(fileName, ios::binary);

	fileName = "PhiDot.bin";
	ofstream PhiDotFile(fileName, ios::binary);

	Complex<double> (*d_Phi)[X];
	Complex<double> (*d_PhiDot)[X];

	cudaMalloc(&d_Phi, (F*X)*sizeof(Complex<double>));
	cudaMalloc(&d_PhiDot, (F*X)*sizeof(Complex<double>));

	const dim3 blockL = {1024};
	const dim3 gridL = {(X+blockL.x-1)/blockL.x};

	const dim3 blockS = {256};
	const dim3 gridS = {(X+blockS.x-1)/blockS.x};

	cudaError_t __err;

	__err = cudaGetLastError();

	if (__err != cudaSuccess) {
		cout << "\nFailed to create variables" << endl;
		cout << cudaGetErrorString(__err) << endl;
		return -1;
	}

	//////////////////
	//Initialization//
	//////////////////

	cout << "Would you like to use " << "initialPhi.txt" << " and " << "initialPhiDot.txt" << "?\n";
	cout << "(Y)es/(N)o:";
	//cin >> inputPhi;

	inputPhi = "y";

	if (inputPhi == "Y" || inputPhi == "y" || inputPhi == "Yes"  || inputPhi == "yes") {
		inputPhi = "initialPhi.bin";
		inputPhiDot = "initialPhiDot.bin";
	} else {
		cout << "\nEnter the name of the input file for Phi:";
		cin >> inputPhi;

		cout << "\nEnter the name of the input file for PhiDot:";
		cin >> inputPhiDot;
	}

	cout << "\n=============================================\n";
	cout << "||                                         ||\n";

	Read(Phi, PhiDot, inputPhi, inputPhiDot, N_i*time_step);

	print(Phi, PhiDot, 0*time_step, PhiFile, PhiDotFile);

	cout << "|| Finished reading files                  ||\n";

	cudaMemcpy(d_Phi, Phi, (F*X)*sizeof(Complex<double>), cudaMemcpyHostToDevice);
	cudaMemcpy(d_PhiDot, PhiDot, (F*X)*sizeof(Complex<double>), cudaMemcpyHostToDevice);

	__err = cudaGetLastError();

	if (__err != cudaSuccess) {
		cout << "||                                         ||\n";
		cout << "|| Failed to initialize                    ||" << endl;
		cout << "|| "  << setw(40) << left << cudaGetErrorString(__err) << "||" << endl;
		cout << "=============================================\n";
		return -1;
	}

	cout << "|| Starting Simulation                     ||" << endl;

	///////////////
	//Calculation//
	///////////////

	for (int n=N_i; n<N_f; n+=S) {

		StepPhiDot<<<gridS, blockS>>>(d_Phi, d_PhiDot, n, time_step/2);
		cudaDeviceSynchronize();

		for (int s = 1; s < S; s++) {
			StepPhi<<<gridL, blockL>>>(d_Phi, d_PhiDot, time_step);
			cudaDeviceSynchronize();

			StepPhiDot<<<gridS, blockS>>>(d_Phi, d_PhiDot, n+s, time_step);
			cudaDeviceSynchronize();
		}

		StepPhi<<<gridL, blockL>>>(d_Phi, d_PhiDot, time_step);
		cudaDeviceSynchronize();

		StepPhiDot<<<gridS, blockS>>>(d_Phi, d_PhiDot, n+S, time_step/2);
		cudaDeviceSynchronize();

		// Print continuous data

		cudaMemcpy(Phi, d_Phi, (F*X)*sizeof(Complex<double>), cudaMemcpyDeviceToHost);
		cudaMemcpy(PhiDot, d_PhiDot, (F*X)*sizeof(Complex<double>), cudaMemcpyDeviceToHost);

		print(Phi, PhiDot, (n+S)*time_step, PhiFile, PhiDotFile);

		cout << "|| Completed timestep "  << setw(20) << left << n+S << " ||" << endl;

	}

	PhiFile.close();
	PhiDotFile.close();

	///////////
	//Cleanup//
	///////////

	cout << "||                                         ||";

	cudaFree(d_Phi);
	cudaFree(d_PhiDot);

	__err = cudaGetLastError();

	if (__err == cudaSuccess) {
		cout << "\n|| Everything seems fine :D                ||" << endl;
	} else {
		cout << "\n|| Something is wrong... :(                ||" << endl;
		cout << "|| "  << setw(40) << left << cudaGetErrorString(__err) << "||" << endl;
	}


	cout << "||                                         ||";
	cout << "\n=============================================\n";

	chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();

	cout << "\nProgram took " << chrono::duration_cast<chrono::seconds>(end - start).count() << " seconds.\n";

	cout << "\n\n////////////////////////////////////////////////////////////////////////////////\n" << endl;

}
