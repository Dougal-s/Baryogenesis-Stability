// Version 20180504-01: Error fixed.

#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <chrono>

#include <1D/parameters.hpp>

const short int E = 4+2*F+2; // number of energies (E, K,G,V, K_phi,K_hu,K_hd,K_l,K_d,K_g1,...,K_g8, G_phi,G_hu,G_hd,G_l,G_d,G_g1,...,G_g8, V_phi,V_AD)

// dependent parameters

const std::complex<double> Amu = std::complex<double>{ modAmu*cos(argAmu), modAmu*sin(argAmu)};
const std::complex<double> Anu = std::complex<double>{ modAnu*cos(argAnu), modAnu*sin(argAnu)};
const std::complex<double> Ad = std::complex<double>{ modAd*cos(argAd), modAd*sin(argAd)};

const double kmu = mu; // phi_0 = 1
const double knu = (( sqrt(norm(Anu)) + sqrt( norm(Anu) - 6.0*(mLsq+mHusq) ) ) / 6.0); // l_0 = 1

// field macros

#define phi Phi[0][x]
#define hu Phi[1][x]
#define hd Phi[2][x]
#define l Phi[3][x]
#define d Phi[4][x]

void Read(std::complex<double> Phi[F][X], std::complex<double> PhiDot[F][X], std::string inputPhiStr, std::string inputPhiDotStr, unsigned int timestep, uint *start) {

	// Variables

	float n;
	uint32_t x;
	bool read = false;
	std::complex<float> tmp[F];

	// read Phi

	std::ifstream inputPhi(inputPhiStr, std::ios::binary);

	inputPhi.seekg(*start, std::ios::beg);

	while (inputPhi.read( reinterpret_cast<char*>(&n), sizeof(n))) {
		if (n != timestep) {
			if (read) {
				break;
			}
			inputPhi.seekg((sizeof(x) + sizeof(tmp)), std::ios::cur);
			continue;
		}
		read = true;
		inputPhi.read( reinterpret_cast<char*>(&x), sizeof(x));
		inputPhi.read( reinterpret_cast<char*>(&tmp), sizeof(tmp));
		for (int f = 0; f < F; f++)
			Phi[f][x] = std::complex<double>{tmp[f].real(), tmp[f].imag()};
	}

	inputPhi.close();

	read = false;

	// read PhiDot

	std::ifstream inputPhiDot(inputPhiDotStr, std::ios::binary);

	inputPhiDot.seekg(*start, std::ios::beg);

	while (inputPhiDot.read( reinterpret_cast<char*>(&n), sizeof(n))) {
		if (n != timestep) {
			if (read) {
				break;
			}
			inputPhiDot.seekg( (sizeof(x) + sizeof(tmp)), std::ios::beg);
			continue;
		}

		read = true;
		inputPhiDot.read( reinterpret_cast<char*>(&x), sizeof(x));
		inputPhiDot.read( reinterpret_cast<char*>(&tmp), sizeof(tmp));
		for (int f = 0; f < F; f++)
			PhiDot[f][x] = std::complex<double>{tmp[f].real(), tmp[f].imag()};
	}

	inputPhiDot.seekg(-sizeof(n) , std::ios::cur);

	*start = inputPhiDot.tellg();

	inputPhiDot.close();
}

void Print(std::ofstream &output, double energies[E], int n, double dt) {

	output << n*dt;

	for (int x = 0; x < E; x++) {
		output << ' ' << energies[x];
	}

	output << std::endl;

}

void CalcEnergies(double energies[E], std::complex<double> Phi[F][X], std::complex<double> PhiDot[F][X]) {
	using namespace std;

	// E, K,G,V, K_phi,K_hu,K_hd,K_l,K_d,K_g1,...,K_g8, G_phi,G_hu,G_hd,G_l,G_d,G_g1,...,G_g8, V_phi,V_AD

	for (int e = 0; e < E; e++) {
		energies[e] = 0;
	}

	for (int f = 0; f < F; f++) {
		for (int x = 0; x < X; x++) {
			energies[4+f] += norm(PhiDot[f][x]); // field kinetic energy
			energies[4+F+f] += norm( Phi[f][(x+1)%X] - Phi[f][x] ) / (dx*dx); // field gradient energy
		}
	}

	for (int x = 0; x < X; x++) {
		energies[4+2*F+0] += 0.5*alpha*mphisq * ( log(norm(phi)+mssq) - 1 ) * norm(phi); // flaton potential energy
		energies[4+2*F+1] += mHusq*norm(hu) + mHdsq*norm(hd) + mLsq*norm(l) + mdsq*norm(d) - 2*real(Amu*kmu*phi*phi*hu*hd) + real(Anu*knu*l*l*hu*hu) - real(Ad*ld*hd*d*d) + norm(kmu*phi*phi*hu+0.5*ld*d*d) + norm(kmu*phi*phi*hd-knu*l*l*hu) + norm(knu*l*hu*hu) + norm(ld*hd*d);
		for (int f = 5; f < 13; f++) {
			energies[4+2*F+1] += (m0/l0)*(m0/l0)*norm(ldg*Phi[f][x]*d) + (m0/l0)*(m0/l0)*(m0/l0)*(m0/l0)*norm(Phi[f][x]*Phi[f][x]); // AD potential energy
		}
	}

	// convert to energy densities
	for (int e = 0; e < E; e++) {
		energies[e] /= X;
	}

	// normalise
	for (int f = 1; f < 5; f++) {
		energies[4+f] *= (l0/phi0)*(l0/phi0); // AD kinetic energy density
		energies[4+F+f] *= (l0/phi0)*(l0/phi0); // AD gradient energy density
	}
	energies[4+2*F+1] *= (l0/phi0)*(l0/phi0); // AD potential energy density
	for (int f = 5; f < F; f++) {
		energies[4+f] *= (m0/phi0)*(m0/phi0); // gluon kinetic energy density
		energies[4+F+f] *= (m0/phi0)*(m0/phi0); // gluon gradient energy density
	}

	// sum over fields
	for (int f = 0; f < F; f++) {
		energies[1] += energies[4+f]; // kinetic energy
		energies[2] += energies[4+F+f]; // gradient energy
	}
	energies[3] = energies[4+2*F+0] + energies[4+2*F+1]; // potential energy

	energies[0] = energies[1] + energies[2] + energies[3]; // total energy

}

///////////////////////////////////////////
//Main Function
///////////////////////////////////////////

std::complex<double> Phi[F][X];
std::complex<double> PhiDot[F][X];
double energies[E];

int main() {

	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

	std::ofstream energiesOutput("energies.txt");

	uint fileStart = 0;

	for (int n = N_i; n <= N_f; n+=S) {
		Read(Phi, PhiDot, "./DataFiles/Phi.bin", "./DataFiles/PhiDot.bin", n*time_step, &fileStart);
		CalcEnergies(energies, Phi, PhiDot);
		Print(energiesOutput, energies, n, time_step);
		std::cout << n << std::endl;
	}

	energiesOutput.close();

	std::cout << "Program took " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start).count() << " seconds\n";

	return 0;
}
