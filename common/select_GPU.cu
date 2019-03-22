// Version 20180131-01: added version number.

#include <iostream>
#include <iomanip>

void SelectDevice () {
	int nDevices;
	int selectedDevice = 0;
	cudaGetDeviceCount(&nDevices);
	
	if (nDevices == 1) {
		cudaSetDevice(selectedDevice);
	} else {
		std::cout << "\n=============================================\n";
		std::cout << "||                                         ||";
		std::cout << "\n|| There are " << nDevices << " CUDA compatible devices.    ||\n";
		std::cout << "||                                         ||\n";
		
		for (int i = 0; i < nDevices; i++) {
			cudaDeviceProp prop;
			cudaGetDeviceProperties(&prop, i);
			std::cout << "||=========================================||\n";
			std::cout << "||                                         ||\n";
			std::cout << "|| Device Number: " << std::setw(25) << std::left << i <<"||\n";
			std::cout << "||   Device name: " << std::setw(25) << std::left << prop.name << "||\n";
			std::cout << "||   Memory Clock Rate (MHz): " << std::setw(13) << std::left << prop.memoryClockRate/1000 << "||\n";
			std::cout << "||   Memory Bus Width (bits): " << std::setw(13) << std::left << prop.memoryBusWidth << "||\n";
			std::cout << "||   Peak Memory Bandwidth (GB/s): " << std::setw(7) << std::left << 2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6 << " ||\n";
			std::cout << "||                                         ||\n";
		}
		
		std::cout << "=============================================\n";
		
		std::cout << "\nPlease enter a Device Number: ";
		std::cin >> selectedDevice;
		
		cudaSetDevice(selectedDevice);
	}
}