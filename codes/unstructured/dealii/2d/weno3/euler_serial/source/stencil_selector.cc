#include "../include/stencil_selector.h" 

// Default constructor 

stencil_selector::stencil_selector() {
	
	// Set all to the worst possible case 
	
	is_N_present = false;
	N_index = -1; 
	
	is_S_present = false;
	S_index = -1;
	
	is_W_present = false;
	W_index = -1;
	
	is_E_present = false;
	E_index = -1; 
	
	NE_index.resize(0); 
	SW_index.resize(0);
	SE_index.resize(0);
	NW_index.resize(0);
}

// Print information about the stencil

void stencil_selector::print_stencil_info() const {
	
	
	std::cout << "--------------------- Face neighbor info ---------------------" << "\n"; 
	
	if (is_W_present) {
		std::cout << "W cell is present" << "\n"; 
		std::cout << "cell index = " << W_index << "\n"; 
	}
	
	else {
		std::cout << "W cell is absent" << "\n"; 
	}
	
	if (is_E_present) {
		std::cout << "E cell is present" << "\n"; 
		std::cout << "cell index = " << E_index << "\n"; 
	}
	
	else {
		std::cout << "E cell is absent" << "\n"; 
	}
	
	if (is_S_present) {
		std::cout << "S cell is present" << "\n"; 
		std::cout << "cell index = " << S_index << "\n"; 
	}
	
	else {
		std::cout << "S cell is absent" << "\n"; 
	}
	
	if (is_N_present) {
		std::cout << "N cell is present" << "\n"; 
		std::cout << "cell index = " << N_index << "\n"; 
	}
	
	else {
		std::cout << "N cell is absent" << "\n"; 
	}
	
	std::cout << "--------------------- Vertex neighbor info ---------------------" << "\n";
	
	std::cout << "No. of SW cells = " << SW_index.size() << "\n"; 
	
	for (unsigned int i = 0; i < SW_index.size(); ++i) {
		std::cout << "  " << SW_index[i] << "\n"; 
	}
	
	std::cout << "No. of SE cells = " << SE_index.size() << "\n"; 
	
	for (unsigned int i = 0; i < SE_index.size(); ++i) {
		std::cout << "  " << SE_index[i] << "\n"; 
	}
	
	std::cout << "No. of NW cells = " << NW_index.size() << "\n"; 
	
	for (unsigned int i = 0; i < NW_index.size(); ++i) {
		std::cout << "  " << NW_index[i] << "\n"; 
	}
	
	std::cout << "No. of NE cells = " << NE_index.size() << "\n"; 
	
	for (unsigned int i = 0; i < NE_index.size(); ++i) {
		std::cout << "  " << NE_index[i] << "\n"; 
	}
}
