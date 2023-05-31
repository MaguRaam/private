/*
 * stencil_selector.h
 *      Author: sunder
 */

#ifndef STENCIL_SELECTOR_H_
#define STENCIL_SELECTOR_H_  

#include "Headers.h"

struct stencil_selector {
	bool is_N_present;
	int N_index; 
	
	bool is_S_present;
	int S_index;
	
	bool is_W_present;
	int W_index;
	
	bool is_E_present;
	int E_index; 
	
	std::vector<unsigned int> NE_index; 
	std::vector<unsigned int> SW_index;
	std::vector<unsigned int> SE_index;
	std::vector<unsigned int> NW_index;
	
	bool is_NN_present;
	int NN_index; 
	
	bool is_SS_present;
	int SS_index;
	
	bool is_WW_present;
	int WW_index;
	
	bool is_EE_present;
	int EE_index; 
	
	stencil_selector();  // Default constructor 
	void print_stencil_info() const;
}; 

#endif /* STENCIL_SELECTOR_H_ */
