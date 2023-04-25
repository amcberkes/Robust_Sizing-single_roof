#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include "simulate_system.h"
#include "common.h"
#include "cheby.h"

using namespace std;

// run_simulations
// load_filename: filename, each line in file contains electricity consumption value
// solar_filename: filename, each line in file contains solar generation value
// id: request id
// metric: 0 for LOLP, 1 for unmet load
// epsilon: number in range [0,1] representing LOLP or unmet load fraction.
// chunk_size: length of time (in days)
// this function will be called at the bottom of the file
SimulationResult run_simulations(vector<double> &load, vector<double> &solar, vector<double> &ev,  int metric, int chunk_size, int number_of_chunks)
{

	// set random seed to a specific value if you want consistency in results
	srand(10);

	// get number of hours in each chunk of 100 houes
	int t_chunk_size = chunk_size*(24/T_u);

	// vector of vectors : vector of all sizing results 
	vector <vector<SimulationResult> > results;

	// get random start times and run simulation on this chunk of data
		//for (int chunk_num = 0; chunk_num < 100; chunk_num += 1){
			for (int chunk_num = 0; chunk_num < 1; chunk_num += 1)
			{

				cout << "-------------chunk number = " << chunk_num << endl;

				// random hour of random day
				int chunk_start_f = rand() % max(solar.size(), load.size());

				// get the weekday
				int one_week = 168;
				int hour_day_sampled = chunk_start_f % one_week;
				int chunk_start = chunk_start_f - hour_day_sampled + one_week;
				int chunk_end = chunk_start + t_chunk_size;

				vector<SimulationResult> sr = simulate(load, solar, ev, chunk_start, chunk_end, 0);

				// DELETE ME: just for debugging
				// for (int i = 0; i < sr.size(); i++)
				for (int i = 50; i < 54; i++)
				{
					cout << "-------------chunk number = " << chunk_num << endl;

					SimulationResult s = sr.at(i);
					cout << "simulation result B " << s.B << endl;
					cout << "simulation result C " << s.C << endl;
					cout << "simulation result Cost " << s.cost << endl;
			}

			// each sr is a sizing curve that I want to print in the graph
			results.push_back(sr);
			}

#ifdef DEBUG
	// print all of the curves
	int chunk_index = 1;
	cout << "DEBUG: sizing_curves" << endl;
	for (vector<vector<SimulationResult>>::iterator it = results.begin(); it != results.end(); ++it, ++chunk_index) {
		cout << "chunk_" << chunk_index << endl;
		for (vector<SimulationResult>::iterator it2 = it->begin() ; it2 != it->end(); ++it2) {
			cout << it2->B << "\t" << it2->C << "\t" << it2->cost << endl;
		}
	}
	cout << "DEBUG: sizing_curves_end" << endl;
#endif

	// calculate the chebyshev curves, find the cheapest system along their upper envelope, and return it
	return calculate_sample_bound(results, epsilon, confidence);
}

int main(int argc, char ** argv) {
	
	int input_process_status = process_input(argv, true);

	// if process_input returns a 1, there was an error in processing the input.
	if (input_process_status) {
		cerr << "Illegal input" << endl;
		return 1;
	}
	
	SimulationResult sr = run_simulations(load, solar, ev, metric, days_in_chunk, number_of_chunks);

	// a temporary fix for "inf" issues.
	// TODO: investigate later why some sr.cost == inf when B or C are low
	
	double cost = sr.B / kWh_in_one_cell * B_inv + sr.C * PV_inv;
	cout << "results battery, solar, cost : "<< sr.B << "\t" << sr.C << "\t" << cost << endl;

	return 0;
}
