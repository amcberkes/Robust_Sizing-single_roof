#include <cmath>
#include <iostream>
#include <vector>
#include "simulate_system.h"

using namespace std;

// parameters specified for an NMC cell with operating range of 1 C charging and discharging

void update_parameters(double n) {

	num_cells = n;

//lower energy content limit
	a1_intercept = 0.0*num_cells;
//upper energy content limit
	a2_intercept = kWh_in_one_cell*num_cells;
//max discharging rate	
	alpha_d = a2_intercept*1.0;
// max charging rate
	alpha_c = a2_intercept*1.0;
	return;
}

//update params for ev
void update_parameters_ev(double n){
	num_cells_ev = n;
	// lower energy content limit
	a1_intercept_ev = 0.0 * num_cells_ev;
	// upper energy content limit
	a2_intercept_ev = kWh_in_one_cell_ev * num_cells_ev;
	// max discharging rate
	alpha_d_ev = a2_intercept_ev * 1.0;
	// max charging rate
	alpha_c_ev = a2_intercept_ev * 1.0;
	return;
}

// decrease the applied (charging) power by increments of (1/30) until the power is 
// low enough to avoid violating the upper energy limit constraint.
double calc_max_charging(double power, double b_prev) {

	double step = power/30.0;

	for (double c = power; c >= 0; c -= step) {
		double upper_lim = a2_slope * (c / nominal_voltage_c) + a2_intercept ;
	
		//cout << "upper_lim = " << upper_lim << endl;
		//eta_c is the charging penalty
		double b = b_prev + c*eta_c*T_u;
		//cout << "b = " << b << endl;
		if (b <= upper_lim) {
			//cout << "GOOD : upper_lim > b " << endl;
			return c;
		}else {
			//cout << "BAD upper_lim < b " << endl;
		}
	}
	return 0;
}


// decrease the applied (discharging) power by increments of (1/30) until the power is
// low enough to avoid violating the lower energy limit constraint.
double calc_max_discharging(double power, double b_prev) {

	double step = power/30.0;

	for (double d = power; d >= 0; d -= step) {
		double lower_lim = a1_slope*(d/nominal_voltage_d) + a1_intercept;
		double b = b_prev - d*eta_d*T_u;
		if (b >= lower_lim) {
			return d;
		}
	}
	return 0;
}


//max charging function for ev battery
double calc_max_charging_ev(double power, double b_prev, bool ev){
	double step = power / 30.0;
	for (double c = power; c >= 0; c -= step){
		double upper_lim = a2_slope_ev * (c / nominal_voltage_c_ev) + a2_intercept_ev;
		double b = b_prev + c * eta_c_ev * T_u;
		if (b <= upper_lim){
			return c;
		}
		else
		{
			// cout << "BAD upper_lim < b " << endl;
		}
	}
	return 0;
}

//max discharging for ev battery
double calc_max_discharging_ev(double power, double b_prev){
	double step = power / 30.0;
	for (double d = power; d >= 0; d -= step){
		double lower_lim = a1_slope_ev * (d / nominal_voltage_d_ev) + a1_intercept_ev;
		double b = b_prev - d * eta_d_ev * T_u;
		if (b >= lower_lim){
			return d;
		}
	}
	return 0;
}

// Note: sim_year calls procedures calc_max_charging and calc_max_discharging.
// You could potentially speed up the computation by expanding these functions into sim_year
// to avoid procedure calls in this inner loop.
double sim(vector <double> &load_trace, vector <double> &solar_trace, int start_index, int end_index, double cells, double pv, double b_0) {

	update_parameters(cells);

	// set the battery
	double b = b_0*cells*kWh_in_one_cell; //0.5*a2_intercept

	int loss_events = 0;

	double load_deficit = 0;
	double load_sum = 0;

	int trace_length_solar = solar_trace.size();
	int trace_length_load = load_trace.size();

	double c = 0.0;
	double d = 0.0;
	double max_c = 0.0;
	double max_d = 0.0;
	int index_t_solar;
	int index_t_load;
	
	bool ev = false;
	double ev_b = 19.5;
	int arrival_time = 18;
	int departure_time = 9;
	int time_length = arrival_time - departure_time;
	bool arr_time[time_length] ;
	double ev_b_state = 0.5;

/*
	for(int i = 0; i < 24; i++){
		if (i == 18 || i == 19|| i == 20|| i == 21|| i == 22|| i == 23|| i == 0|| i == 1|| i == 2|| i == 3|| i == 4|| i == 5|| i == 6|| i == 7|| i == 8){
			arr_time[i] = true;
			cout << "set to true" << endl;
		} else {
			arr_time[i] = false;
			cout << "set to false" << endl;
		}
		
	}*/
	
	for (int t = start_index; t < end_index; t++) {
		//cout << "current index" << t << endl;
		// wrap around to the start of the trace if we hit the end.
		index_t_solar = t % trace_length_solar;
		index_t_load = t % trace_length_load;

		load_sum += load_trace[index_t_load];
		int tt = t % 24;
		cout << "time: " << tt << endl;
		if(arr_time[tt] == true){
			ev = true;
			//cout << "time = " << tt << endl;
			//cout << "ev_-b = " << ev_b << endl;
		}
		
		if(tt == arrival_time){
			ev_b = 19.5;
		}

		if (tt == 23 || tt == 0 || tt == 1 || tt == 2 || tt == 3 || tt == 4 || tt == 5 || tt == 6){
			c = fmax(solar_trace[index_t_solar] * pv - load_trace[index_t_load] - 1, 0);
			d = fmax(load_trace[index_t_load] + 1 - solar_trace[index_t_solar] * pv, 0);
			cout << "ev_b before increase: " << ev_b << endl;
			// charge the ev battery by 1 kwh
			ev_b = ev_b + 1;
			cout << "ev_b after increase: " << ev_b << endl;
		} else {
			c = fmax(solar_trace[index_t_solar] * pv - load_trace[index_t_load], 0);
			d = fmax(load_trace[index_t_load] - solar_trace[index_t_solar] * pv, 0);
		}
		cout << "c: " << c << endl;
		cout << "d: " << d << endl;
		cout << "ev_b: " << ev_b << endl;
		cout << "b: " << b << endl;

		// alpha_c = kWh_in_one_cell*num_cells = max solar power generated
		max_c = fmin(calc_max_charging(c,b), alpha_c);
		max_d = fmin(calc_max_discharging(d,b), alpha_d);

		cout << "max_c: " << max_c << endl;
		cout << "max_d: " << max_d << endl;

		//differenet charging policies here
		bool stat_prioritized = true;
		if (tt == 18 || tt == 19 || tt == 20 || tt == 21 || tt == 22 || tt == 23 || tt == 0 || tt == 1 || tt == 2 || tt == 3 || tt == 4 || tt == 5 || tt == 6 || tt == 7 || tt == 8)
		{
			cout << "ev = true with time: " << tt << endl;
			double max_c_ev = fmin(calc_max_charging_ev(c, ev_b, ev), alpha_c_ev);
			double max_d_ev = fmin(calc_max_discharging_ev(d, ev_b), alpha_d_ev);
			if(stat_prioritized){
				if(c > max_c){
					
					double rest_c = c - max_c;
					cout << "charge ev with: " << rest_c << endl;
					//we generated more charge than we can store in stationary battery and will store the excess charge in the ev
					if(rest_c < max_c_ev){
						ev_b = ev_b + rest_c * eta_c_ev * T_u;
						b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
					}else{
						ev_b = ev_b + max_c * eta_c_ev * T_u;
						b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
					}
				} else if(max_d < d){
					// we need to discharge the ev battery
					double rest_d = d - max_d;
					cout << "discharge ev with: " << rest_d << endl;
					cout << "ev_b before discharge with: " << ev_b << endl;
					cout << "b before evb discharge with: " << b << endl;
					if (rest_d > max_d_ev){
						loss_events += 1;
						load_deficit += (rest_d - max_d_ev);
						
						ev_b = ev_b - max_d_ev * eta_d_ev * T_u;
						b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
					}else{
						cout << "rest_d_ev before discharge: " << rest_d << endl;
						ev_b = ev_b - rest_d * eta_d_ev * T_u;
						b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
					}
				}
				else {
					cout << "ev not needed, b before =  " << b << endl;
					b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
					cout << "ev not needed, b after =  " << b << endl;
				}
			}
		}
		else
		{
			cout << "ev = false with time: " << tt << endl;

			// at each time step either c or d is 0, which is why either max_c or max_d is 0 and it works out
			cout << "b before update: " << b << endl;
			b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
			cout << "b after update: " << b << endl;
			//  if we didnt get to discharge as much as we wanted, there is a loss
			if (max_d < d){
				loss_events += 1;
				load_deficit += (d - max_d);
			}
		}
	}

	if (metric == 0) {
		// lolp
		return loss_events/((end_index - start_index)*1.0);
	} else {
		// metric == 1, eue
		return load_deficit/(load_sum*1.0);
	}
}


// Run simulation for provides solar and load trace to find cheapest combination of
// load and solar that can meet the epsilon target
vector <SimulationResult> simulate(vector <double> &load_trace, vector <double> &solar_trace, int start_index, int end_index, double b_0) {

	// first, find the lowest value of cells that will get us epsilon loss when the PV is maximized
	// use binary search
	double cells_U = cells_max;
	double cells_L = cells_min;
	double mid_cells = 0.0;
	double loss = 0.0;

	while (cells_U - cells_L > cells_step) {

		mid_cells = (cells_L + cells_U) / 2.0;

		loss = sim(load_trace, solar_trace, start_index, end_index, mid_cells, pv_max, b_0);

		//cout << "sim result with " << a2_intercept << " kWh and " << pv_max << " pv: " << loss << endl;
		if (loss > epsilon) {
			cells_L = mid_cells;
		} else {
		 	// (loss <= epsilon)
			cells_U = mid_cells;
		}
	}

	// set the starting number of battery cells to be the upper limit that was converged on
	double starting_cells = cells_U;
	double starting_cost = B_inv*starting_cells + PV_inv * pv_max;
	double lowest_feasible_pv = pv_max;


	double lowest_cost = starting_cost;
	double lowest_B = starting_cells*kWh_in_one_cell;
	double lowest_C = pv_max;

	vector <SimulationResult> curve;
	curve.push_back(SimulationResult(starting_cells*kWh_in_one_cell, lowest_feasible_pv, starting_cost));
	//cout << "starting cells: " << starting_cells << endl;

	for (double cells = starting_cells; cells <= cells_max; cells += cells_step) {

		// for each value of cells, find the lowest pv that meets the epsilon loss constraint
		double loss = 0;
		while (true) {
			
			loss = sim(load_trace, solar_trace, start_index, end_index, cells, lowest_feasible_pv - pv_step, b_0);

			if (loss < epsilon) {
				lowest_feasible_pv -= pv_step;
			} else {
				break;
			}

			// this only happens if the trace is very short, since the battery starts half full
			// and can prevent loss without pv for a short time
			if (lowest_feasible_pv <= 0) {
				lowest_feasible_pv = 0;
				break;
			}
		}

		double cost = B_inv*cells + PV_inv*lowest_feasible_pv;

		curve.push_back(SimulationResult(cells*kWh_in_one_cell,lowest_feasible_pv, cost));

		if (cost < lowest_cost) {
			lowest_cost = cost;
			lowest_B = cells*kWh_in_one_cell;
			lowest_C = lowest_feasible_pv;
		}

	} 

	//return SimulationResult(lowest_B, lowest_C, lowest_cost);
	return curve;
}




// now I have modeled that ev battery can be used as backup when stationary battery is full or empty



/*still need to model that ev is fully charged when it leaves : muss load part verändern 
	- 7h für eine volle ladung 
	- 39 kwh
	- kommt halb voll an --> 3.5h ladung zwischen von 19,5kwh 

*/


// and need to model stochastic arrival and departure time 

// then add other operating policies 

// then add flexible initial ev battery state and load schedule 