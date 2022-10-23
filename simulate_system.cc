#include <cmath>
#include <iostream>
#include <vector>
#include "simulate_system.h"

using namespace std;

// parameters specified for an NMC cell with operating range of 1 C charging and discharging

void update_parameters(double n) {

	num_cells = n;

	a1_intercept = 0.0*num_cells;
	
	a2_intercept = kWh_in_one_cell*num_cells;
	
	alpha_d = a2_intercept*1.0;
	alpha_c = a2_intercept*1.0;
	return;
}

// decrease the applied (charging) power by increments of (1/30) until the power is 
// low enough to avoid violating the upper energy limit constraint.
double calc_max_charging(double power, double b_prev, bool ev) {

	double step = power/30.0;

	for (double c = power; c >= 0; c -= step) {
		double upper_lim = 0.0;
		
		if(ev){
			//theoretisches upper limit
			 upper_lim = a2_slope * (c / nominal_voltage_c) + a2_intercept + 18;
		}else{
			 upper_lim = a2_slope * (c / nominal_voltage_c) + a2_intercept ;
		}
		

		//upper_lim = a2_slope * (c / nominal_voltage_c) + a2_intercept;
		//cout << "upper_lim = " << upper_lim << endl;
		
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
	// ev_b says how much the battery in ev is charged
	double ev_b = 0.0;
	for (int t = start_index; t < end_index; t++) {
		//cout << "current index" << t << endl;
		// wrap around to the start of the trace if we hit the end.
		index_t_solar = t % trace_length_solar;
		index_t_load = t % trace_length_load;

		load_sum += load_trace[index_t_load];

		// first, calculate how much power is available for charging, and how much is needed to discharge
		//c = how much power is available to charge battery
		if (t % 24 == 23 || t % 24 == 0 || t % 24 == 1 || t % 24 == 2 || t % 24 == 3 || t % 24 == 4 || t % 24 == 5 || t % 24 == 6)
		{
			c = fmax(solar_trace[index_t_solar] * pv - load_trace[index_t_load] - 1, 0);
			// d = how much energy we need to take out of battery
			d = fmax(load_trace[index_t_load] + 1 - solar_trace[index_t_solar] * pv, 0);
		} else {
			c = fmax(solar_trace[index_t_solar] * pv - load_trace[index_t_load], 0);
			// d = how much energy we need to take out of battery
			d = fmax(load_trace[index_t_load] - solar_trace[index_t_solar] * pv, 0);
		}
		//How does ev_b vary with the time ? 
		int tt = t %24;
		
		
		//cout << "time = " << tt << endl;
		switch(tt){
			case 0:
				ev_b = 12;
				ev = true;
				break;
			case 1:
				ev_b = 13;
				ev = true;
				break;
			case 2:
				ev_b = 14;
				ev = true;
				break;
			case 3:
				ev_b = 15;
				ev = true;
				break;
			case 4:
				ev_b = 16;
				ev = true;
				break;
			case 5:
				ev_b = 17;
				ev = true;
				break;
			case 6:
				ev_b = 18;
				ev = true;
				break;
			case 7:
				ev_b = 18;
				ev = true;
				break;
			case 8:
				ev_b = 18;
				ev = true;
				break;
			case 9:
				ev_b = 0;
				ev = false;
				break;
			case 10:
				ev_b = 0;
				ev = false;
				break;
			case 11:
				ev_b = 0;
				ev = false;
				break;
			case 12:
				ev_b = 0;
				ev = false;
				break;
			case 13:
				ev_b = 0;
				ev = false;
				break;
			case 14:
				ev_b = 0;
				ev = false;
				break;
			case 15:
				ev_b = 0;
				ev = false;
				break;
			case 16:
				ev_b = 0;
				ev = false;
				break;
			case 17:
				ev_b = 0;
				ev = false;
				break;
			case 18:
				ev_b = 10;
				ev = true;
				break;
			case 19:
				ev_b = 10;
				ev = true;
				break;
			case 20:
				ev_b = 10;
				ev = true;
				break;
			case 21:
				ev_b = 10;
				ev = true;
				break;
			case 22:
				ev_b = 10;
				ev = true;
				break;
			case 23:
				ev_b = 11;
				ev = true;
				break;
		}



		// constrain the power
		// alpha_c = kWh_in_one_cell*num_cells = max solar power generated
		// max_c is the max amount that we can charge b
		//!!!! b here is the current total battery charged (sum of charged in b and ev_b)
		//cout << "ev_b = " << ev_b << endl;
		//cout << "b before update = " << b << endl;
		
		// DIESE LINE IST FALSCH 
		double b_new = b + ev_b;
		
		//cout << "b after update = " << b << endl;
		max_c = fmin(calc_max_charging(c,b_new, ev), alpha_c);
		//cout << "c = " << c << endl;
		//cout << "max_c = " << max_c << endl;
		// alpha_d = alpha_cs
		// max_d is max amount that we can decharge
		
		max_d = fmin(calc_max_discharging(d,b), alpha_d);
		//cout << "d = " << d << endl;
		//cout << "max_d = " << max_d << endl;

		// at each time step either c or d is 0, which is why either max_c or max_d is 0 and it works out
		b = b + max_c*eta_c*T_u - max_d*eta_d*T_u ;
		//b = b - ev_b;
		// if we didnt get to discharge as much as we wanted, there is a loss
		if (max_d < d) {
			loss_events += 1;
			load_deficit += (d - max_d);
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
