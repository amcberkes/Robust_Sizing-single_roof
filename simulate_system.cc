#include <cmath>
#include <iostream>
#include <vector>
#include "simulate_system.h"
#include <random>
#include <iostream>
#include <cstdlib>
#include <ctime> // For time()

using namespace std;
double current_ev_kWh;
double static ev_goal_kWh = 0.8 * 60.0;
double static charging_rate = 7.4;
int loss_events = 0;
double load_deficit = 0;
double load_sum = 0;
double ev_b = 0.0;
double b = 0.0;
double static t_ch = 3;

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

/*
decrease the applied (charging) power by increments of (1/30) until the power is
low enough to avoid violating the upper energy limit constraint.

returns the max amount of power that we can charge the stationary battery with
*/
double calc_max_charging(double power, double b_prev) {

	double step = power/30.0;

	for (double c = power; c >= 0; c -= step) {
		double upper_lim = a2_slope * (c / nominal_voltage_c) + a2_intercept ;
		double b = b_prev + c*eta_c*T_u;
		if (b <= upper_lim) {
			return c;
		}
	}
	return 0;
}

/*
decrease the applied (discharging) power by increments of (1/30) until the power is
 low enough to avoid violating the lower energy limit constraint.

 returns the max amount of power that we can discharge the stationary battery with

*/
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
double calc_max_charging_ev(double power, double ev_b, double ev_goal_kWh)
{
	double step = power / 30.0;
	for (double c = power; c >= 0; c -= step){
		// 80% SOC is lower limit
		double upper_lim = ev_goal_kWh;
		//cout << "ev_b before update l.84: " << ev_b << endl;
		ev_b = ev_b + c * eta_c_ev * T_u;
		//cout << "ev_b after update: " << ev_b << endl;

		if (ev_b <= upper_lim){
			return c;
		}
	}
	return 0;
}

//max discharging for ev battery
double calc_max_discharging_ev(double power, double ev_b){
	double step = power / 30.0;
	for (double d = power; d >= 0; d -= step){
		// 20% SOC is lower limit
		double lower_lim = 0.2 * 60.0;
		//cout << "ev_b before update l.101: " << ev_b << endl;

		ev_b = ev_b - d * eta_d_ev * T_u;
		//cout << "ev_b after update: " << ev_b << endl;

		if (ev_b >= lower_lim){
			return d;
		}
	}
	return 0;
}


// EV charging control 
int naive(int t, double ev_b, int next_dept){
	return t;
}

int lastp(int t, double ev_b, int next_dept){
	// compute pink line: ensure that ev is 80% charged when it leaves the house
	double diff = ev_goal_kWh - ev_b;
	int hours_needed = ceil(diff / charging_rate);
	int latest_t = next_dept - hours_needed;
	return latest_t;
}

int mincost(int t, double ev_b, int next_dept){
	return t_ch;
}


//Real Time Management
void unidirectional(bool z, double ev_b, double c, double d, double max_c, double max_d, double max_c_ev, double max_d_ev, double b){
	//charge: 1=stationary, 2= ev, 3 = verloren
	// discharge: 1= stationary, 2 = grid
	if(c > 0){
		if (c <= max_c){
			b = b + c * eta_c * T_u;
		} else {
			double res = c - max_c;
			b = b + max_c * eta_c * T_u;
			if(res <= max_c_ev){
				ev_b = ev_b + res * eta_c_ev * T_u;
			}else{
				ev_b = ev_b + max_c_ev * eta_c_ev * T_u;
				}
		}
	} 
	if(d > 0){
		if(d <= max_d){
			b = b - d * eta_d * T_u;
		} else {
			double res = d - max_d;
			b = b - max_d * eta_d * T_u;
			loss_events += 1;
			load_deficit += res;
		}

	}
}

void r_degradation(bool z,double ev_b,double c,double d,double max_c,double max_d,double max_c_ev,double max_d_ev, double b){
	// charge: 1=stationary, 2= ev, 3 = verloren
	//  discharge: 1= stationary, 2 = ev, 3= grid
	if (c > 0){
		if (c <= max_c){
			b = b + c * eta_c * T_u;
		}
		else{
			double res = c - max_c;
			b = b + max_c * eta_c * T_u;
			if (res <= max_c_ev){
				ev_b = ev_b + res * eta_c_ev * T_u;
			}else{
				ev_b = ev_b + max_c_ev * eta_c_ev * T_u;
			}
		}
	}
	if (d > 0){
		if (d <= max_d){
			b = b - d * eta_d * T_u;
		} else{
			double res = d - max_d;
			b = b - max_d * eta_d * T_u;
			if (res <= max_d_ev){
				ev_b = ev_b - res * eta_d_ev * T_u;
			} else {
				ev_b = ev_b - max_d_ev * eta_c_ev * T_u;
				loss_events += 1;
				double back = res - max_d_ev;
				load_deficit += back;
			}

			
		}
	}
}

void min_storage(bool z,double ev_b,double c,double d,double max_c,double max_d,doube max_c_ev,double max_d_ev, double b){
	// charge: 1=ev, 2= stationary, 3 = verloren
	//  discharge: 1= ev, 2 = stationary, 3= grid
	if (c > 0){
		if (c <= max_c_ev){
			ev_b = ev_b + c * eta_c_ev * T_u;
		} else{
			double res = c - max_c_ev;
			ev_b = ev_b + max_c_ev * eta_c_ev * T_u;
			if (res <= max_c){
				b = b + res * eta_c * T_u;
			}else{
				b = b + max_c * eta_c * T_u;
			}
		}
	}
	if (d > 0){
		if (d <= max_d_ev){
			ev_b = ev_b - d * eta_d_ev * T_u;
		} else{
			double res = d - max_d_ev;
			ev_b = ev_b - max_d_ev * eta_d_ev * T_u;
			if (res <= max_d){
				b = b - res * eta_d * T_u;
			} else{
				b = b - max_d * eta_c * T_u;
				loss_events += 1;
				double back = res - max_d;
				load_deficit += back;
			}
		}
	}
}

void most_sustainable(bool z, double ev_b, double c, double d, double max_c, double max_d, double max_c_ev, double max_d_ev, double b){
	// charge: 1=ev, 2= stationary, 3 = verloren
	//  discharge: 1= stationary, 2 = ev, 3= grid
	if (c > 0){
		if (c <= max_c_ev){
			ev_b = ev_b + c * eta_c_ev * T_u;
		} else{
			double res = c - max_c_ev;
			ev_b = ev_b + max_c_ev * eta_c_ev * T_u;
			if (res <= max_c){
				b = b + res * eta_c * T_u;
			} else{
				b = b + max_c * eta_c * T_u;
			}
		}
	}
	if (d > 0){
		if (d <= max_d){
			b = b - d * eta_d * T_u;
		} else{
			double res = d - max_d;
			b = b - max_d * eta_d * T_u;
			if (res <= max_d_ev){
				ev_b = ev_b - res * eta_d_ev * T_u;
			} else{
				ev_b = ev_b - max_d_ev * eta_c_ev * T_u;
				loss_events += 1;
				double back = res - max_d_ev;
				load_deficit += back;
			}
		}
	}
}

// Note: sim_year calls procedures calc_max_charging and calc_max_discharging.
// You could potentially speed up the computation by expanding these functions into sim_year
// to avoid procedure calls in this inner loop.
double sim(vector<double> &load_trace, vector<double> &solar_trace, vector<double> &ev_trace, int start_index, int end_index, double cells, double pv, double b_0){

	update_parameters(cells);

	// set the battery
	 b = b_0*cells*kWh_in_one_cell; //0.5*a2_intercept
	
	int num_trips = 0;
	int t_dept;
	int t_arr;
	double soc_arr;
	bool ev_at_home[24];
	double ev_soc[24];

	int trace_length_solar = solar_trace.size();
	int trace_length_load = load_trace.size();
	int ev_trace_index = 0;

	int next_dept;
	double c = 0.0;
	double d = 0.0;
	double max_c = 0.0;
	double max_d = 0.0;
	double max_c_ev = 0.0;
	double max_d_ev = 0.0;
	int index_t_solar;
	int index_t_load;
	int counter = 0;
	int t_charge = 0;
	bool z = false;

	loss_events = 0;
	load_deficit = 0;
	load_sum = 0;

	//for each of the (100) days in the sample
	//for (int i = start_index; i < start_index + days_in_chunk; i++){
	for (int i = start_index ; i < start_index + 5; i++){
		cout << " DAY NUMBER  : " << i % start_index << endl;
		ev_trace_index = i % start_index + counter;

		for (int k = 0; k < 24; k++){
			ev_at_home[k] = true;
		}

		// ----------------------------------------------------------------read EV inputs-----------------------------------
		num_trips = ev_trace[ev_trace_index];
		cout << "num_trips : " << num_trips << endl;
		// next departures
		int array_length = num_trips + 1;
		int next_dept_arr[array_length];
		int next_dept_size = sizeof(next_dept_arr) / sizeof(int);
		// cout << "1 - length of next_dept_arr  is : " << next_dept_size << endl;

		for (int j = 0; j < num_trips; j++){
			ev_trace_index = ev_trace_index + 1;
			counter = counter + 1;
			t_dept = ev_trace[ev_trace_index];
			cout << "t_dept : " << t_dept << endl;

			ev_trace_index = ev_trace_index + 1;
			counter = counter + 1;
			t_arr = ev_trace[ev_trace_index];
			cout << "t_arr : " << t_arr << endl;

			ev_trace_index = ev_trace_index + 1;
			counter = counter + 1;
			soc_arr = ev_trace[ev_trace_index];
			cout << "soc_Arr : " << soc_arr << endl;

			// next departure
			next_dept_arr[j] = t_dept;

			for (int h = t_dept; h < t_arr; h++){
				ev_at_home[h] = false;
			}
			ev_soc[t_arr] = soc_arr * 60;
		}

		// if next num_dept is 0: counts number of consecutive days with no trips
		int pad = 0;
		//  get num trips from the next day
		int ev_trace_index_2 = i % start_index + counter + 1;
		while (ev_trace[ev_trace_index_2 + pad] == 0){
			//  on the next day there are no trips
			pad = pad + 1;
		}
		next_dept_arr[num_trips] = ev_trace[ev_trace_index + 2 + pad];
		pad = 0;

		//----------------------------------------------------------------------------------- Start hourly EMS --------------

		for(int t = 0; t<24; t++){
			int index = t + 24*i;
			index_t_solar = index % trace_length_solar;
			index_t_load = index % trace_length_load;
			load_sum += load_trace[index_t_load];
			double load = load_trace[index_t_load];

			int trips_counter = 0;
			if (next_dept == t){
				trips_counter = trips_counter + 1;
			}
			next_dept = next_dept_arr[trips_counter];

			if (ev_soc[t] == 0 && t == 0 && i == 0){
				// only runs for first monday 0h : assume soc is 60% charged initially
				ev_b = 36;
			} else{
				ev_b = ev_soc[t];
			}

			//----------------------------------------------------   EV Charging Control --------------------------------
			if (ev_at_home[t]){
				t_charge = naive(t, ev_b, next_dept);
				// t_charge = lastp(t, ev_b, next_dept);
				// t_charge = mincost(t, ev_b, next_dept);
				if(t == t_charge){
					z = true;
				}
				if (z){
					c = fmax(solar_trace[index_t_solar] * pv - load - 7.4, 0);
					d = fmax(load + 7.4 - solar_trace[index_t_solar] * pv, 0);
				} else {
					c = fmax(solar_trace[index_t_solar] * pv - load, 0);
					d = fmax(load - solar_trace[index_t_solar] * pv, 0);
				}
				max_c = fmin(calc_max_charging(c, b), alpha_c);
				max_d = fmin(calc_max_discharging(d, b), alpha_d);
				max_c_ev = fmin(calc_max_charging_ev(c, ev_b, ev_goal_kWh), alpha_c_ev);
				max_d_ev = fmin(calc_max_discharging_ev(d, ev_b), alpha_d_ev);
			}
			else{
				z = false;
				ev_b = 0.0;
				c = fmax(solar_trace[index_t_solar] * pv - load, 0);
				d = fmax(load - solar_trace[index_t_solar] * pv, 0);
				max_c = fmin(calc_max_charging(c, b), alpha_c);
				max_d = fmin(calc_max_discharging(d, b), alpha_d);
				max_c_ev = 0.0;
				max_d_ev = 0.0;
			}

			//----------------------------------------------------   Real time Management ----------------------------------------
			unidirectional(z, ev_b, c, d, max_c, max_d, max_c_ev, max_c_ev, b);
			//r_degradation(z, ev_b, c, d, max_c, max_d, max_c_ev, max_c_ev, b);
			//min_storage(z, ev_b, c, d, max_c, max_d, max_c_ev, max_c_ev, b);
			//most_sustainable(z, ev_b, c, d, max_c, max_d, max_c_ev, max_c_ev, b);

			if(z){
				ev_b = ev_b + 7.4;
			}
			ev_soc[(t + 1) % 24] = ev_b;
		}
	}

	if (metric == 0) {
		// lolp
		double result = loss_events / ((100) * 1.0) ;
		//cout << "RESULT NUMBER OF LOSS EVENTS" << result << endl;
		return result;
	} else {
		// metric == 1, eue
		//cout << "RESULT load deficit = " << load_deficit << endl;
		double result = load_deficit / (load_sum * 1.0);
		//cout << "RESULT LOSS" << result <<endl;
		return result;
	}
}

// Run simulation for provides solar and load trace to find cheapest combination of
// load and solar that can meet the epsilon target

// this function is called for each sample: constructs one sizing curve
vector<SimulationResult> simulate(vector<double> &load_trace, vector<double> &solar_trace, vector<double> &ev_trace, int start_index, int end_index, double b_0)
{

	// for Cmax, find Bmin such that system still satisfies the target performance epsilon
	double cells_U = cells_max;
	double cells_L = cells_min;
	double mid_cells = 0.0;
	double loss = 0.0;

	while (cells_U - cells_L > cells_step) {

		mid_cells = (cells_L + cells_U) / 2.0;
		loss = sim(load_trace, solar_trace, ev_trace, start_index, end_index, mid_cells, pv_max, b_0);

		// binary search
		if (loss > epsilon) {
			cells_L = mid_cells;
		} else {
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
	// first point of the sizing curve at Cmax
	curve.push_back(SimulationResult(starting_cells*kWh_in_one_cell, lowest_feasible_pv, starting_cost));

	for (double cells = starting_cells; cells <= cells_max; cells += cells_step) {
	//for (double cells = starting_cells; cells <= starting_cells + cells_step; cells += cells_step){

		// for each value of cells, find the lowest pv that meets the epsilon loss constraint
		double loss = 0;
		while (true)
		{
			//cout << "-------call sim() with the following number of cells = "<< cells<<endl;
			// compute loss of current sizing. 
			loss = sim(load_trace, solar_trace, ev_trace, start_index, end_index, cells, lowest_feasible_pv - pv_step, b_0);
			//remove this
			break;
			cout << "completed simulation with loss =  " << loss << endl;

			if (loss < epsilon) {
				// we can try an ev en smaller pv size, since epsilon not violated yet
				lowest_feasible_pv -= pv_step;
			} else {
				// break exits the innermost loop containing it: breaks out of while loop
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

// why do we need this?
		if (cost < lowest_cost) {
			lowest_cost = cost;
			lowest_B = cells*kWh_in_one_cell;
			lowest_C = lowest_feasible_pv;
		}
		}

	return curve;
}
