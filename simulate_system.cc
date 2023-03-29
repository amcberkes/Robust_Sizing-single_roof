#include <cmath>
#include <iostream>
#include <vector>
#include "simulate_system.h"
#include <random>
#include <iostream>
#include <cstdlib>
#include <ctime> // For time()

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
void update_parameters_ev(){
	num_cells_ev = 192;
	// lower energy content limit
	a1_intercept_ev = 0.0 * num_cells_ev;
	// upper energy content limit
	a2_intercept_ev = kWh_in_one_cell_ev * num_cells_ev;
	//a2_intercept_ev = 60;
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
		double b = b_prev + c*eta_c*T_u;
		if (b <= upper_lim) {
			return c;
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
	update_parameters_ev();

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
	double ev_b = 62.0;
	int arrival_time = 18;
	int departure_time = 9;

	// generate stochastic arrival time 
		default_random_engine e(18);
		std::normal_distribution<double> distribution(18.0, 2.0);
		srand(time(0)); // Initialize random number generator.
		int index = (rand() % 10) + 1;
		for (int i = 0; i < 10; i++)
		{
			int arr_time = distribution(e);
			if (i == index)
			{
				arrival_time = arr_time;
				//  cout << "arrival_time: " << arrival_time << endl;
			}
		}
	
	
	
	// generate stochastic departure time

		std::normal_distribution<double> distribution2(9, 1.0);
		srand(time(0)); // Initialize random number generator.
		int index2 = (rand() % 10) + 1;
		for (int i = 0; i < 10; i++)
		{
			int dept_time = distribution2(e);
			if (i == index2)
			{
				departure_time = dept_time;
				// cout << "departure_time: " << departure_time << endl;
			}
		}
	
	

	for (int t = start_index; t < end_index; t++) {
		// wrap around to the start of the trace if we hit the end.
		index_t_solar = t % trace_length_solar;
		index_t_load = t % trace_length_load;

		load_sum += load_trace[index_t_load];
		int tt = t % 24;
		cout << "HALLOOOO time: " << tt << endl;

		if (tt == arrival_time)
		{
			// assume that ev arrives 70% charged
			ev_b = 43.4;
		}
		if (tt == departure_time)
		{
			// cout << "ev_b before departure: " << ev_b << endl;
			ev_b = 0;
		}
		if (tt >= arrival_time || tt < 3)
		{
			ev = true;
		}
		else{
			ev = false;
		}
		bool new_person = false;
		bool heat_pump = false;
		bool second_ev = false;
		double coef = 1.0;
		double coef_ev = 1.0;
		double hp_load = 0.0;


		double load = load_trace[index_t_load] * coef;
		double ev_load = coef_ev*6.6;

		// different charging policies here
		bool stat_prioritized = false;
		bool ev_prioritized = false;
		bool round_robin = false;
		// store charge in ev when it is there and discharge battery
		bool c_ev_d_b = false;
		// uni-directional ev, just addiotional load
		bool c_ev_d_b_only = false;



		// pink line: ensure that ev is 80% charged when it leaves the house
		// make sure that it is the time in hours 
		double pink_line = 6.6*tt;
		cout << "pink linr time: " << pink_line << endl;

		if(ev_b > pink_line){
			// we do not need to charge the ev yet. normal bi-directional behaviour
			ev_prioritized = true;
			c = fmax(solar_trace[index_t_solar] * pv - load, 0);
			d = fmax(load - solar_trace[index_t_solar] * pv, 0);
			
		}else{
			// need to charge the ev now. charge it at the max rate and do not discharge it anymore
			ev = false;
			c_ev_d_b_only = true;
			c = fmax(solar_trace[index_t_solar] * pv - load - ev_load, 0);
			d = fmax(load + ev_load - solar_trace[index_t_solar] * pv, 0);
			double max_c_ev = fmin(calc_max_charging_ev(6.6, ev_b, ev), alpha_c_ev);
			while(ev_b < 49.6){
				ev_b = ev_b + max_c_ev;
			}
		
		}


		//cout << "c: " << c << endl;
		//cout << "d: " << d << endl;
		//cout << "ev_b: " << ev_b << endl;
		//cout << "b: " << b << endl;

		// alpha_c = kWh_in_one_cell*num_cells = max solar power generated
		max_c = fmin(calc_max_charging(c,b), alpha_c);
		max_d = fmin(calc_max_discharging(d,b), alpha_d);

		//cout << "max_c: " << max_c << endl;
		//cout << "max_d: " << max_d << endl;

	
		if (ev)
		{
			// cout << "inside ev true blcok" << endl;
			// cout << "ev = true with time: " << tt << endl;
			double max_c_ev = fmin(calc_max_charging_ev(c, ev_b, ev), alpha_c_ev);
			double max_d_ev = fmin(calc_max_discharging_ev(d, ev_b), alpha_d_ev);
			if (round_robin)
			{
				// cout << "inside round robin block " << endl;
				stat_prioritized = false;
				ev_prioritized = false;
				if (tt % 2 == 0)
				{
					stat_prioritized = true;
				}
				else
				{
					ev_prioritized = true;
				}
			}
			if (stat_prioritized)
			{
				if (c > max_c)
				{

					double rest_c = c - max_c;
					// cout << "charge ev with: " << rest_c << endl;
					// we generated more charge than we can store in stationary battery and will store the excess charge in the ev
					if (rest_c < max_c_ev)
					{
						ev_b = ev_b + rest_c * eta_c_ev * T_u;
						b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
					}
					else
					{
						ev_b = ev_b + max_c_ev * eta_c_ev * T_u;
						b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
					}
				}
				else if (max_d < d)
				{
					// we need to discharge the ev battery
					double rest_d = d - max_d;
					// cout << "discharge ev with: " << rest_d << endl;
					// cout << "ev_b before discharge with: " << ev_b << endl;
					// cout << "b before evb discharge with: " << b << endl;
					if (rest_d > max_d_ev)
					{
						loss_events += 1;
						load_deficit += (rest_d - max_d_ev);

						ev_b = ev_b - max_d_ev * eta_d_ev * T_u;
						b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
						// cout << "ev_b after not full discharge: " << ev_b << endl;
					}
					else
					{
						// cout << "rest_d_ev before discharge: " << rest_d << endl;
						ev_b = ev_b - rest_d * eta_d_ev * T_u;
						// cout << "ev_b after  full discharge: " << ev_b << endl;

						b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
					}
				}
				else
				{
					// cout << "ev not needed, b before =  " << b << endl;
					// cout << "only update stationary storage - before: " << b << endl;

					b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
					// cout << "only update stationary storage - after: " << b << endl;

					//  cout << "ev not needed, b after =  " << b << endl;
				}
			}
			if (ev_prioritized)
			{
				// cannot store charge in ev as it is already full
				if (c > max_c_ev)
				{

					double rest_c = c - max_c_ev;
					// cout << "charge ev with: " << rest_c << endl;
					// we generated more charge than we can store in stationary battery and will store the excess charge in the ev
					// can fully charge stationary with the rest
					if (rest_c < max_c)
					{
						ev_b = ev_b + max_c_ev * eta_c_ev * T_u;
						b = b + rest_c * eta_c * T_u - max_d * eta_d * T_u;
					}
					// some charge is lost
					else
					{
						ev_b = ev_b + max_c_ev * eta_c_ev * T_u;
						b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
					}
				}
				// cannot discharge ev, as it is empty
				else if (max_d_ev < d)
				{
					// we need to discharge the ev battery
					double rest_d = d - max_d_ev;
					// cout << "discharge ev with: " << rest_d << endl;
					// cout << "ev_b before discharge with: " << ev_b << endl;
					// cout << "b before evb discharge with: " << b << endl;
					// cannot fully discharge stationary eother
					if (rest_d > max_d)
					{
						loss_events += 1;
						load_deficit += (rest_d - max_d);

						ev_b = ev_b - max_d_ev * eta_d_ev * T_u;
						b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
						// cout << "ev_b after not full discharge: " << ev_b << endl;
					}
					else
					{
						// cout << "rest_d_ev before discharge: " << rest_d << endl;
						ev_b = ev_b - max_d_ev * eta_d_ev * T_u;
						// cout << "ev_b after  full discharge: " << ev_b << endl;

						b = b + max_c * eta_c * T_u - rest_d * eta_d * T_u;
					}
				}
				// only use ev
				else
				{
					// cout << "ev not needed, b before =  " << b << endl;
					// cout << "only update stationary storage - before: " << b << endl;
					ev_b = ev_b + max_c_ev * eta_c * T_u - max_d_ev * eta_d * T_u;
					// cout << "only update stationary storage - after: " << b << endl;

					// cout << "ev not needed, b after =  " << b << endl;
				}
			}
			if (c_ev_d_b_only)
			{
				// cannot store charge in ev as it is already full
				if (c > max_c_ev)
				{

					double rest_c = c - max_c_ev;
					// cout << "charge ev with: " << rest_c << endl;
					// we generated more charge than we can store in stationary battery and will store the excess charge in the ev
					// can fully charge stationary with the rest
					if (rest_c < max_c)
					{
						ev_b = ev_b + max_c_ev * eta_c_ev * T_u;
						b = b + rest_c * eta_c * T_u - max_d * eta_d * T_u;
					}
					// some charge is lost
					else
					{
						ev_b = ev_b + max_c_ev * eta_c_ev * T_u;
						b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
					}
				}
				// cannot discharge b, as it is empty
				else if (max_d < d)
				{
					// we need to discharge the ev battery
					double rest_d = d - max_d;
					// cout << "discharge ev with: " << rest_d << endl;
					// cout << "ev_b before discharge with: " << ev_b << endl;
					// cout << "b before evb discharge with: " << b << endl;

					loss_events += 1;
					load_deficit += (rest_d - max_d_ev);

					b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
				}
				// only use ev
				else
				{
					// cout << "ev not needed, b before =  " << b << endl;
					// cout << "only update stationary storage - before: " << b << endl;
					ev_b = ev_b + max_c_ev * eta_c * T_u;
					b = b - max_d * eta_d * T_u;
					// cout << "only update stationary storage - after: " << b << endl;

					// cout << "ev not needed, b after =  " << b << endl;
				}
			}
			if (c_ev_d_b)
			{
				// cannot store charge in ev as it is already full
				if (c > max_c_ev)
				{

					double rest_c = c - max_c_ev;
					// cout << "charge ev with: " << rest_c << endl;
					// we generated more charge than we can store in stationary battery and will store the excess charge in the ev
					// can fully charge stationary with the rest
					if (rest_c < max_c)
					{
						ev_b = ev_b + max_c_ev * eta_c_ev * T_u;
						b = b + rest_c * eta_c * T_u - max_d * eta_d * T_u;
					}
					// some charge is lost
					else
					{
						ev_b = ev_b + max_c_ev * eta_c_ev * T_u;
						b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
					}
				}
				// cannot discharge b, as it is empty
				else if (max_d < d)
				{
					// we need to discharge the ev battery
					double rest_d = d - max_d;
					// cout << "discharge ev with: " << rest_d << endl;
					// cout << "ev_b before discharge with: " << ev_b << endl;
					// cout << "b before evb discharge with: " << b << endl;
					if (rest_d > max_d_ev)
					{
						loss_events += 1;
						load_deficit += (rest_d - max_d_ev);

						ev_b = ev_b - max_d_ev * eta_d_ev * T_u;
						b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
						// cout << "ev_b after not full discharge: " << ev_b << endl;
					}
					else
					{
						// cout << "rest_d_ev before discharge: " << rest_d << endl;
						ev_b = ev_b - rest_d * eta_d_ev * T_u;
						// cout << "ev_b after  full discharge: " << ev_b << endl;

						b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
					}
				}
				// only use ev
				else
				{
					// cout << "ev not needed, b before =  " << b << endl;
					// cout << "only update stationary storage - before: " << b << endl;
					ev_b = ev_b + max_c_ev * eta_c * T_u;
					b = b - max_d * eta_d * T_u;
					// cout << "only update stationary storage - after: " << b << endl;

					// cout << "ev not needed, b after =  " << b << endl;
				}
			}
		}
		else
		{
			//cout << "ev = false with time: " << tt << endl;
			// at each time step either c or d is 0, which is why either max_c or max_d is 0 and it works out
			//cout << "b before update: " << b << endl;
			b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
			//cout << "b after update: " << b << endl;
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





/*still need to model that ev is fully charged when it leaves : muss load part verändern 
	- 7h für eine volle ladung 
	- 39 kwh
	- kommt halb voll an --> 3.5h ladung zwischen von 19,5kwh 

*/



// then add flexible initial ev battery state and load schedule

/*
Priority of improvements:

- better charging schedule
- variable batterie stand when you arrive
- lifestyle changes
- more sophisticated operating policy
- Maybe add loss metric for battery level of ev: how fully charged does it have to be when you leave the house ?


*/