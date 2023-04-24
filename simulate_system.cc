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

// returns true if the ev needs to be charged at the current time step
bool latest_charge(int hour, double current_ev_soc, int t_next_dept){
	// compute pink line: ensure that ev is 80% charged when it leaves the house
	double diff = ev_goal_kWh - current_ev_soc;
	int hours_needed = ceil(diff / charging_rate);
	int latest_t = t_next_dept - hours_needed;
	if (hour == latest_t){
		return true;
	}
	else{
		return false;
	}
}

//operating policies
void giordano(){
}

/*
void stat_prioritized(){
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


void ev_prioritized(){
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

void round_robin(){
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

*/

// previously c_ev_d_b
void bidirectional(double ev_b, double b, double c, double d, double max_c, double max_d, double max_c_ev, double max_d_ev, bool needs_charge){
	/*
	we can charge the EV battery when there is extra charge left (and prioritise charging ev)
	we can discharge the EV battery when the stationary battery is empty 

	entweder is c positiv und wir chargen ev - sonst stationary wenn ev voll
	oder d is positiv und wir dischargen stationary - sonst ev - sonst grid 
	*/
	if (c > 0 && c > max_c_ev){
		// cannot store charge in ev as it is already full
		double rest_c = c - max_c_ev;
		if (rest_c < max_c){
			// can fully charge stationary with the rest
			//cout << "ev_b before update l .273: " << ev_b << endl;
			ev_b = ev_b + max_c_ev * eta_c_ev * T_u;
			//cout << "ev_b after update: " << ev_b << endl;
			b = b + rest_c * eta_c * T_u - max_d * eta_d * T_u;
		} else{
			// some charge is lost
			//cout << "ev_b before update l.281: " << ev_b << endl;
			ev_b = ev_b + max_c_ev * eta_c_ev * T_u;
			//cout << "ev_b after update: " << ev_b << endl;
			b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
		}
	}
	
	else if (d > 0 &&  max_d < d){
		// cannot discharge b, as it is empty. we need to discharge the ev battery
		double rest_d = d - max_d;
		if (rest_d > max_d_ev){
			// ev is empty too. need grid.
			loss_events += 1;
			load_deficit += (rest_d - max_d_ev);
			//cout << "ev_b before update l .297: " << ev_b << endl;
			ev_b = ev_b - max_d_ev * eta_d_ev * T_u;
			//cout << "ev_b after update: " << ev_b << endl;
			b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
		}
		else{
			// can discharge ev battery
			//cout << "ev_b before update l. 306: " << ev_b << endl;
			ev_b = ev_b - rest_d * eta_d_ev * T_u;
			//cout << "ev_b after update: " << ev_b << endl;
			b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
		}
	} else {
		// normal behaviour : charge ev or discharge stationary
		//cout << "ev_b before update l.315: " << ev_b << endl;
		ev_b = ev_b + max_c_ev * eta_c * T_u;
		//cout << "ev_b after update: " << ev_b << endl;
		b = b - max_d * eta_d * T_u;
	}
}

// previously c_ev_d_b_only
void unidirectional(double ev_b, double b, double c, double d, double max_c, double max_d, double max_c_ev, double max_d_ev, bool needs_charge){
	/*
	we can charge the EV battery when there is extra charge left (and prioritise charging ev)
	we can not discharge the EV battery

	entweder is c positiv und wir chargen ev - sonst stationary wenn ev voll 
	oder d is positiv und wir dischargen stationary - sonst grid
	*/
	if (c > 0 && c > max_c_ev){
		// cannot store charge in ev as it is already full
		double rest_c = c - max_c_ev;
		// we generated more charge than we can store in ev battery and will store the excess charge in the stationary
		if (rest_c < max_c){
			// can fully charge stationary with the rest
			ev_b = ev_b + max_c_ev * eta_c_ev * T_u;
			b = b + rest_c * eta_c * T_u - max_d * eta_d * T_u;
		} else{
			// some charge is lost: stationary and ev are both full
			ev_b = ev_b + max_c_ev * eta_c_ev * T_u;
			b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
		}
	}
	else if (d > 0 && max_d < d){
		// cannot discharge b, as it is empty. since uni-directional, we are not going to discharge the ev.
		loss_events += 1;

		load_deficit += (d - max_d);
		//cout << "increased load deficit = " << load_deficit << endl;

		b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
	} else{
		// normal behaviour : charge ev or discharge stationary
		ev_b = ev_b + max_c_ev * eta_c * T_u;
		b = b - max_d * eta_d * T_u;
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
	double current_ev_soc;
	bool needs_charge = false;

	int trace_length_solar = solar_trace.size();
	int trace_length_load = load_trace.size();
	int trace_length_ev = ev_trace.size();
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
	bool ev_presence = false;
	int counter = 0;

	loss_events = 0;
	load_deficit = 0;
	load_sum = 0;

	//for each of the (100) days in the sample
	for (int i = start_index; i < start_index + days_in_chunk; i++){

		//	cout << " DAY NUMBER  : " << i % start_index << endl;
			ev_trace_index = i % start_index + counter;

			// initalise arrays
			for (int k = 0; k < 24; k++){
				ev_at_home[k] = true;
			}

			// ----------------------------------------------------------------read EV inputs-----------------------------------
			num_trips = ev_trace[ev_trace_index];
			//cout << "num_trips : " << num_trips << endl;
			// next departures
			int array_length = num_trips + 1;
			int next_dept_arr[array_length];
			int next_dept_size = sizeof(next_dept_arr) / sizeof(int);
			// cout << "1 - length of next_dept_arr  is : " << next_dept_size << endl;

			for (int j = 0; j < num_trips; j++){
				ev_trace_index = ev_trace_index + 1;
				counter = counter + 1;
				t_dept = ev_trace[ev_trace_index];
				//cout << "t_dept : " << t_dept << endl;

				ev_trace_index = ev_trace_index + 1;
				counter = counter + 1;
				t_arr = ev_trace[ev_trace_index];
				//cout << "t_arr : " << t_arr << endl;

				ev_trace_index = ev_trace_index + 1;
				counter = counter + 1;
				soc_arr = ev_trace[ev_trace_index];
				//cout << "soc_Arr : " << soc_arr << endl;

				// next departure
				next_dept_arr[j] = t_dept;

				for (int h = t_dept; h < t_arr; h++){
					ev_at_home[h] = false;
				}
				// initalise soc array
				ev_soc[t_arr] = soc_arr * 60;
			}
			//----------------------------------------------------------------------------------------------------------------

			// if next num_dept is 0: counts number of consecutive days with no trips
			int pad = 0;
			// how do I get num trips from the next day?
			int ev_trace_index_2 = i % start_index + counter + 1;
			// cout << "ev_trace_index is = " << ev_trace_index  << endl;
			// cout << "next num trips that we read in = " << ev_trace[ev_trace_index_2 + pad ] << endl;

			while (ev_trace[ev_trace_index_2 + pad] == 0){
				//cout << "on the next day there are no trips" << endl;
				// on the next day there are no trips
				pad = pad + 1;
				// cout << "pad = " << pad<<endl;
			}

			next_dept_arr[num_trips] = ev_trace[ev_trace_index + 2 + pad];
			//cout << "start to loop through hours" << endl;
			pad = 0;
			// DELETE ME - DEBUGGING
			/*
			for (int p = 0; p < 24; p++){
				cout << "ev_at_home value at time" << p<< "is : " << ev_at_home[p] <<endl;
				cout << "ev_soc value at time" << p << "is : " << ev_soc[p] << endl;
			}

			cout << "length of next_dept_arr  is : " << next_dept_size << endl;
			for (int p = 0; p < next_dept_size; p++)
			{
				cout << "next_dept_arr value is : " << next_dept_arr[p] << endl;
			}

*/
			// ev is at home the whole day today
			if (num_trips == 0){
				// should ignore next_dept today
				//cout << "NO TRIP TODAY ev at home: " << endl;

				for (int t = 0; t < 24; t++){
					// wrap around to the start of the trace if we hit the end.

					index_t_solar = t % trace_length_solar;
					index_t_load = t % trace_length_load;
					load_sum += load_trace[index_t_load];
					//cout << "hour : " << t << endl;
					//cout << "ev_soc value at time" << t << "is : " << ev_soc[t] << endl;
					double load = load_trace[index_t_load];
					next_dept = next_dept_arr[0];

					//cout << "next departure is : " << next_dept << endl;
					if (ev_soc[t] == 0 && t == 0 && i == 0){
						// only runs for first monday 0h : assume soc is 60% charged initially
						ev_b = 36;
					}
					else{
						//	cout << "ev_b before update l .504: " << ev_b << endl;
						ev_b = ev_soc[t];
						// cout << "ev_b after update l.509: " << ev_b << endl;
					}

					//cout << "ev battery is " << ev_b << endl;
					ev_soc[(t + 1) % 24] = ev_b;

					
					// no need to charge the EV
					
					c = fmax(solar_trace[index_t_solar] * pv - load, 0);
					d = fmax(load - solar_trace[index_t_solar] * pv, 0);

					//cout << "c value is = " << c << endl;
					//cout << "d value is = " << d << endl;

					max_c = fmin(calc_max_charging(c, b), alpha_c);
					max_d = fmin(calc_max_discharging(d, b), alpha_d);

					//cout << "max_c value is = " << max_c << endl;
					//cout << "max_d value is = " << max_d << endl;
					max_c_ev = fmin(calc_max_charging_ev(c, ev_b, ev_goal_kWh), alpha_c_ev);
					max_d_ev = fmin(calc_max_discharging_ev(d, ev_b), alpha_d_ev);

					//cout << "max_c_ev value is = " << max_c_ev << endl;
					//cout << "max_d_ev value is = " << max_d_ev << endl;

					// call operating policy that we want
					// giordano();
					// bidirectional(ev_b, b, c, d, max_c, max_d, max_c_ev, max_d_ev, needs_charge);

					//cout << "stationary b before calling unidrectional() is = " << b << endl;

					 unidirectional(ev_b, b, c, d, max_c, max_d, max_c_ev, max_d_ev, needs_charge);

					//cout << "stationary b after calling unidrectional() is = " << b << endl;
					//cout << "ev_b after calling unidrectional() is = " << ev_b << endl;				
			}
		}

		// iterate through each hour to simulate charging behaviour
		if(num_trips > 0) {
			for (int t = 0; t < 24; t++){
			// wrap around to the start of the trace if we hit the end.
			
			index_t_solar = t % trace_length_solar;
			index_t_load = t % trace_length_load;
			load_sum += load_trace[index_t_load] ;
			//cout << "hour : " << t << endl;
			//cout << "ev_soc value at time" << t << "is : " << ev_soc[t] << endl;
			double load = load_trace[index_t_load];

			int trips_counter = 0;
			// TODO: this does not work, as we always have one next value for trips with multiple days
			if (next_dept == t){
				trips_counter = trips_counter + 1;
			}
			next_dept = next_dept_arr[trips_counter];
			//cout << "next departure is : " << next_dept << endl;
			
			if (ev_at_home[t] == true){
				// ev is at home
				if (ev_soc[t] == 0 && t == 0 && i == 0){
					// only runs for first monday 0h : assume soc is 60% charged initially
					ev_b = 36;
				}
				else{
					ev_b = ev_soc[t];
				}
				// cout << "ev battery is " << ev_b << endl;
				needs_charge = latest_charge(t, ev_b, next_dept);
				// cout << "needs charge : " << needs_charge << endl;

				if (needs_charge){
					// charge the EV
					//cout << "ev needs to be charged " << endl;
					c = fmax(solar_trace[index_t_solar] * pv - load - 7.4, 0);
					//TODO : macht die zeile sinn bei discharge?
					d = fmax(load + 7.4 - solar_trace[index_t_solar] * pv, 0);
				}
				else {
					// no need to charge the EV
					//cout << "NO CHARGE FOR EV " << endl;

					c = fmax(solar_trace[index_t_solar] * pv - load, 0);
					d = fmax(load - solar_trace[index_t_solar] * pv, 0);
				}
				
				//cout << "c value is = " << c << endl;
				//cout << "d value is = " << d << endl;

				max_c = fmin(calc_max_charging(c, b), alpha_c);
				max_d = fmin(calc_max_discharging(d, b), alpha_d);

				//cout << "max_c value is = " << max_c << endl;
				//cout << "max_d value is = " << max_d << endl;

				max_c_ev = fmin(calc_max_charging_ev(c, ev_b, ev_goal_kWh), alpha_c_ev);
				max_d_ev = fmin(calc_max_discharging_ev(d, ev_b), alpha_d_ev);

				//cout << "max_c_ev value is = " << max_c_ev << endl;
				//cout << "max_d_ev value is = " << max_d_ev << endl;

				// call operating policy that we want
				//giordano();
				//bidirectional(ev_b, b, c, d, max_c, max_d, max_c_ev, max_d_ev, needs_charge);

				//cout << "stationary b before calling unidrectional() is = " << b << endl;

				unidirectional(ev_b, b, c, d, max_c, max_d, max_c_ev, max_d_ev, needs_charge);

				//cout << "stationary b after calling unidrectional() is = " << b << endl;
				//cout << "ev_b after calling unidrectional() is = " << ev_b << endl;

				if (needs_charge){
					ev_b = ev_b + 7.4;
				}
				ev_soc[(t + 1) % 24] = ev_b;
			}
			else
			{
				// EV NOT AT HOME
				ev_b = 0.0;
				c = fmax(solar_trace[index_t_solar] * pv - load, 0);
				d = fmax(load - solar_trace[index_t_solar] * pv, 0);
				/*
				alpha_c = max charging rate of stationary battery;
				alpha_d = max discharging rate of stationary battery;
				max_c = max charging that we can apply to stationary battery this hour
				max_d = max discharging that we can apply to stationary battery this hour
				*/
				max_c = fmin(calc_max_charging(c, b), alpha_c);
				max_d = fmin(calc_max_discharging(d, b), alpha_d);
				b = b + max_c * eta_c * T_u - max_d * eta_d * T_u;
				
				if (max_d < d){
					loss_events += 1;
					load_deficit += (d - max_d);
					//cout << "increased load deficit = " << load_deficit << endl;
				}
			}

		}
		}
		}

	if (metric == 0) {
		// lolp
		double result = loss_events / ((100) * 1.0) ;
		cout << "RESULT NUMBER OF LOSS EVENTS" << result << endl;
		return result;
	} else {
		// metric == 1, eue
		cout << "RESULT load deficit = " << load_deficit << endl;
		cout << "RESULT total load = " << load_sum << endl;

		double result = load_deficit / (load_sum * 1.0);
		cout << "RESULT LOSS" << result <<endl;
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

		// for each value of cells, find the lowest pv that meets the epsilon loss constraint
		double loss = 0;
		while (true) {
			//cout << "-------call sim() with the following number of cells = "<< cells<<endl;
			// compute loss of current sizing. 
			loss = sim(load_trace, solar_trace, ev_trace, start_index, end_index, cells, lowest_feasible_pv - pv_step, b_0);
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
