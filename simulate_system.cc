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
double calc_max_charging(double power, double b_prev) {

	double step = power/30.0;

	for (double c = power; c >= 0; c -= step) {
		double upper_lim = a2_slope*(c/nominal_voltage_c) + a2_intercept;
		double b = b_prev + c*eta_c*T_u;
		if (b <= upper_lim) {
			// wenn battery nach chargen von c noch unter der upper limit ist, können wir die battery mit amount c chargen
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
		//ANAIS: decrease battery value because we are discharging
		if (b >= lower_lim) {
			//wenn battery nach discharging noch über dem lower limit wäre, können wir d amount dischargen, sonst muss d kleiner gemacht werden
			return d;
		}
	}
	return 0;
}


// Note: sim_year calls procedures calc_max_charging and calc_max_discharging.
// You could potentially speed up the computation by expanding these functions into sim_year
// to avoid procedure calls in this inner loop.
// FUNCTION RETURNS LOLP OR EUE LOSS FOR A SPECIFIC C AND B
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
	bool increased_storage= false;
	//EV battery. how much it is charged 
	double b_ev = 0.0;
	//for each time step in the sample/Scenario
	cout << " startindex = " << start_index << endl;
	cout << " start_time = " << start_index%24 << endl;
	cout << " end_index = " << end_index << endl;
	cout << " end_time = " << end_index%24 << endl;
	for (int t = start_index; t < end_index; t++) {
	//ANAIS: NEED TO MODIFY HERE!!!!!
	cout << " t in loop = " << t << endl;
	// wrap around to the start of the trace if we hit the end.
	index_t_solar = t % trace_length_solar;
	index_t_load = t % trace_length_load;
	// load sum is a cumulative sum of the past load demands in the sample
	load_sum += load_trace[index_t_load];
	// EV arrives at the house
	if (t % 24 == 23 || t % 24 == 0 || t % 24 == 1 || t % 24 == 2 || t % 24 == 3 || t % 24 == 4 || t % 24 == 5 || t % 24 == 6){
		cout << "in the increase load loop with t = " << t << endl;
		cout << "uhrzeit = " << t%24 << endl;
		// increase load by 1 kwh each hour between 23:00 - 07:00 to charge 8 Kwh in total
		//  first, calculate how much power is available for charging, and how much is needed to discharge
		//  c ist wie viel von der solar generation übrig bliebt nachdem load abgedeckt ist = wie viel wir batterie chargen können
		c = fmax(solar_trace[index_t_solar] * pv - load_trace[index_t_load] - 1, 0);
		// d ist wie viel vom load wir nicht mit der generierten solar energie abdecken können =  wie viel wir aus batterie ziehen müssen
		d = fmax(load_trace[index_t_load] + 1 - solar_trace[index_t_solar] * pv, 0);
	}

	else {
			// first, calculate how much power is available for charging, and how much is needed to discharge
			// c ist wie viel von der solar generation übrig bliebt nachdem load abgedeckt ist = wie viel wir batterie chargen können
			c = fmax(solar_trace[index_t_solar] * pv - load_trace[index_t_load], 0);
			// d ist wie viel vom load wir nicht mit der generierten solar energie abdecken können =  wie viel wir aus batterie ziehen müssen
			d = fmax(load_trace[index_t_load] - solar_trace[index_t_solar] * pv, 0);
	}
		//batterry arrives with 10kwh full, 8kwh need to be charged
	if (t % 24 == 18 || t % 24 == 19 || t % 24 == 20|| t % 24 == 21|| t % 24 == 22){
		if(!increased_storage){
			b_ev += 10;
			increased_storage = true;
		}
		
		// constrain the power
		// max_c ist die max an solar charge die wir in batterie storen können
		max_c = fmin(calc_max_charging(c, b + b_ev), alpha_c);
		cout << "max charge an solar die wir in batterie storen können = " << max_c << endl;
		// max_d is max an energie die wir aus batterie ziehen können
		max_d = fmin(calc_max_discharging(d, b + b_ev), alpha_d);
		cout << "max charge an energie die wir aus batterie ziehen können = " << max_d << endl;

		// update energy content of battery at time slot t
		//  mind. eins von beiden ist 0 weil wir entweder zu viel strom oder zu wenig strom produziert haben um zeitpunkt t
		cout << "batterie zustand vor update = " << b << endl;
		cout << "batterie zustand ev vor update = " << b_ev << endl;
		b = b + max_c * eta_c * T_u - max_d * eta_d * T_u + b_ev;
		cout << "batterie zustand nach update = " << b << endl;
	
		
	}
	else if (t % 24 == 0 || t % 24 == 1 || t % 24 == 2 || t % 24 == 3 || t % 24 == 4 || t % 24 == 5 || t % 24 == 6 || t % 24 == 7 ){
		if(increased_storage){
			b_ev = b_ev - 1;;
		} else{
	
			if (t % 24 == 0)
			{
				b_ev = 11;
				increased_storage = true;
			}
			if (t % 24 == 1)
			{
				b_ev = 12;
				increased_storage = true;
			}
			if (t % 24 == 2)
			{
				b_ev = 13;
				increased_storage = true;
			}
			if (t % 24 == 3)
			{
				b_ev = 14;
				increased_storage = true;
			}
			if (t % 24 == 4)
			{
				b_ev = 15;
				increased_storage = true;
			}
			if (t % 24 == 5)
			{
				b_ev = 16;
				increased_storage = true;
			}
			if (t % 24 == 6)
			{
				b_ev = 17;
				increased_storage = true;
			}
			if (t % 24 == 7)
			{
				b_ev= 18;
				increased_storage = true;
			}
		}
		// constrain the power
		// max_c ist die max an solar charge die wir in batterie storen können
		max_c = fmin(calc_max_charging(c, b + b_ev), alpha_c);
		cout << "max charge an solar die wir in batterie storen können = " << max_c << endl;
		// max_d is max an energie die wir aus batterie ziehen können
		max_d = fmin(calc_max_discharging(d, b + b_ev), alpha_d);
		cout << "max charge an energie die wir aus batterie ziehen können = " << max_d << endl;

		// update energy content of battery at time slot t
		//  mind. eins von beiden ist 0 weil wir entweder zu viel strom oder zu wenig strom produziert haben um zeitpunkt t
		cout << "batterie zustand vor update = " << b << endl;
		cout << "batterie zustand ev vor update = " << b_ev << endl;
		b = b + max_c * eta_c * T_u - max_d * eta_d * T_u + b_ev;
		cout << "batterie zustand nach update = " << b << endl;
		// if we didnt get to discharge as much as we wanted, there is a loss
		}
	else {
		// constrain the power
		// max_c ist die max an solar charge die wir in batterie storen können
		max_c = fmin(calc_max_charging(c, b + b_ev), alpha_c);
		cout << "max charge an solar die wir in batterie storen können = " << max_c << endl;
		// max_d is max an energie die wir aus batterie ziehen können
		max_d = fmin(calc_max_discharging(d, b + b_ev), alpha_d);
		cout << "max charge an energie die wir aus batterie ziehen können = " << max_d << endl;

		// update energy content of battery at time slot t
		//  mind. eins von beiden ist 0 weil wir entweder zu viel strom oder zu wenig strom produziert haben um zeitpunkt t
		cout << "batterie zustand vor update = " << b << endl;
		cout << "batterie zustand ev vor update = " << b_ev << endl;
		b = b + max_c * eta_c * T_u - max_d * eta_d * T_u + b_ev;
		cout << "batterie zustand nach update = " << b << endl;
		// if we didnt get to discharge as much as we wanted, there is a loss
	}
	

		
		if (max_d < d) {
			loss_events += 1;
			load_deficit += (d - max_d);
		}
	}
	cout << "finished looping over sample " << endl;
	cout << "number of loss events for this sample= " << loss_events << endl;
	//this function returns the loss 

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
// RETURNS THE SIZING CURVE FOR THIS SCENARIO 
vector <SimulationResult> simulate(vector <double> &load_trace, vector <double> &solar_trace, int start_index, int end_index, double b_0) {

	// first, find the lowest value of cells that will get us epsilon loss when the PV is maximized
	// use binary search
	double cells_U = cells_max;
	double cells_L = cells_min;
	double mid_cells = 0.0;
	double loss = 0.0;
//PART 1 OF COMPUTING SIZING CURVE: FOR CURRENT C, FIND SMALLEST B THAT SATISFIES QOS

	while (cells_U - cells_L > cells_step) {

		mid_cells = (cells_L + cells_U) / 2.0;
		// call function from above
		loss = sim(load_trace, solar_trace, start_index, end_index, mid_cells, pv_max, b_0);

		//cout << "sim result with " << a2_intercept << " kWh and " << pv_max << " pv: " << loss << endl;
		// we want loss to be < epsilon 
		if (loss > epsilon) {
			// das wollen wir nicht deshalb muss battery size stark increased werden
			cells_L = mid_cells;
		} else {
		 	// (loss <= epsilon)
			// können schauen ob wir weniger cells als upper limit benutzen können
			cells_U = mid_cells;
		}
	}

	// set the starting number of battery cells to be the upper limit that was converged on
	double starting_cells = cells_U;
	// double B_inv; // cost per cell and double PV_inv; // cost per unit (kW) of PV
	double starting_cost = B_inv*starting_cells + PV_inv * pv_max;
	double lowest_feasible_pv = pv_max;


	double lowest_cost = starting_cost;
	double lowest_B = starting_cells*kWh_in_one_cell;
	double lowest_C = pv_max;


//curve is a vector of (b,c,cost) triplets
	vector <SimulationResult> curve;
	curve.push_back(SimulationResult(starting_cells*kWh_in_one_cell, lowest_feasible_pv, starting_cost));
	//cout << "starting cells: " << starting_cells << endl;

//PART 2 OF COMPUTING SIZING CURVE: DECREASE C SO MUCH DASS ES GERADE SO NOCH KLAPPT MIT DER CUURENT CELLS NUMBER

//cells_max is max amount of battery in kwh abailable, specified by user in terminal
	for (double cells = starting_cells; cells <= cells_max; cells += cells_step) {

		// for each value of cells, find the lowest pv that meets the epsilon loss constraint
		double loss = 0;
		while (true) {
			
			loss = sim(load_trace, solar_trace, start_index, end_index, cells, lowest_feasible_pv - pv_step, b_0);

			if (loss < epsilon) {
				//decrease pv if we are still in the green area with the current cells available
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
