import numpy as np
import matplotlib.pyplot as plt




'''
We need to generate a list of the following type:
t_start, t_end, soc_arrival

for a period of 4 years (365 days * 4).

1 - Generate functions for 24h and let it run 1460 days. distinguish between weekdays and weekends, 
so have different fucntions for that.

2 - for each WFH type, run the functions 1460 times to get a time series, by putting in the corresponding parameters.


Other remarks:
- hourly time resolution : assume a minimum trip duration of one hour 



'''


# function to generate uniformly distributed time (int format)
def uniform_t(a,b):
    t = round(np.random.uniform(a,b))
    return t


# function for SOC
def soc_arrival(a,b):
    soc = np.random.uniform(a,b)
    return soc

def t1_weekday():
     # departure
    dep_t = uniform_t(7.0,10.0) 
    # arrival
    arr_t = uniform_t(dep_t + 7.0,dep_t + 10.0) 
    # SOC
    soc = soc_arrival(0.6, 0.75)
    return dep_t, arr_t, soc


def t3_weekday():
    # departure
    dep_t = uniform_t(8.0,18.0)
    # arrival
    arr_t = uniform_t(dep_t + 1.0,dep_t + 2.0)
    # SOC
    soc = soc_arrival(0.70, 0.75)
    return dep_t, arr_t, soc

def tx_weekend():
    # departure
    dep_t = uniform_t(9.0,18.0)
    # arrival
    arr_t = uniform_t(dep_t + 1.0,dep_t + 2.0)
    # SOC
    soc = soc_arrival(0.5, 0.75)
    return dep_t, arr_t, soc

################### WFH T1 - commute to work every day #############################

def wfh_t1():
    weekday = 0
    while weekday < 1460:
        print("weekday = " + str(weekday))

        #weekend
        if weekday%7 == 5 or weekday%7 == 6: 
            num_trips = uniform_t(0.0,2.0)
            double_trips = num_trips * 1.0
            if num_trips == 0:
                with open("wfht1.txt", "a") as file:
                    file.write(str(double_trips)+ '\n')
                file.close()
            if num_trips == 1:
                dep_t, arr_t, soc =  tx_weekend()
                 # open file
                with open("wfht1.txt", "a") as file:
                 # set a list of lines to add:
                    lines = [str(num_trips), str(dep_t), str(arr_t), str(soc)]
                  # write to file and add a separator
                    file.writelines(s + '\n' for s in lines)
                file.close()
            if num_trips > 1:
                dep_t, arr_t, soc =  tx_weekend()
                with open("wfht1.txt", "a") as file:
                    lines = [str(num_trips),str(dep_t), str(arr_t), str(soc)]
                    file.writelines(s + '\n' for s in lines)

                    for i in range(num_trips - 1):
                        dep_t_next, arr_t_next, soc_next =  tx_weekend()
                        if dep_t_next > arr_t:
                            lines = [str(dep_t_next), str(arr_t_next), str(soc_next)]
                            file.writelines(s + '\n' for s in lines)
                        else:
                            if i == 0:
                                dep_t_next, arr_t_next, soc_next =  min(21, arr_t + 1), min(22, arr_t + 2), soc_arrival(0.5, 0.75)
                            if i == 1:
                                dep_t_next, arr_t_next, soc_next =  min(22, arr_t + 1), min(23, arr_t + 2), soc_arrival(0.5, 0.75)
                            lines = [str(dep_t_next), str(arr_t_next), str(soc_next)]
                            file.writelines(s + '\n' for s in lines)
                        arr_t = arr_t_next
                file.close()
        # weekday commute
        if weekday%7 == 0 or weekday%7 == 1 or weekday%7 == 2 or weekday%7 == 3 or weekday%7 == 4:
            num_trips = 1
            dep_t, arr_t, soc =  t1_weekday()
             # open file
            with open("wfht1.txt", "a") as file:
                # set a list of lines to add:
                lines = [str(num_trips),str(dep_t), str(arr_t), str(soc)]
                # write to file and add a separator
                file.writelines(s + '\n' for s in lines)
            file.close()
        weekday = weekday + 1

################### WFH T2 - commutes to work on Tuesdays and Thursdays #############################
# behave like T1 on weekday 2 and 4, and like t3 on the other days 

def wfh_t2():
    weekday = 0
    
    while weekday < 1460:
        print("weekday = " + str(weekday))
        # commute to work - Tuesday and Thursday
        if weekday%7 == 1 or weekday%7 == 3: 
            num_trips = 1
            dep_t, arr_t, soc =  t1_weekday()
             # open file
            with open("wfht2.txt", "a") as file:
                # set a list of lines to add:
                lines = [str(num_trips),str(dep_t), str(arr_t), str(soc)]
                # write to file and add a separator
                file.writelines(s + '\n' for s in lines)
            file.close()
        #WFH - Monday, Wednesday, Friday
        if weekday%7 == 0 or weekday%7 == 2 or weekday%7 == 4:
            num_trips = uniform_t(0.0,1.0)
            if num_trips == 0:
                with open("wfht2.txt", "a") as file:
                    file.write(str(num_trips)+ '\n')
                file.close()
            if num_trips == 1:
                dep_t, arr_t, soc =  t3_weekday()
                 # open file
                with open("wfht2.txt", "a") as file:
                 # set a list of lines to add:
                    lines = [str(num_trips),str(dep_t), str(arr_t), str(soc)]
                  # write to file and add a separator
                    file.writelines(s + '\n' for s in lines)
                file.close()
            #minmum 2 trips this day
            if num_trips > 1:
                dep_t, arr_t, soc =  t3_weekday()
                with open("wfht2.txt", "a") as file:
                    lines = [str(num_trips),str(dep_t), str(arr_t), str(soc)]
                    file.writelines(s + '\n' for s in lines)

                    for i in range(num_trips - 1):
                        dep_t_next, arr_t_next, soc_next =  t3_weekday()
                        if dep_t_next > arr_t:
                            lines = [str(dep_t_next), str(arr_t_next), str(soc_next)]
                            file.writelines(s + '\n' for s in lines)
                        else:
                            if i == 0:
                                dep_t_next, arr_t_next, soc_next =  min(21, arr_t + 1), min(22, arr_t + 2), soc_arrival(0.5, 0.75)
                            lines = [str(dep_t_next), str(arr_t_next), str(soc_next)]
                            file.writelines(s + '\n' for s in lines)
                        arr_t = arr_t_next
                file.close()
                    
        # weekend
        if weekday%7 == 5 or weekday%7 == 6:
            num_trips = uniform_t(0.0,2.0)
            if num_trips == 0:
                with open("wfht2.txt", "a") as file:
                    file.write(str(num_trips)+ '\n')
                file.close()
            if num_trips == 1:
                dep_t, arr_t, soc =  tx_weekend()
                 # open file
                with open("wfht2.txt", "a") as file:
                 # set a list of lines to add:
                    lines = [str(num_trips), str(dep_t), str(arr_t), str(soc)]
                  # write to file and add a separator
                    file.writelines(s + '\n' for s in lines)
                file.close()
            if num_trips > 1:
                dep_t, arr_t, soc =  tx_weekend()
                with open("wfht2.txt", "a") as file:
                    lines = [str(num_trips),str(dep_t), str(arr_t), str(soc)]
                    file.writelines(s + '\n' for s in lines)

                    for i in range(num_trips - 1):
                        dep_t_next, arr_t_next, soc_next =  tx_weekend()
                        if dep_t_next > arr_t:
                            lines = [str(dep_t_next), str(arr_t_next), str(soc_next)]
                            file.writelines(s + '\n' for s in lines)
                        else:
                            if i == 0:
                                dep_t_next, arr_t_next, soc_next =  min(21, arr_t + 1), min(22, arr_t + 2), soc_arrival(0.5, 0.75)
                            if i == 1:
                                dep_t_next, arr_t_next, soc_next =  min(22, arr_t + 1), min(23, arr_t + 2), soc_arrival(0.5, 0.75)
                            lines = [str(dep_t_next), str(arr_t_next), str(soc_next)]
                            file.writelines(s + '\n' for s in lines)
                        arr_t = arr_t_next
                file.close()
        weekday = weekday + 1


################### WFH T3 #############################

def wfh_t3():
    weekday = 0

    while weekday < 1460:
        print("weekday = " + str(weekday))

        #weekend
        if weekday%7 == 5 or weekday%7 == 6: 
            num_trips = uniform_t(0.0,2.0)
            if num_trips == 0:
                with open("wfht3.txt", "a") as file:
                    file.write(str(num_trips)+ '\n')
                file.close()
            if num_trips == 1:
                dep_t, arr_t, soc =  tx_weekend()
                 # open file
                with open("wfht3.txt", "a") as file:
                 # set a list of lines to add:
                    lines = [str(num_trips),str(dep_t), str(arr_t), str(soc)]
                  # write to file and add a separator
                    file.writelines(s + '\n' for s in lines)
                file.close()
            if num_trips > 1:
                dep_t, arr_t, soc =  tx_weekend()
                with open("wfht3.txt", "a") as file:
                    lines = [str(num_trips),str(dep_t), str(arr_t), str(soc)]
                    file.writelines(s + '\n' for s in lines)

                    for i in range(num_trips - 1):
                        dep_t_next, arr_t_next, soc_next =  tx_weekend()
                        if dep_t_next > arr_t:
                            lines = [str(dep_t_next), str(arr_t_next), str(soc_next)]
                            file.writelines(s + '\n' for s in lines)
                        else:
                            if i == 0:
                                dep_t_next, arr_t_next, soc_next =  min(21, arr_t + 1), min(22, arr_t + 2), soc_arrival(0.5, 0.75)
                            if i == 1:
                                dep_t_next, arr_t_next, soc_next =  min(22, arr_t + 1), min(23, arr_t + 2), soc_arrival(0.5, 0.75)
                            lines = [str(dep_t_next), str(arr_t_next), str(soc_next)]
                            file.writelines(s + '\n' for s in lines)
                        arr_t = arr_t_next
                file.close()
        
        #WFH
        if weekday%7 == 0 or weekday%7 == 1 or weekday%7 == 2 or weekday%7 == 3 or weekday%7 == 4 :
            num_trips = uniform_t(0.0,1.0)
            if num_trips == 0:
                with open("wfht3.txt", "a") as file:
                    file.write(str(num_trips)+ '\n')
                file.close()
            if num_trips == 1:
                dep_t, arr_t, soc =  t3_weekday()
                 # open file
                with open("wfht3.txt", "a") as file:
                 # set a list of lines to add:
                    lines = [str(num_trips), str(dep_t), str(arr_t), str(soc)]
                  # write to file and add a separator
                    file.writelines(s + '\n' for s in lines)
                file.close()
            # min 2 trips
            if num_trips > 1:
                dep_t, arr_t, soc =  t3_weekday()
                with open("wfht3.txt", "a") as file:
                    lines = [str(num_trips),str(dep_t), str(arr_t), str(soc)]
                    file.writelines(s + '\n' for s in lines)

                    for i in range(num_trips - 1):
                        dep_t_next, arr_t_next, soc_next =  t3_weekday()
                        if dep_t_next > arr_t:
                            lines = [str(dep_t_next), str(arr_t_next), str(soc_next)]
                            file.writelines(s + '\n' for s in lines)
                        else:
                            if i == 0:
                                dep_t_next, arr_t_next, soc_next =  min(21, arr_t + 1), min(22, arr_t + 2), soc_arrival(0.5, 0.75)
                            lines = [str(dep_t_next), str(arr_t_next), str(soc_next)]
                            file.writelines(s + '\n' for s in lines)
                        arr_t = arr_t_next
                file.close()
        weekday = weekday + 1
        


################################# call whatever function you want to run 
wfh_t1()
