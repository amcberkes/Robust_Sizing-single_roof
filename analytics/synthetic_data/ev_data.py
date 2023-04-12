import numpy as np
import matplotlib.pyplot as plt




'''
We need to generate a list of the following type:
t_start, t_end, soc_arrival

for a period of one year (365 days).

1 - Generate functions for 24h and let it run 365 days. distinguish between weekdays and weekends, 
so have different fucntions for that.

2 - for each WFH type, run the functions 365 times to get a time series, by putting in the corresponding parameters.


Other remarks:
- hourly time resolution : assume a minimum trip duration of one hour 
- charging rate: 


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
    dep_t = uniform_t(7,10)
    # arrival
    arr_t = uniform_t(dep_t + 6,dep_t + 9)
    # SOC
    soc = soc_arrival(0.6, 0.75)
    return dep_t, arr_t, soc

def t1_weekend():
    # departure
    dep_t = uniform_t(10,16)
    # arrival
    arr_t = uniform_t(dep_t + 1,dep_t + 5)
    # SOC
    soc = soc_arrival(0.5, 0.75)
    return dep_t, arr_t, soc

def t3_weekday():
    # departure
    dep_t = uniform_t(10,16)
    # arrival
    arr_t = uniform_t(dep_t + 1,dep_t + 5)
    # SOC
    soc = soc_arrival(0.5, 0.75)
    return dep_t, arr_t, soc

def t3_weekend():
    # departure
    dep_t = uniform_t(10,16)
    # arrival
    arr_t = uniform_t(dep_t + 1,dep_t + 5)
    # SOC
    soc = soc_arrival(0.5, 0.75)
    return dep_t, arr_t, soc

################### WFH T1 - commute to work every day #############################

def wfh_t1():
    weekday = 0

    while weekday < 365:

        if weekday%7 != 5 or weekday%7 != 6: 
            dep_t, arr_t, soc =  t1_weekday()
        else:
           dep_t, arr_t, soc =  t1_weekend()

        # open file
        with open("wfht1.txt", "a") as file:
            # set a list of lines to add:
            lines = [str(dep_t), str(arr_t), str(soc)]
            # write to file and add a separator
            file.writelines(s + '\n' for s in lines)
    weekday = weekday + 1

################### WFH T2 - commutes to work on Tuesdays and Thursdays #############################
# behave like T1 on weekday 2 and 4, and like t3 on the other days 

def wfh_t2():
    weekday = 0

    while weekday < 365:

        if weekday%7 == 1 or weekday%7 == 3: 
            dep_t, arr_t, soc =  t1_weekday()
        if weekday%7 == 0 or weekday%7 == 2 or weekday%7 == 4:
            dep_t, arr_t, soc =  t3_weekday()
        else:
            dep_t, arr_t, soc =  t3_weekend()
        # open file
        with open("wfht2.txt", "a") as file:
            # set a list of lines to add:
            lines = [str(dep_t), str(arr_t), str(soc)]
            # write to file and add a separator
            file.writelines(s + '\n' for s in lines)
    weekday = weekday + 1


################### WFH T3 #############################

def wfh_t3():
    weekday = 0

    while weekday < 365:

        if weekday%7 != 5 or weekday%7 != 6: 
            dep_t, arr_t, soc =  t3_weekday()
        else:
            dep_t, arr_t, soc =  t3_weekend()


        # open file
        with open("wfht3.txt", "a") as file:
            # set a list of lines to add:
            lines = [str(dep_t), str(arr_t), str(soc)]
            # write to file and add a separator
            file.writelines(s + '\n' for s in lines)
    weekday = weekday + 1

################################# call whatever function you want to run 
wfh_t2()
