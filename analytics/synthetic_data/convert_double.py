
file = open("wfht2_v2.txt")
for line in file:
    double = float(line)
    print(double)

    with open("wfht2_v2_d2.txt", "a") as file:
        file.write(str(double)+ '\n')
    file.close()