
file = open("wfht1.txt")
for line in file:
    double = float(line)
    print(double)

    with open("wfht1_2.txt", "a") as file:
        file.write(str(double)+ '\n')
    file.close()