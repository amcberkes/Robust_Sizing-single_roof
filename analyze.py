#import matplotlib.pyplot as plt

b = []
c = []
results = []
with open('chebyonc_woev.txt') as f:
    lines = f.readlines()
    
    for x in lines:
        splitted = x.split()
        print(splitted)
        b.append(splitted[0])
        c.append(splitted[1])
#plt.plot(b, c)  # Plot the chart
#plt.xlim(50, 140)
#plt.xlim(0, 140)
#plt.show() 
    for x in c:
        y =  x.replace(".",",")
        print(y)
        results.append(y)

        with open('chebyonc_woevcomma.txt.txt', 'a') as f:
            f.write(y  +"\n")

