x = [2,3,4,6,4,56,78,4,4,6,3,2,6,7,4,3,1,34,56,7,3,7,9,3,12,43,5,7,9]
lenght_window = 8
for i in range(0,len(x)-lenght_window+1):
    print(x[i:i+lenght_window],len(x[i:i+lenght_window]))