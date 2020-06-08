import numpy as np

n = input("Enter the amount of primes: ")
index = input("Enter indices: ").split(" ")
index = np.array([int(x) for x in index])

max = np.max(index)

Primes, Numbers = [], []
count = 2

while len(Primes) <= max:
    if len(Primes)!=0:
        F = np.array([count%prime==0 for prime in Primes])
        if not True in F: Primes.append(count)
    else: Primes.append(count)
    #Numbers.append(count)
    count+=1
result = [Primes[x-1] for x in index]
print(result)
