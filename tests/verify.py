import numpy as np
from functools import cmp_to_key

N = int(input())
arr = [[0 for _ in range(N)] for _ in range(N)]
for i in range(N):
    nums = input().split()
    for j in range(0, 2*N, 2):
        arr[i][int(j/2)] = complex(float(nums[j]), float(nums[j+1]))

arr = np.array(arr)

def compare(a, b):
    return abs(b) - abs(a)

eigenvalues, trash = np.linalg.eig(arr)

for eig in sorted(eigenvalues, key=cmp_to_key(compare)):
    print("{:.6f} {:.6f}".format(eig.real, eig.imag))