# Calculate odds ratio

Neither = 8336
A.notB = 803
B.notA = 1754
Both = 52

total = sum(Neither, A.notB, B.notA, Both)

pA = (A.notB + Both) / total
pB = (B.notA + Both) / total

pA.notB = A.notB / total
pB.notA = B.notA / total

OR = (pA/pB.notA) / (pB/pA.notB)

log2(OR)
