def v_2(x):
    if x == 0:
        return 99999999999
    out = 0
    while x % 2 == 0:
        out += 1
        x //= 2
    return out


records = []
max_v2 = 0
n0, m0 = 0, 0
M = 1000
print(v_2(9**599 + 7))
# for n in range(M): #7k -> 30
#     if n % 10 == 0:
#         pass
#         #print(n/M)
#     for m in range(M):
#         v = v_2(9**m + 7**n)
#         if v > max_v2:
#             records.append([v, m, n])
#             max_v2 = v
#             n0, m0 = n, m
#
# for r in records:
#     print(r)
# print('max = ' + str(max_v2) + ' m = ' + str(m0) + ' n = ' + str(n0))
