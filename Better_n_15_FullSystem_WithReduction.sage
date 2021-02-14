import scipy.special
from itertools import combinations 
import timeit
import random

start = timeit.timeit()
#define the paramaters

#how many clocks are required
no_clocks = 5000

#the size of the register
state_size = 15

#how many output bits are needed
no_output_bits=5000

#degree of the system of equations to be solved
degree = 2

#the number of monomials of degree 1->max degree
#the variables 'no_monomials' will be the number of linearisation variables
#'scipy.special.comb' provides the binomial coefficient

s_out = [0]*state_size
b_out = [0]*state_size

for i in range(state_size):
    s_out[i] = random.randint(0,1)
    b_out[i] = random.randint(0,1)
    
print(s_out)
print(b_out)

add=[0]*no_clocks
S_out = s_out+add
B_out = b_out+add



O = [0]*(no_clocks-state_size)
for i in range(state_size,no_clocks):
    O[i-state_size] =mod((S_out[i-(state_size-13)]+S_out[i-(state_size-1)]*S_out[i-(state_size-9)]+S_out[i-(state_size-4)]*S_out[i-(state_size-9)]+S_out[i-(state_size-9)]*S_out[i-(state_size-13)])*B_out[i-(state_size-4)]+(S_out[i-(state_size-4)]+S_out[i-(state_size-1)]*S_out[i-(state_size-13)]+S_out[i-(state_size-9)]*S_out[i-(state_size-13)]+S_out[i-(state_size-1)]*S_out[i-(state_size-4)]*S_out[i-(state_size-9)]+S_out[i-(state_size-1)]*S_out[i-(state_size-9)]*S_out[i-(state_size-13)]),2)
    B_out[i] = mod(S_out[i-(state_size-0 )]+ B_out[i-(state_size-0 )]+ B_out[i-(state_size-13 )]+ B_out[i-(state_size-2)]*B_out[i-(state_size-14 )]+ B_out[i-(state_size-3)]*B_out[i-(state_size-5 )],2)
    S_out[i] = mod((S_out[i-15]+S_out[i-8]),2)


no_monomials = 0
for i in range(1,degree+1):
    no_monomials = no_monomials + binomial(state_size,i)

#define polynomial ring in which the the INITIAL STATE, OUTPUT and LINEARISATION variables will be defined
R=BooleanPolynomialRing(state_size+no_monomials,'x')
#Inject these variables for use
R.inject_variables()

#define initial state
s = [0]*state_size
for i in range(state_size):
    s[i] = eval("x" + str(i))

# print(s)
#define output
# Y = [0]*no_output_bits
# for i in range(no_output_bits):
#     Y[i] = eval("x" + str(state_size+i))

# print(Y)
#We can consider 's' as a register that will be shifted for a specified number of times. Thus, there will exist more entries than simply those defined in 's'.
#We therefore initialise a vector 'add' that will be appended to s to hold all the update bits (one vector for equations with respect to keystream bits S and one with respect the the linear approximation$
Alpha = [1940 , 1938 , 1935 , 1933 , 1927 , 1922 , 1921 , 1916 , 1914 , 1913 , 1912 , 1910 , 1909 , 1908 , 1906 , 1903 , 1902 , 1899 , 1895 , 1892 , 1891 , 1890 , 1888 , 1882 , 1880 , 1877 , 1875 , 1871 , 1867 , 1866 , 1865 , 1864 , 1863 , 1862 , 1860 , 1859 , 1854 , 1851 , 1850 , 1848 , 1847 , 1845 , 1844 , 1843 , 1842 , 1840 , 1839 , 1833 , 1831 , 1830 , 1828 , 1827 , 1823 , 1822 , 1821 , 1820 , 1819 , 1814 , 1812 , 1811 , 1810 , 1809 , 1806 , 1805 , 1804 , 1803 , 1800 , 1799 , 1797 , 1796 , 1795 , 1791 , 1789 , 1788 , 1784 , 1778 , 1777 , 1773 , 1770 , 1769 , 1768 , 1765 , 1764 , 1762 , 1761 , 1759 , 1755 , 1753 , 1750 , 1749 , 1747 , 1746 , 1744 , 1743 , 1741 , 1740 , 1739 , 1737 , 1734 , 1727 , 1726 , 1724 , 1722 , 1721 , 1720 , 1716 , 1715 , 1712 , 1711 , 1710 , 1709 , 1708 , 1706 , 1705 , 1704 , 1702 , 1701 , 1695 , 1692 , 1691 , 1686 , 1685 , 1678 , 1673 , 1669 , 1666 , 1665 , 1664 , 1663 , 1659 , 1657 , 1654 , 1650 , 1649 , 1648 , 1646 , 1644 , 1643 , 1642 , 1638 , 1637 , 1636 , 1633 , 1632 , 1631 , 1628 , 1627 , 1625 , 1623 , 1621 , 1620 , 1618 , 1617 , 1616 , 1614 , 1612 , 1610 , 1609 , 1607 , 1606 , 1605 , 1604 , 1601 , 1597 , 1596 , 1595 , 1590 , 1585 , 1583 , 1582 , 1581 , 1579 , 1577 , 1575 , 1574 , 1573 , 1572 , 1570 , 1569 , 1568 , 1566 , 1565 , 1563 , 1561 , 1556 , 1555 , 1554 , 1551 , 1548 , 1543 , 1541 , 1540 , 1539 , 1536 , 1535 , 1534 , 1530 , 1527 , 1526 , 1524 , 1523 , 1522 , 1521 , 1520 , 1517 , 1513 , 1512 , 1509 , 1507 , 1506 , 1504 , 1501 , 1498 , 1496 , 1495 , 1492 , 1491 , 1484 , 1483 , 1482 , 1481 , 1479 , 1478 , 1476 , 1475 , 1474 , 1473 , 1469 , 1468 , 1467 , 1465 , 1464 , 1463 , 1461 , 1460 , 1459 , 1457 , 1456 , 1455 , 1454 , 1453 , 1450 , 1448 , 1447 , 1443 , 1441 , 1440 , 1438 , 1435 , 1433 , 1431 , 1429 , 1418 , 1416 , 1415 , 1414 , 1413 , 1412 , 1411 , 1408 , 1407 , 1405 , 1404 , 1401 , 1400 , 1399 , 1396 , 1395 , 1394 , 1393 , 1388 , 1387 , 1383 , 1382 , 1380 , 1379 , 1376 , 1375 , 1374 , 1373 , 1371 , 1370 , 1369 , 1368 , 1360 , 1359 , 1358 , 1357 , 1356 , 1355 , 1354 , 1351 , 1349 , 1348 , 1343 , 1342 , 1341 , 1339 , 1336 , 1334 , 1332 , 1328 , 1327 , 1325 , 1324 , 1323 , 1320 , 1319 , 1317 , 1316 , 1313 , 1311 , 1310 , 1307 , 1306 , 1305 , 1302 , 1300 , 1299 , 1297 , 1296 , 1295 , 1293 , 1287 , 1286 , 1282 , 1280 , 1278 , 1277 , 1274 , 1272 , 1271 , 1270 , 1269 , 1268 , 1266 , 1265 , 1264 , 1263 , 1260 , 1258 , 1256 , 1255 , 1251 , 1250 , 1249 , 1246 , 1245 , 1244 , 1243 , 1241 , 1239 , 1236 , 1235 , 1234 , 1233 , 1231 , 1226 , 1224 , 1223 , 1222 , 1221 , 1220 , 1219 , 1218 , 1216 , 1215 , 1214 , 1211 , 1210 , 1208 , 1207 , 1206 , 1205 , 1203 , 1202 , 1201 , 1199 , 1198 , 1189 , 1188 , 1187 , 1183 , 1182 , 1180 , 1179 , 1177 , 1175 , 1174 , 1171 , 1170 , 1169 , 1168 , 1167 , 1166 , 1164 , 1163 , 1161 , 1159 , 1158 , 1157 , 1156 , 1155 , 1154 , 1151 , 1148 , 1146 , 1143 , 1142 , 1141 , 1135 , 1134 , 1133 , 1132 , 1126 , 1123 , 1122 , 1120 , 1119 , 1116 , 1113 , 1111 , 1108 , 1101 , 1097 , 1095 , 1093 , 1092 , 1091 , 1090 , 1088 , 1086 , 1083 , 1082 , 1080 , 1078 , 1075 , 1074 , 1073 , 1070 , 1066 , 1064 , 1059 , 1057 , 1056 , 1052 , 1049 , 1046 , 1042 , 1041 , 1038 , 1034 , 1033 , 1027 , 1023 , 1021 , 1020 , 1017 , 1016 , 1015 , 1014 , 1013 , 1012 , 1008 , 1004 , 1001 , 998 , 997 , 996 , 995 , 993 , 991 , 989 , 988 , 987 , 986 , 982 , 981 , 976 , 973 , 972 , 970 , 969 , 966 , 961 , 952 , 951 , 948 , 945 , 944 , 943 , 941 , 939 , 937 , 936 , 931 , 930 , 926 , 925 , 917 , 916 , 914 , 913 , 912 , 906 , 905 , 904 , 903 , 900 , 899 , 898 , 892 , 891 , 887 , 885 , 884 , 883 , 882 , 879 , 875 , 874 , 870 , 868 , 866 , 863 , 858 , 857 , 856 , 855 , 850 , 845 , 842 , 840 , 838 , 837 , 835 , 834 , 832 , 831 , 828 , 827 , 823 , 822 , 821 , 820 , 818 , 816 , 815 , 814 , 806 , 804 , 803 , 802 , 801 , 799 , 798 , 796 , 795 , 793 , 791 , 789 , 788 , 787 , 786 , 785 , 783 , 781 , 780 , 779 , 778 , 776 , 775 , 774 , 773 , 768 , 767 , 765 , 764 , 763 , 760 , 757 , 756 , 755 , 754 , 752 , 750 , 749 , 748 , 745 , 743 , 742 , 739 , 738 , 737 , 736 , 735 , 734 , 733 , 728 , 723 , 722 , 719 , 717 , 716 , 714 , 712 , 711 , 707 , 704 , 703 , 702 , 701 , 700 , 699 , 696 , 695 , 694 , 692 , 688 , 683 , 681 , 680 , 675 , 673 , 671 , 668 , 664 , 663 , 662 , 660 , 659 , 656 , 655 , 652 , 650 , 649 , 645 , 641 , 638 , 637 , 633 , 632 , 631 , 628 , 626 , 624 , 621 , 619 , 615 , 613 , 612 , 610 , 607 , 605 , 604 , 600 , 597 , 596 , 595 , 591 , 589 , 588 , 587 , 586 , 585 , 584 , 583 , 582 , 580 , 579 , 578 , 574 , 572 , 571 , 568 , 564 , 562 , 560 , 559 , 558 , 556 , 553 , 550 , 549 , 548 , 544 , 542 , 536 , 528 , 521 , 520 , 517 , 516 , 514 , 508 , 507 , 506 , 505 , 504 , 503 , 501 , 499 , 496 , 493 , 491 , 490 , 489 , 486 , 484 , 482 , 481 , 478 , 475 , 474 , 473 , 472 , 467 , 466 , 463 , 461 , 460 , 459 , 458 , 452 , 451 , 449 , 448 , 446 , 445 , 444 , 442 , 441 , 438 , 437 , 433 , 429 , 427 , 426 , 425 , 424 , 423 , 421 , 419 , 418 , 417 , 416 , 413 , 411 , 409 , 407 , 406 , 404 , 398 , 397 , 395 , 394 , 392 , 391 , 390 , 389 , 388 , 385 , 382 , 379 , 374 , 371 , 370 , 368 , 367 , 366 , 364 , 363 , 362 , 357 , 355 , 352 , 351 , 350 , 349 , 346 , 345 , 343 , 340 , 338 , 336 , 335 , 333 , 332 , 330 , 329 , 327 , 326 , 325 , 322 , 320 , 316 , 313 , 311 , 310 , 308 , 307 , 303 , 302 , 301 , 296 , 289 , 286 , 284 , 283 , 280 , 276 , 270 , 269 , 268 , 265 , 263 , 262 , 261 , 260 , 257 , 255 , 251 , 249 , 248 , 247 , 245 , 244 , 243 , 242 , 241 , 238 , 234 , 232 , 228 , 227 , 225 , 224 , 223 , 222 , 221 , 220 , 219 , 213 , 212 , 210 , 206 , 204 , 203 , 201 , 200 , 199 , 195 , 194 , 193 , 192 , 190 , 188 , 187 , 182 , 181 , 180 , 177 , 176 , 175 , 173 , 169 , 167 , 165 , 163 , 162 , 161 , 160 , 159 , 154 , 152 , 150 , 147 , 144 , 143 , 142 , 140 , 138 , 137 , 130 , 128 , 127 , 126 , 124 , 119 , 118 , 117 , 115 , 114 , 113 , 111 , 109 , 108 , 107 , 105 , 104 , 98 , 97 , 96 , 91 , 89 , 87 , 85 , 84 , 83 , 80 , 79 , 77 , 76 , 63 , 59 , 57 , 56 , 55 , 54 , 53 , 52 , 51 , 50 , 48 , 47 , 46 , 43 , 40 , 39 , 38 , 36 , 35 , 25 , 24 , 23 , 20 , 18 , 17 , 15 , 13 , 12 , 11 , 9 , 5 , 0]

add=[0]*no_clocks
S = s+add
Final_System = [0]*(no_clocks-state_size-max(Alpha))

#Generate equations for update using output bits
# S[15] = S[0]+S[7]
for i in range(len(s),no_clocks):
     S[i] = (S[i-15]+S[i-8])

for j in range(no_clocks-state_size-max(Alpha)):
    for i in Alpha:
    #s1s2z=s1s2s9
    #Here we are going to build two systems; one containing z, the other not
           #Left[j]  = Left[j]+(S[(i+state_size)-(state_size-3)]*S[(i+state_size)-(state_size-1)]*S[(i+state_size)-(state_size-6 )]+ S[(i+state_size)-(state_size-1)]*S[(i+state_size)-(state_size-10)]*S[(i+state_size)-(state_size-6)])
            Final_System [j] = Final_System[j]+S[(i+state_size+j)-(state_size-4)]*S[(i+state_size+j)-(state_size-1)]*S[(i+state_size+j)-(state_size-13)]*S[(i+state_size+j)-(state_size-9 )]+ S[(i+state_size+j)-(state_size-4)]*S[(i+state_size+j)-(state_size-1)]*S[(i+state_size+j)-(state_size-9 )]+ S[(i+state_size+j)-(state_size-4)]*S[(i+state_size+j)-(state_size-13)]*S[(i+state_size+j)-(state_size-9 )]+ S[(i+state_size+j)-(state_size-4)]*S[(i+state_size+j)-(state_size-13 )]+ S[(i+state_size+j)-(state_size-4)]*S[(i+state_size+j)-(state_size-9 )]+ S[(i+state_size+j)-(state_size-4 )]+ S[(i+state_size+j)-(state_size-1)]*S[(i+state_size+j)-(state_size-13)]*S[(i+state_size+j)-(state_size-9 )]+ S[(i+state_size+j)-(state_size-1)]*S[(i+state_size+j)-(state_size-13 )]+ S[(i+state_size+j)-(state_size-1)]*S[(i+state_size+j)-(state_size-9 )]+ S[(i+state_size+j)-(state_size-1)]+ (S[(i+state_size+j)-(state_size-13)]+S[(i+state_size+j)-(state_size-1)]*S[(i+state_size+j)-(state_size-9)]+S[(i+state_size+j)-(state_size-4)]*S[(i+state_size+j)-(state_size-9)]+S[(i+state_size+j)-(state_size-9)]*S[(i+state_size+j)-(state_size-13)]+1)*O[i+j]#     for j in range(len(O)):
#         m = eval("x" + str(state_size+j))
#         n = eval('O['+str(j)+']')
#         Final_System[i]=Final_System[i].subs({m:n})
#Define linear substitution variables
U = [0]*no_monomials
for i in range(no_monomials):
    U[i] = eval("x" + str(state_size+i))
# print(len(U))
# print(s)
#create monomial list and count the number of combinations
monomials =[1]*1
for i in range(1,len(U)):
    C = Combinations(s,i)
    no = C.cardinality()
    monomials_temp = [1]*no
    for i in range(no):
        for x in C[i]:
            monomials_temp[i]=monomials_temp[i]*x
    monomials = monomials+monomials_temp

#initialise the dictionary
d = {}
#fill the dictionary with the substituion variables
for i in range(1,len(U)+1):
    d[monomials[i]] = U[i-1]

Final_System_SubWithoutY = [0]*len(Final_System)
Final_System_SubWithY = [0]*len(Final_System)
Final_System_Sub = [0]*len(Final_System)

#For each equation in the system
for i in range(len(Final_System)):

    #NOTE:if you sub for the output bits first, you then only need to do the first line of code as there will no longer be output variables
    #make substitution for monomials without output bits
    Final_System_Sub[i] = sum(Final_System[i].monomial_coefficient(m) * v for m,v in d.items())



F = GF(2)
M = matrix(F,len(Final_System_Sub),len(U))
#Fill coefficient matrix with 1s in posisitons corresponding to initial
# state bits and odd keystream bits
for j in range(len(Final_System_Sub)):
    for i in range(len(U)):
        M[j,i] = Final_System_Sub[j].monomial_coefficient(U[i])
w= vector(F,[0]*len(Final_System_Sub))

for i in range(len(Final_System)):
#     if (EQ_1[i].monomial_coefficient(1)==1):
    if (Final_System[i].monomial_coefficient(1)==1):      
        w[i]=1
        
# A = print(M.solve_right(w))
A = M.right_kernel(w)
print(A.dimension())
print(A[1])

stop = timeit.timeit()

print(stop-start)


original_output = sys.stdout
fp = open("Better_n_15.txt","w")
sys.stdout = fp

for i in range(A.dimension()):
    print(A[i])


