for i in range(5):
    filename = "InFile" + str(i) + ".txt"
    f = open(filename, "w")
    f.write(
"""1 0 0  
5  300 0.3 15000 2000 2000""")

    f.write(" " + str(200 * i) + "\n")
    
    f.write(
"""0.0001    0.1 0.20   0.01
2   1149.00   1149.00
1.0000    0.0000
0.0000    0.0000
0.0000    0.0000
0.0000    0.0000
5    0.05  1000.00000

7  1.0
10300  1 0 0 0
1149  10 0 0 0
508  10 0 0 0
291   10 0 40 40
474   10 0 0 0
843   10 0 0 0
333  10  0 0 0
7   0.5000
10300  1 0 0 0
1149  10 0 0 0
508   10 0 0 0
291   10 0 40 40
474   10 0 0 0
843   10 0 0 0
333  10  0 0 0

0 3 2 2 1 1 1 
1 3 2 2 1 1 1""")