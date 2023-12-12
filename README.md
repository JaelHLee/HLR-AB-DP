# HLR-AB-DP
Heterogeneous local reconfiguration adaptive behaviors to defend percolation in networks.
"""
The code uses the libraries `networkx` and `seaborn`, which can be installed by running the following instractions in the terminal:
``` bash
pip install networkx
```
``` bash
pip install seaborn
```
"""
#5 Scenario：
Original scenario, network not responding.

LRAB1：Even distribution. reference paper: Rapisardi, Giacomo, Ivan Kryven, and Alex Arenas. "Percolation in networks with local homeostatic plasticity." Nature Communications 13.1 (2022): 122.

LRAB2：Edge weight distribution. 

LRAB3：Source node total out weight flow distribution. 

LRAB4：Maximum edge distribution.

LRAB5：Maximum entropy distribution. 
"""

"""
#3 Indicators：
Code file:
    Indicator 1 code file: 'meank.py'
    
    Indicator 2 code file: 'connectivity.py'
    
    Indicator 3 code file: 'UD.py'
"""
Enter a description of the parameter:
    ins：Data file path.
        For example：ins='datasets/data/***.txt'
        
    dist： weight distribution 'gauss', a theoretical network with edge weights following an N(0.5, 0.12) distribution.
        For example：dist='gauss'
        
    outstring：The main result output document for graphing.
        For example：outstring='out_connectivity'
        
    plot：If or not output image. Optional sets:{'noplot','plot'}, corresponding to: do not output image, output image respectively
       For example：plot='plot'
       
    scenario：Scenario settings. Optional set: a subset of {1, 2, 3, 4, 5}.
        For example：scenario=[1,2,3,4,5]
        
    type:Network type. Optional sets: {'ba','sw','random'}.
        For example：type='ba'
        
    yy:Threshold range. yy = np.linspace(min=0.02,max,num), min in parentheses is the minimum value of 0.02 (i.e., step size), max is the maximum value, and num is the number. Ensure that the step size of 0.02, according to the trend of the line to adjust the maximum value, so that the line tends to stabilize. Step = max/num 
        Example: yy = np.linspace(0.02,1,50) # the minimum value is 0.02, the maximum value is 1, the step size of 0.02
"""
