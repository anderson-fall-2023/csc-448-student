import numpy as np
import pandas as pd

from sklearn.linear_model import LinearRegression

def exercise_1():
    a = None
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    a = np.ones((6,4),dtype=int)*2
    ## END SOLUTION
    return a

def exercise_2():
    b = None
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    b = np.ones((6,4))
    b[range(4),range(4)] = 3
    ## END SOLUTION
    return b

def exercise_3(c):
    m = None
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    m = c.mean()
    ## END SOLUTION
    return m

def exercise_4(c):
    row_means,col_means = None,None
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    row_means,col_means = c.mean(axis=1),c.mean(axis=0)
    ## END SOLUTION
    return row_means,col_means

def exercise_5(a):
    c1 = 0
    c2 = 0
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            if a[i,j] == 1:
                c1 += 1
    c2 = (a == 1).sum()
    ## END SOLUTION
    
    return c1,c2

def exercise_6():
    a = None
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    a = pd.DataFrame(np.ones((6,4),dtype=int)*2)
    ## END SOLUTION
    return a

def exercise_7():
    b = None
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    b = np.ones((6,4))
    b[range(4),range(4)] = 3
    b = pd.DataFrame(b)
    ## END SOLUTION
    return b

def exercise_8(df):
    sel = None
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    sel = df['age']
    ## END SOLUTION
    return sel

def exercise_9(df):
    sel = None
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    sel = df.set_index('sex').loc['Male']
    ## END SOLUTION
    return sel

def exercise_10(df):
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    df.reset_index(inplace=True)
    ## END SOLUTION
    return

def exercise_11(X,y):
    # YOUR SOLUTION HERE
    model = None
    ## BEGIN SOLUTION
    model = LinearRegression().fit(X,y)
    ## END SOLUTION
    return model
