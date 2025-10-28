# this code can only be used for the minimal configuration of the bilocality scenario. Note that the cardinality of the probability distribution should be (2, 2, 2, 2, 2) (corresponding to (a, b, c, x, z))
import numpy as np
import gurobipy as gp
from gurobipy import GRB

def testbilocal(p):
    with gp.Env(empty=True) as env:
        env.setParam('OutputFlag', 0)  # Supress output from solver
        env.start()
        # Create a new model
        m = gp.Model("bilinear",env=env)
        
        # Create variables
        alpha = m.addVar(name="alpha")
        beta = m.addVar(name="beta")
        S = np.zeros(16,dtype=object)
        for j in range(16):
            S[j] = m.addVar(name="S"+str(j))
        
        # Set objective: maximize x
        m.setObjective(1, GRB.MAXIMIZE)
        
        #marginal probabilities                   
        p_ac = np.sum(p,axis=1) 
        p_a = np.sum(p_ac,axis=1) 
        p_c = np.sum(p_ac,axis=0) 
        
        #to obtain the rectangle's surface:                    
        x = [0,0,0,0]
        z = [0,0,0,0] 
        
        x[0] = alpha
        x[1] = p_a[1,0,0] - alpha
        x[2] = p_a[1,1,0] - alpha
        x[3] = p_a[0,1,0] - p_a[1,0,0] + alpha
        z[0] = beta
        z[1] = p_c[1,0,0] - beta
        z[2] = p_c[1,0,1] - beta
        z[3] = p_c[0,0,1] - p_c[1,0,0] + beta
        
        R = []
        for i in range(4):
            for k in range(4):
                R=np.append(R,z[i]*x[k])
        
        # Add constraints: 
        m.addConstrs((S[j]>=0 for j in range(16)), "cc1") 
        m.addConstrs((S[j]<=R[j] for j in range(16)), "cc2") 
        
        m.addConstr((S[0]+S[1]+S[4]+S[5]==p[1,1,1,0,0]),"c1")
        m.addConstr((S[2]+S[3]+S[6]+S[7]==p[0,1,1,0,0]),"c2")
        m.addConstr((S[8]+S[9]+S[12]+S[13]==p[1,1,0,0,0]),"c3")
        m.addConstr((S[10]+S[11]+S[14]+S[15]==p[0,1,0,0,0]),"c4")
        m.addConstr((S[0]+S[2]+S[4]+S[6]==p[1,1,1,1,0]),"c5")
        m.addConstr((S[8]+S[10]+S[12]+S[14]==p[1,1,0,1,0]),"c6")
        m.addConstr((S[0]+S[1]+S[8]+S[9]==p[1,1,1,0,1]),"c7")
        m.addConstr((S[2]+S[3]+S[10]+S[11]==p[0,1,1,0,1]),"c8")
        m.addConstr((S[0]+S[2]+S[8]+S[10]==p[1,1,1,1,1]),"c9")
        #m.addConstr((S.sum()<=1),"c99") #this condition is redundant
        m.addConstr((alpha<=p_a[1,1,0]),"ca")
        m.addConstr((alpha<=p_a[1,0,0]),"caa")
        m.addConstr((alpha>=np.max([0,p_a[1,1,0]-p_a[0,0,0]])),"caaa")
        m.addConstr((beta<=p_c[1,0,1]),"cb")
        m.addConstr((beta<=p_c[1,0,0]),"cbb")
        m.addConstr((beta>=np.max([0,p_c[1,0,1]-p_c[0,0,0]])),"cbbb")
        
        
        # Solve bilinear model
        m.Params.NonConvex = 2
        
        try:
            m.optimize()
            m.getAttr('x')  # To trigger an error in case it is not solved, because I made it silent
            return 1
        except:
            return 0
