# this code can only be used for the minimal configuration of the bilocality scenario. Note that the cardinality of the probability distribution should be (2, 2, 2, 2, 2) (corresponding to (a, b, c, x, z))
import numpy as np
import gurobipy as gp
from gurobipy import GRB

def test_oneclassicalsource_BC(p):
    with gp.Env(empty=True) as env:
        env.setParam('OutputFlag', 0)  # Supress output from solver
        env.start()
        with gp.Model("bilinear",env=env) as m:
            Qabc0c1x_ = np.zeros((2,2,2,2,2), dtype=object) 
            for a,b,c0,c1,x in np.ndindex(2, 2, 2, 2, 2):
                Qabc0c1x_[a,b,c0,c1,x] = m.addVar(name="Q_%d_%d_%d_%d_%d" % (a, b, c0, c1, x))
            #Now we define the marginals 
            Qbc0c1x_ = np.array([m.addVar(name="Q_%d_%d_%d_%d" % (b,c0,c1,x)) for b,c0,c1,x in np.ndindex(2, 2, 2, 2)], dtype=object).reshape((2,2,2,2))
            Qac0c1x_ = np.array([m.addVar(name="Q_%d_%d_%d_%d" % (a,c0,c1,x)) for a,c0,c1,x in np.ndindex(2, 2, 2, 2)], dtype=object).reshape((2,2,2,2))
            Qabc1x_ = np.array([m.addVar(name="Q_%d_%d_%d_%d" % (a,b,c1,x)) for a,b,c1,x in np.ndindex(2, 2, 2, 2)], dtype=object).reshape((2,2,2,2))
            Qabc0x_ = np.array([m.addVar(name="Q_%d_%d_%d_%d" % (a,b,c0,x)) for a,b,c0,x in np.ndindex(2, 2, 2, 2)], dtype=object).reshape((2,2,2,2))
            Qc0c1_ = np.array([m.addVar(name="Q_%d_%d" % (c0,c1)) for c0,c1 in np.ndindex(2, 2)], dtype=object).reshape((2,2))
            Qax_ = np.array([m.addVar(name="Q_%d_%d" % (a,x)) for a,x in np.ndindex(2, 2)], dtype=object).reshape((2,2))
            
            #Constraints for positivity and less than 1
            for Qbc0c1x in Qbc0c1x_.flatten():
                m.addConstr(Qbc0c1x >= 0)
                m.addConstr(Qbc0c1x <= 1)
            for Qac0c1x in Qac0c1x_.flatten():
                m.addConstr(Qac0c1x >= 0)
                m.addConstr(Qac0c1x <= 1)
            for Qabc1x in Qabc1x_.flatten():
                m.addConstr(Qabc1x >= 0)
                m.addConstr(Qabc1x <= 1)
            for Qabc0x in Qabc0x_.flatten():
                m.addConstr(Qabc0x >= 0)
                m.addConstr(Qabc0x <= 1)
            for Qc0c1 in Qc0c1_.flatten():
                m.addConstr(Qc0c1 >= 0)
                m.addConstr(Qc0c1 <= 1)
            for Qax in Qax_.flatten():
                m.addConstr(Qax >= 0)
                m.addConstr(Qax <= 1)
            for Qabc0c1x in Qabc0c1x_.flatten():
                m.addConstr(Qabc0c1x >= 0)
                m.addConstr(Qabc0c1x <= 1)
                
            sum_a_Qabc0c1x = Qabc0c1x_.sum(axis=0)
            for b, c0, c1, x in np.ndindex(2, 2, 2, 2):
                m.addConstr(Qbc0c1x_[b,c0,c1,x] == sum_a_Qabc0c1x[b,c0,c1,x])
                
            sum_b_Qabc0c1x = Qabc0c1x_.sum(axis=1)
            for a, c0, c1, x in np.ndindex(2, 2, 2, 2):
                m.addConstr(Qac0c1x_[a,c0,c1,x] == sum_b_Qabc0c1x[a,c0,c1,x])
                
            sum_c0_Qabc0c1x = Qabc0c1x_.sum(axis=2)
            for a, b, c1, x in np.ndindex(2, 2, 2, 2):
                m.addConstr(Qabc1x_[a,b,c1,x] == sum_c0_Qabc0c1x[a,b,c1,x])
            
            sum_c1_Qabc0c1x = Qabc0c1x_.sum(axis=3)
            for a, b, c0, x in np.ndindex(2, 2, 2, 2):
                m.addConstr(Qabc0x_[a,b,c0,x] == sum_c1_Qabc0c1x[a,b,c0,x])
                
            sum_ab_Qabc0c1x = Qabc0c1x_.sum(axis=(0,1))[:,:,0]
            for c0, c1 in np.ndindex(2, 2):
                m.addConstr(Qc0c1_[c0,c1] == sum_ab_Qabc0c1x[c0,c1])
            
            sum_bc0c1_Qabc0c1x = Qabc0c1x_.sum(axis=(1,2,3))
            for a,x in np.ndindex(2, 2):
                m.addConstr(Qax_[a,x] == sum_bc0c1_Qabc0c1x[a,x])
                
            for x in range(2):
                for a, b, c0, c1 in np.ndindex(2, 2, 2, 2):
                    m.addConstr(Qabc0c1x_.sum(axis=(0,1,2,3))[x]==1)
            
            #Constraints NS
            for b,c0,c1 in np.ndindex(2, 2, 2):
                m.addConstr(Qbc0c1x_[b,c0,c1,0]==Qbc0c1x_[b,c0,c1,1])
            
            #Constraints independence
            for a,c0,c1,x in np.ndindex(2, 2, 2, 2):
                m.addConstr(Qac0c1x_[a,c0,c1,x]==Qax_[a,x]*Qc0c1_[c0,c1])
            
            #Third type of constraints
            for a,b,c0,c1,x in np.ndindex(2, 2, 2, 2, 2):
                m.addConstr(Qabc0x_[a,b,c0,x]==p[a,b,c0,x,0])
                m.addConstr(Qabc1x_[a,b,c1,x]==p[a,b,c1,x,1])
             
            # Feasibility problem
            m.setObjective(0.0, GRB.MAXIMIZE)

            # Solve bilinear model
            m.Params.NonConvex = 2

            try:
                m.optimize()
                m.getAttr('x')  # To trigger an error in case it is not solved, because I made it silent
                return 1
            except:
                return 0
