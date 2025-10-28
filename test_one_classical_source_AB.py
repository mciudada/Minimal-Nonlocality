import numpy as np
import gurobipy as gp
from gurobipy import GRB

def test_oneclassicalsource_AB(p):
    with gp.Env(empty=True) as env:
        env.setParam('OutputFlag', 0)  # Supress output from solver
        env.start()
        with gp.Model("bilinear",env=env) as m:
            Qa0a1bcz_ = np.zeros((2,2,2,2,2), dtype=object) 
            for a0,a1,b,c,z in np.ndindex(2, 2, 2, 2, 2):
                Qa0a1bcz_[a0,a1,b,c,z] = m.addVar(name="Q_%d_%d_%d_%d_%d" % (a0, a1, b, c, z))
            #Now we define the marginals 
            Qa0a1bz_ = np.array([m.addVar(name="Q_%d_%d_%d_%d" % (a0,a1,b,z)) for a0,a1,b,z in np.ndindex(2, 2, 2, 2)], dtype=object).reshape((2,2,2,2))
            Qa0a1cz_ = np.array([m.addVar(name="Q_%d_%d_%d_%d" % (a0,a1,c,z)) for a0,a1,c,z in np.ndindex(2, 2, 2, 2)], dtype=object).reshape((2,2,2,2))
            Qa0bcz_ = np.array([m.addVar(name="Q_%d_%d_%d_%d" % (a0,b,c,z)) for a0,b,c,z in np.ndindex(2, 2, 2, 2)], dtype=object).reshape((2,2,2,2))
            Qa1bcz_ = np.array([m.addVar(name="Q_%d_%d_%d_%d" % (a1,b,c,z)) for a1,b,c,z in np.ndindex(2, 2, 2, 2)], dtype=object).reshape((2,2,2,2))
            Qa0a1_ = np.array([m.addVar(name="Q_%d_%d" % (a0,a1)) for a0,a1 in np.ndindex(2, 2)], dtype=object).reshape((2,2))
            Qcz_ = np.array([m.addVar(name="Q_%d_%d" % (c,z)) for c,z in np.ndindex(2, 2)], dtype=object).reshape((2,2))
            
            #Constraints for positivity and less than 1
            for Qa0a1bz in Qa0a1bz_.flatten():
                m.addConstr(Qa0a1bz >= 0)
                m.addConstr(Qa0a1bz <= 1)
            for Qa0a1cz in Qa0a1cz_.flatten():
                m.addConstr(Qa0a1cz >= 0)
                m.addConstr(Qa0a1cz <= 1)
            for Qa0bcz in Qa0bcz_.flatten():
                m.addConstr(Qa0bcz >= 0)
                m.addConstr(Qa0bcz <= 1)
            for Qa1bcz in Qa1bcz_.flatten():
                m.addConstr(Qa1bcz >= 0)
                m.addConstr(Qa1bcz <= 1)
            for Qa0a1 in Qa0a1_.flatten():
                m.addConstr(Qa0a1 >= 0)
                m.addConstr(Qa0a1 <= 1)
            for Qcz in Qcz_.flatten():
                m.addConstr(Qcz >= 0)
                m.addConstr(Qcz <= 1)
            for Qa0a1bcz in Qa0a1bcz_.flatten():
                m.addConstr(Qa0a1bcz >= 0)
                m.addConstr(Qa0a1bcz <= 1)
                
            sum_c_Qa0a1bcz = Qa0a1bcz_.sum(axis=3)
            for a0, a1, b, z in np.ndindex(2, 2, 2, 2):
                m.addConstr(Qa0a1bz_[a0,a1,b,z] == sum_c_Qa0a1bcz[a0,a1,b,z])
                
            sum_b_Qa0a1bcz = Qa0a1bcz_.sum(axis=2)
            for a0, a1, c, z in np.ndindex(2, 2, 2, 2):
                m.addConstr(Qa0a1cz_[a0,a1,c,z]==sum_b_Qa0a1bcz[a0,a1,c,z])
            
            sum_a1_Qa0a1bcz = Qa0a1bcz_.sum(axis=1)
            for a0, b, c, z in np.ndindex(2, 2, 2, 2):
                m.addConstr(Qa0bcz_[a0,b,c,z]==sum_a1_Qa0a1bcz[a0,b,c,z])
            
            sum_a0_Qa0a1bcz = Qa0a1bcz_.sum(axis=0)
            for a1, b, c, z in np.ndindex(2, 2, 2, 2):
                m.addConstr(Qa1bcz_[a1,b,c,z]==sum_a0_Qa0a1bcz[a1,b,c,z])
            
            sum_bc_Qa0a1bcz = Qa0a1bcz_.sum(axis=(2,3))[:,:,0]  
            for a0, a1 in np.ndindex(2, 2):
                m.addConstr(Qa0a1_[a0,a1]==sum_bc_Qa0a1bcz[a0,a1])
            
            sum_a0a1b_Qa0a1bcz = Qa0a1bcz_.sum(axis=(0,1,2))
            for c,z in np.ndindex(2, 2):
                m.addConstr(Qcz_[c,z]==sum_a0a1b_Qa0a1bcz[c,z])
            
            for z in range(2):
                for a0, a1, b, c in np.ndindex(2, 2, 2, 2):
                    m.addConstr(Qa0a1bcz_.sum(axis=(0,1,2,3))[z]==1)
            
            #Constraints NS
            for a0, a1, b in np.ndindex(2, 2, 2):
                m.addConstr(Qa0a1bz_[a0,a1,b,0]==Qa0a1bz_[a0,a1,b,1])
            
            #Constraints independence
            for a0, a1, c, z in np.ndindex(2, 2, 2, 2):
                m.addConstr(Qa0a1cz_[a0,a1,c,z]==Qa0a1_[a0,a1]*Qcz_[c,z])
            
            #Third type of constraints
            for a0, a1, b, c, z in np.ndindex(2, 2, 2, 2, 2):
                m.addConstr(Qa0bcz_[a0,b,c,z]==p[a0,b,c,0,z])
                m.addConstr(Qa1bcz_[a1,b,c,z]==p[a1,b,c,1,z])
             
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
            
