# this code can be used for any bilocal scenario. Note that if any of the parties has no input choice, the cardinality should be one.
import numpy as np
import itertools
import gurobipy as gp
from gurobipy import GRB

def det_strats_1party(nr_outs, nr_ins):
    '''
    Returns all deterministic strategies for a single party as an array
    called as det_strats[output,input,strategy_index]
    '''
    nrdets = nr_outs**nr_ins
    det_strats = np.zeros((nr_outs, nr_ins, nrdets))
    for det_idx in range(nrdets):
        lambda_string = np.unravel_index(det_idx, [nr_outs]*nr_ins)
        for x, ax in enumerate(lambda_string):
            det_strats[ax, x, det_idx] = 1
    return det_strats

def testbilocal_2(prob):    
    dims = prob.shape
    nrparties = int(len(dims)/2)
    outs = dims[:nrparties]
    ins  = dims[nrparties:]
    nrdetstrats = np.power(outs,ins)

    with gp.Env(empty=True) as env:
        env.setParam('OutputFlag', 0)  # Supress output from solver
        env.start()
        with gp.Model("bilinear",env=env) as m:
            Si_ = np.array([m.addVar(name="Si_%d" % i) for i in range(nrdetstrats[0])], dtype=object)
            Sk_ = np.array([m.addVar(name="Sk_%d" % i) for i in range(nrdetstrats[2])], dtype=object)
            Sijk_ = np.zeros(nrdetstrats, dtype=object)
            for alpha, beta, gamma in itertools.product(*[range(x) for x in nrdetstrats]):
                Sijk_[alpha, beta, gamma] = m.addVar(name="Sijk_%d_%d_%d" % (alpha, beta, gamma))

            for Si in Si_:
                m.addConstr(Si >= 0)
                m.addConstr(Si <= 1)
            for Sk in Sk_:
                m.addConstr(Sk >= 0)
                m.addConstr(Sk <= 1)
            for Sijk in Sijk_.flatten():
                m.addConstr(Sijk >= 0)
                m.addConstr(Sijk <= 1)

            sum_jk_Sijk = Sijk_.sum(axis=(1,2))
            for i, Si in enumerate(Si_):
                m.addConstr(sum_jk_Sijk[i] == Si)

            sum_ij_Sijk = Sijk_.sum(axis=(0,1))
            for k, Sk in enumerate(Sk_):
                m.addConstr(sum_ij_Sijk[k] == Sk)

            m.addConstr(Sijk_.sum() == 1, "Sijk_sum")
            
            sum_j_Sijk  = Sijk_.sum(axis=1)
            Sik_ = np.zeros((nrdetstrats[0],nrdetstrats[2]), dtype=object)
            for i, k in itertools.product(*[range(x) for x in (nrdetstrats[0],nrdetstrats[2])]):
                Sik_[i, k] = m.addVar(name="Sik_%d_%d" % (alpha, gamma))
                m.addConstr(Sik_[i, k] == sum_j_Sijk[i, k])

            for i, Si in enumerate(Si_):
                for k, Sk in enumerate(Sk_):
                    m.addConstr(Sik_[i, k] == Si * Sk)

            D_A = det_strats_1party(outs[0], ins[0])
            D_B = det_strats_1party(outs[1], ins[1])
            D_C = det_strats_1party(outs[2], ins[2])
            for a,b,c,x,y,z in itertools.product(*[range(i) for i in (*outs,*ins)]):
                summ = 0
                for alpha, beta, gamma in itertools.product(*[range(x) for x in nrdetstrats]):
                    summ = summ + ( Sijk_[alpha, beta, gamma] * D_A[a, x, alpha]
                                                              * D_B[b, y, beta]
                                                              * D_C[c, z, gamma] )
                m.addConstr(summ == prob[a,b,c,x,y,z])

            m.setObjective(0.0, GRB.MAXIMIZE)
            m.Params.NonConvex = 2
            try:
                m.optimize()
                m.getAttr('x')  # To trigger an error in case it is not solved, because I made it silent
                return 1
            elif m.Status == GRB.INFEASIBLE:
                return 0
            else:
                raise RuntimeError(f"Gurobi stopped with status {m.Status}")
