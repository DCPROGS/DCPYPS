import dcpyps

def CH82():
    
    mectitle = 'CH82'
    ratetitle = 'CH82 numerical example'

    A2RS = dcpyps.State('A', 'A2R*', 60e-12)
    ARS  = dcpyps.State('A', 'AR*', 60e-12)
    A2R  = dcpyps.State('B', 'A2R', 0.0)
    AR   = dcpyps.State('B', 'AR', 0.0)
    R    = dcpyps.State('C', 'R', 0.0)

    RateList = [
         dcpyps.Rate(15.0, AR, ARS, name='beta1', limits=[1e-15,1e+7]),
         dcpyps.Rate(15000.0, A2R, A2RS, name='beta2', limits=[1e-15,1e+7]),
         dcpyps.Rate(3000.0, ARS, AR, name='alpha1', limits=[1e-15,1e+7]),
         dcpyps.Rate(500.0, A2RS, A2R, name='alpha2', limits=[1e-15,1e+7]),
         dcpyps.Rate(2000.0, AR, R, name='k(-1)', limits=[1e-15,1e+7]),
         dcpyps.Rate(2 * 2000.0, A2R, AR, name='2k(-2)', limits=[1e-15,1e+7]),
         dcpyps.Rate(2 * 5.0e07, R, AR, name='2k(+1)', eff='c', limits=[1e-15,1e+10]),
         dcpyps.Rate(5.0e08, ARS, A2RS, name='k*(+2)', eff='c', fixed=True, limits=[1e-15,1e+10]),
         dcpyps.Rate(5.0e08, AR, A2R, name='k(+2)', eff='c', limits=[1e-15,1e+10]),
         #dcpyps.Rate(2 * 1.0 / 3.0, A2RS, ARS, name='k*(-2)', limits=[1e-15,1e+7])
         dcpyps.Rate(0.66667, A2RS, ARS, name='2k*(-2)', mr=True, limits=[1e-15,1e+7])
         ]

    CycleList = [dcpyps.Cycle(['A2R*', 'AR*', 'AR', 'A2R'], ['A2R*', 'AR*'])]

    fastblk = False
    KBlk = 0.001

    return  dcpyps.Mechanism(RateList, CycleList, mtitle=mectitle, rtitle=ratetitle) #, fastblk, KBlk)

def AChR_diamond():
    
    mectitle = 'diamond'
    ratetitle = 'from CHH 2003'

    A2RS = dcpyps.State('A', 'A2R*', 60e-12)
    ARSa  = dcpyps.State('A', 'AR*a', 60e-12)
    ARSb  = dcpyps.State('A', 'AR*b', 60e-12)
    A2R  = dcpyps.State('B', 'A2R', 0.0)
    ARa   = dcpyps.State('B', 'ARa', 0.0)
    ARb   = dcpyps.State('B', 'ARb', 0.0)
    R    = dcpyps.State('C', 'R', 0.0)

    RateList = [
         dcpyps.Rate(50.0, ARa, ARSa, name='beta1a', limits=[1e-2,1e+6]),
         dcpyps.Rate(6000.0, ARSa, ARa, name='alpha1a', limits=[1e-2,1e+6]),
         dcpyps.Rate(150.0, ARb, ARSb, name='beta1b', limits=[1e-2,1e+6]),
         dcpyps.Rate(50000.0, ARSb, ARb, name='alpha1b', limits=[1e-2,1e+6]),
         dcpyps.Rate(52000.0, A2R, A2RS, name='beta2', limits=[1e-2,1e+6]),
         dcpyps.Rate(2000.0, A2RS, A2R, name='alpha2', limits=[1e-2,1e+6]),
         
         dcpyps.Rate(1500.0, A2R, ARb, name='k(-2a)', limits=[1e-2,1e+6]),
         dcpyps.Rate(2.0e08, ARb, A2R, name='k(+2a)', eff='c', limits=[1e-2,1e+10]),
         dcpyps.Rate(10000.0, A2R, ARa, name='k(-2b)', limits=[1e-2,1e+6]),
         dcpyps.Rate(4.0e08, ARa, A2R, name='k(+2b)', eff='c', limits=[1e-2,1e+10]),
         
         dcpyps.Rate(1500.0, ARa, R, name='k(-1a)', limits=[1e-2,1e+6]),
         dcpyps.Rate(2.0e08, R, ARa, name='k(+1a)', eff='c', limits=[1e-2,1e+10]),
         dcpyps.Rate(10000.0, ARb, R, name='k(-1b)', limits=[1e-2,1e+6]),
         dcpyps.Rate(4.0e08, R, ARb, name='k(+1b)', eff='c', limits=[1e-2,1e+10])         
         ]

    CycleList = [dcpyps.Cycle(['A2R', 'ARa', 'R', 'ARb'])]

    return  dcpyps.Mechanism(RateList, CycleList, mtitle=mectitle, rtitle=ratetitle)

def load_AChR_diamond_independent_binding(rates=None):
    mec = AChR_diamond()
    # PREPARE RATE CONSTANTS.
    if rates is not None:
        mec.set_rateconstants(rates)
    mec.Rates[11].is_constrained = True
    mec.Rates[11].constrain_func = dcpyps.mechanism.constrain_rate_multiple
    mec.Rates[11].constrain_args = [7, 1]             
    mec.Rates[13].is_constrained = True
    mec.Rates[13].constrain_func = dcpyps.mechanism.constrain_rate_multiple
    mec.Rates[13].constrain_args = [9, 1]
    mec.Rates[10].is_constrained = True
    mec.Rates[10].constrain_func = dcpyps.mechanism.constrain_rate_multiple
    mec.Rates[10].constrain_args = [6, 1]
    mec.Rates[12].is_constrained = True
    mec.Rates[12].constrain_func = dcpyps.mechanism.constrain_rate_multiple
    mec.Rates[12].constrain_args = [8, 1]
    mec.update_constrains()
    return mec, mec.get_free_parameter_names()

def load_GlyR_flip_independent_binding(rates=None):
    # LOAD FLIP MECHANISM USED Burzomato et al 2004
    #mecfn = "/mec/demomec.mec"
    #version, meclist, max_mecnum = dcpyps.dcio.mec_get_list(mecfn)
    #mec = dcpyps.dcio.mec_load(mecfn, meclist[2][0])
    mec = GlyR_flip()
    # PREPARE RATE CONSTANTS.
    if rates is not None:
        mec.set_rateconstants(rates)
    # Fixed rates.
    for i in range(len(mec.Rates)):
        mec.Rates[i].fixed = False
    # Constrained rates.
    mec.Rates[21].is_constrained = True
    mec.Rates[21].constrain_func = dcpyps.mechanism.constrain_rate_multiple
    mec.Rates[21].constrain_args = [17, 3]
    mec.Rates[19].is_constrained = True
    mec.Rates[19].constrain_func = dcpyps.mechanism.constrain_rate_multiple
    mec.Rates[19].constrain_args = [17, 2]
    mec.Rates[16].is_constrained = True
    mec.Rates[16].constrain_func = dcpyps.mechanism.constrain_rate_multiple
    mec.Rates[16].constrain_args = [20, 3]
    mec.Rates[18].is_constrained = True
    mec.Rates[18].constrain_func = dcpyps.mechanism.constrain_rate_multiple
    mec.Rates[18].constrain_args = [20, 2]
    mec.Rates[8].is_constrained = True
    mec.Rates[8].constrain_func = dcpyps.mechanism.constrain_rate_multiple
    mec.Rates[8].constrain_args = [12, 1.5]
    mec.Rates[13].is_constrained = True
    mec.Rates[13].constrain_func = dcpyps.mechanism.constrain_rate_multiple
    mec.Rates[13].constrain_args = [9, 2]
    mec.update_constrains()
    mec.set_mr(True, 7, 1)
    mec.set_mr(True, 15, 0)
    mec.update_constrains()
    return mec, mec.get_free_parameter_names()


def GlyR_flip():
    
    mectitle = '3 site flip'
    ratetitle = 'Burzomato2004 heteromericGlyR'

    A3FS = dcpyps.State('A', 'A3F*', 60e-12)
    A2FS = dcpyps.State('A', 'A2F*', 60e-12)
    AFS  = dcpyps.State('A', 'AF*', 60e-12)
    A3F  = dcpyps.State('B', 'A3F', 0.0)
    A2F  = dcpyps.State('B', 'A2F', 0.0)
    AF   = dcpyps.State('B', 'AF', 0.0)
    A3R  = dcpyps.State('B', 'A3R', 0.0)
    A2R  = dcpyps.State('B', 'A2R', 0.0)
    AR   = dcpyps.State('B', 'AR', 0.0)
    R    = dcpyps.State('C', 'R', 0.0)

    RateList = [
        dcpyps.Rate(3400.0, AFS, AF, name='alpha1', limits=[1e-15,1e+7]),
        dcpyps.Rate(4200.0, AF, AFS, name='beta1', limits=[1e-15,1e+7]),
        dcpyps.Rate(2100.0, A2FS, A2F, name='alpha2', limits=[1e-15,1e+7]),
        dcpyps.Rate(28000.0, A2F, A2FS, name='beta2', limits=[1e-15,1e+7]),
        dcpyps.Rate(7000.0, A3FS, A3F, name='alpha3', limits=[1e-15,1e+7]),
        dcpyps.Rate(129000.0, A3F, A3FS, name='beta3', limits=[1e-15,1e+7]),
        dcpyps.Rate(900.0, A3F, A3R, name='gamma3', limits=[1e-15,1e+7]),
        dcpyps.Rate(20900.0, A3R, A3F, name='delta3', limits=[1e-15,1e+7]),
        dcpyps.Rate(3 * 1200, A3F, A2F, name='3kf(-3)', limits=[1e-15,1e+7]),
        dcpyps.Rate(150e06, A2F, A3F, name='kf(+3)', eff='c', limits=[1e-15,1e+10]),
        dcpyps.Rate(18000.0, A2F, A2R, name='gamma2', limits=[1e-15,1e+7]),
        dcpyps.Rate(6800.0, A2R, A2F, name='delta2', limits=[1e-15,1e+7]),
        dcpyps.Rate(2 * 1200, A2F, AF, name='2kf(-2)', limits=[1e-15,1e+7]),
        dcpyps.Rate(2 * 150e06, AF, A2F, name='2kf(+2)', eff='c', limits=[1e-15,1e+10]),
        dcpyps.Rate(29000.0, AF, AR, name='gamma1', limits=[1e-15,1e+7]),
        dcpyps.Rate(180.0, AR, AF, name='delta1', limits=[1e-15,1e+7]),
        dcpyps.Rate(3 * 300.0, A3R, A2R, name='3k(-3)', limits=[1e-15,1e+7]),
        dcpyps.Rate(0.59e06, A2R, A3R, name='k(+3)', eff='c', limits=[1e-15,1e+10]),
        dcpyps.Rate(2 * 300.0, A2R, AR, name='2k(-2)', limits=[1e-15,1e+7]),
        dcpyps.Rate(2 * 0.59e06, AR, A2R, name='2k(+2)', eff='c', limits=[1e-15,1e+10]),
        dcpyps.Rate(300.0, AR, R, name='k(-1)', limits=[1e-15,1e+7]),
        dcpyps.Rate(3 * 0.59e06, R, AR, name='3k(+1)', eff='c', limits=[1e-15,1e+10])
        ]

    CycleList = [dcpyps.Cycle(['A2F', 'AF', 'AR', 'A2R']), 
        dcpyps.Cycle(['A3F', 'A2F', 'A2R', 'A3R'])]
    return  dcpyps.Mechanism(RateList, CycleList, mtitle=mectitle, rtitle=ratetitle)

def fully_connected_cycle():
    
    mectitle = 'Fully connected cycle'
    ratetitle = 'JP numbers'

    O2 = dcpyps.State('A', 'O2', 60e-12)
    O1 = dcpyps.State('A', 'O1', 60e-12)
    C2 = dcpyps.State('B', 'C2', 0.0)
    C1 = dcpyps.State('C', 'C1', 0.0)
    
    KC, tKO, KO = 1.0, 0.05, 5.0
    k12, k13, k14 = 1.2, 1.3, 1.4
    k23, k24, k34 = 2.3, 2.4, 3.4

    RateList = [
         dcpyps.Rate(k12,      C1, C2, name='k12', eff='c', limits=[1e-15,1e+7]),
         dcpyps.Rate(k12 / KC, C2, C1, name='k21', limits=[1e-15,1e+7]),
         dcpyps.Rate(k24 / KC, C2, O2, name='k24', limits=[1e-15,1e+7]),
         dcpyps.Rate(k24 / KO, O2, C2, name='k42', limits=[1e-15,1e+7]),
         dcpyps.Rate(k34 / KO, O2, O1, name='k34', limits=[1e-15,1e+7]),
         dcpyps.Rate(k34 /tKO, O1, O2, name='k43', eff='c', limits=[1e-15,1e+7]),
         dcpyps.Rate(k13 /tKO, O1, C1, name='k13', limits=[1e-15,1e+10]),
         dcpyps.Rate(k13,      C1, O1, name='k31', limits=[1e-15,1e+10]),
         dcpyps.Rate(k14,      C1, O2, name='k14', eff='c', limits=[1e-15,1e+10]),
         dcpyps.Rate(k14 / KO, O2, C1, name='k41', limits=[1e-15,1e+7]),
         dcpyps.Rate(k23 / KC, C2, O1, name='k23', limits=[1e-15,1e+7]),
         dcpyps.Rate(k23 /tKO, O1, C2, name='k32', eff='c', limits=[1e-15,1e+10]),
         ]

    CycleList = [
        dcpyps.Cycle(['C1', 'C2', 'O2'], ['C2', 'O2']),
        dcpyps.Cycle(['C2', 'O2', 'O1'], ['O2', 'O1']),
        dcpyps.Cycle(['C1', 'C2', 'O1'], ['C1', 'O1']),
        ]

    return  dcpyps.Mechanism(RateList, CycleList, mtitle=mectitle, rtitle=ratetitle)

def CCO():
    
    mectitle = 'C-C-O'
    ratetitle = 'quasi random numbers'

    ARS  = dcpyps.State('A', 'AR*', 50e-12)
    AR   = dcpyps.State('B', 'AR', 0.0)
    R    = dcpyps.State('C', 'R', 0.0)

    RateList = [
         dcpyps.Rate(15000.0, AR, ARS, name='beta', limits=[1e-15,1e+7]),
         dcpyps.Rate(500.0, ARS, AR, name='alpha', limits=[1e-15,1e+7]),
         dcpyps.Rate(2000.0, AR, R, name='koff', limits=[1e-15,1e+7]),
         dcpyps.Rate(5.0e08, R, AR, name='kon', eff='c', limits=[1e-15,1e+10]),
         ]

    return  dcpyps.Mechanism(RateList, mtitle=mectitle, rtitle=ratetitle)

def CCOD():
    
    mectitle = 'C-C-O-D'
    ratetitle = 'quasi random numbers'

    AD   = dcpyps.State('B', 'AD', 0.0)
    ARS  = dcpyps.State('A', 'AR*', 50e-12)
    AR   = dcpyps.State('B', 'AR', 0.0)
    R    = dcpyps.State('C', 'R', 0.0)

    RateList = [
         dcpyps.Rate(15000.0, AR, ARS, name='beta', limits=[1e-15,1e+7]),
         dcpyps.Rate(500.0, ARS, AR, name='alpha', limits=[1e-15,1e+7]),
         dcpyps.Rate(200.0, AD, ARS, name='doff', limits=[1e-15,1e+7]),
         dcpyps.Rate(10.0, ARS, AD, name='don', limits=[1e-15,1e+7]),
         dcpyps.Rate(2000.0, AR, R, name='koff', limits=[1e-15,1e+7]),
         dcpyps.Rate(5.0e08, R, AR, name='kon', eff='c', limits=[1e-15,1e+10]),
         ]

    return  dcpyps.Mechanism(RateList, mtitle=mectitle, rtitle=ratetitle)

def CCOB():
    
    mectitle = 'C-C-O-B'
    ratetitle = 'quasi random numbers'

    ARSB = dcpyps.State('B', 'AR*B', 0.0)
    ARS  = dcpyps.State('A', 'AR*', 50e-12)
    AR   = dcpyps.State('B', 'AR', 0.0)
    R    = dcpyps.State('C', 'R', 0.0)

    RateList = [
         dcpyps.Rate(15000.0, AR, ARS, name='beta', limits=[1e-15,1e+7]),
         dcpyps.Rate(500.0, ARS, AR, name='alpha', limits=[1e-15,1e+7]),
         dcpyps.Rate(2000.0, AR, R, name='koff', limits=[1e-15,1e+7]),
         dcpyps.Rate(5.0e08, R, AR, name='kon', eff='c', limits=[1e-15,1e+10]),
         dcpyps.Rate(9.0e07, ARS, ARSB, name='kBon', eff='c', limits=[1e-15,1e+10]),
         dcpyps.Rate(90000.0, ARSB, ARS, name='kBoff', limits=[1e-15,1e+7])
         ]

    return  dcpyps.Mechanism(RateList, mtitle=mectitle, rtitle=ratetitle)

def CO():
    
    mectitle = 'C-O'
    ratetitle = 'quasi random numbers'

    RS  = dcpyps.State('A', 'O', 50e-12)
    R   = dcpyps.State('B', 'C', 0.0)

    RateList = [
         dcpyps.Rate(20.0, R, RS, name='beta', limits=[1e-15,1e+7]),
         dcpyps.Rate(50.0, RS, R, name='alpha', limits=[1e-15,1e+7]),
         ]

    return  dcpyps.Mechanism(RateList, mtitle=mectitle, rtitle=ratetitle)

def six_cycles_mec():
    
    mectitle = 'six cycles'
    ratetitle = 'quasi random numbers'

    A = dcpyps.State('A', 'A', 60e-12)
    B = dcpyps.State('A', 'B', 60e-12)
    C = dcpyps.State('A', 'C', 60e-12)
    D = dcpyps.State('B', 'D', 60e-12)
    E = dcpyps.State('B', 'E', 0.0)
    F = dcpyps.State('B', 'F', 0.0)
    G = dcpyps.State('B', 'G', 0.0)
    H = dcpyps.State('B', 'H', 0.0)
    I = dcpyps.State('B', 'I', 0.0)
    J = dcpyps.State('B', 'J', 0.0)
    K = dcpyps.State('B', 'K', 0.0)
    L = dcpyps.State('C', 'L', 0.0)

    RateList = [
         dcpyps.Rate(15.0, A, B, name='ab'),
         dcpyps.Rate(25.0, B, A, name='ba'),
         dcpyps.Rate(15.0, A, E, name='ae'),
         dcpyps.Rate(25.0, E, A, name='ea'),
         dcpyps.Rate(15.0, B, C, name='bc'),
         dcpyps.Rate(25.0, C, B, name='cb'),
         dcpyps.Rate(15.0, B, F, name='bf'),
         dcpyps.Rate(25.0, F, B, name='fb'),
         dcpyps.Rate(15.0, C, D, name='cd'),
         dcpyps.Rate(25.0, D, C, name='dc'),
         dcpyps.Rate(15.0, C, G, name='cg'),
         dcpyps.Rate(25.0, G, C, name='gc'),
         dcpyps.Rate(15.0, D, H, name='dh'),
         dcpyps.Rate(25.0, H, D, name='hd'),
         dcpyps.Rate(15.0, E, F, name='ef'),
         dcpyps.Rate(25.0, F, E, name='fe'),
         dcpyps.Rate(15.0, E, I, name='ei'),
         dcpyps.Rate(25.0, I, E, name='ie'),
         dcpyps.Rate(15.0, F, G, name='fg'),
         dcpyps.Rate(25.0, G, F, name='gf'),
         dcpyps.Rate(15.0, F, J, name='fj'),
         dcpyps.Rate(25.0, J, F, name='jf'),
         dcpyps.Rate(15.0, G, H, name='gh'),
         dcpyps.Rate(25.0, H, G, name='hg'),
         dcpyps.Rate(15.0, G, K, name='gk'),
         dcpyps.Rate(25.0, K, G, name='kg'),
         dcpyps.Rate(15.0, H, L, name='hl'),
         dcpyps.Rate(25.0, L, H, name='lh'),
         dcpyps.Rate(15.0, I, J, name='ij'),
         dcpyps.Rate(25.0, J, I, name='ji'),
         dcpyps.Rate(15.0, J, K, name='jk'),
         dcpyps.Rate(25.0, K, J, name='kj'),
         dcpyps.Rate(15.0, K, L, name='kl'),
         dcpyps.Rate(25.0, L, K, name='lk')
        ]

    return  dcpyps.Mechanism(RateList, mtitle=mectitle, rtitle=ratetitle)
