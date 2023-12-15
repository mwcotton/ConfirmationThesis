import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

# Labels for plotting axis
labels = {'del_mu': r'$\beta\Delta \mu$', 'del_e': '$\Delta e$', 'alpha': r'$\alpha$', 'enz_star': '$\phi_e$', 'Mep': '$M_{ep}$', 'Mes': '$M_{se}$', 'k_rat': '$k_{rat}$', 'd': '$d$',
          'Mee': '$M_{ee}$', 'Mpp': '$M_{pp}$', 'Mss': '$M_{ss}$', 'v': '$v_s$', 'vv': '$v_e$', 'Dep': '$D_{ep}$', 'Dpe': '$D_{pe}$', 'Des': '$D_{se}$', 'Dse': '$D_{se}$', 'Dee': '$D_{ee}$', 'Dpp': '$D_{pp}$', 'Dss': '$D_{ss}$', 'v_rat': r'$\frac{v_e}{v_s}$', 'k_spo': '$k_{spo}$', 'k_cat': '$k_{cat}$'}


def get_steady_bi(del_mu=1, del_e=1, enz_star=0.1, k_spo=1, k_cat=1, Dpe=1, Dse=1, Dps=10, v_rat=10):

    # expect v_rat > 1

    alpha = 1 - enz_star

    RP = Rp(del_mu, del_e, enz_star, k_spo, k_cat)
    RS = Rs(del_e, enz_star, k_spo, k_cat)
    RE = Re(del_mu, del_e, enz_star, 1-enz_star, k_spo, k_cat)

    sub_star = alpha*RP/(RS+RP)
    prod_star = alpha*RS/(RS+RP)

    Mpe = -Dpe*enz_star*prod_star
    Mse = -Dse*enz_star*sub_star
    Mps = -Dps*prod_star*sub_star

    Mes = Mse*v_rat
    Mep = Mpe*v_rat
    Msp = Mps

    Mss = -(Mes+Mps)
    Mpp = -(Mep+Msp)
    Mee = -(Mse+Mpe)

    return sub_star, prod_star, RS, RP, RE, Mee, Mes, Mep, Mse, Mss, Msp, Mpe, Mps, Mpp


def steady_vol(del_mu, del_e, enz_star, alpha, k_spo, k_cat):

    RP = Rp(del_mu, del_e, enz_star, k_spo, k_cat)
    RS = Rs(del_e, enz_star, k_spo, k_cat)

    sub_star = alpha*RP/(RS+RP)
    prod_star = alpha*RS/(RS+RP)

    return sub_star, prod_star


def mu_eff(del_mu=1, del_e=1, enz_star=0.1, k_spo=1, k_cat=1, Dpe=1, Dse=1, v_rat=10):

    k_rat = k_cat/k_spo
    # v_rat = v_e/v_s

    sub_star, prod_star = steady_vol(
        del_mu, del_e, enz_star, 1-enz_star, k_spo, k_cat)

    return np.log(enz_star) - v_rat*np.log(Dse*sub_star+Dpe*prod_star)


def f_eff(del_mu=1, del_e=1, enz_star=0.1, k_spo=1, k_cat=1, Dpe=1, Dse=1, v_rat=10):

    k_rat = k_cat/k_spo
    ink_rat = k_spo/k_cat
    # v_rat = v_e/v_s

    f = (v_rat - 1)*enz_star
    f += v_rat*np.log(1-enz_star)
    f += v_rat*ink_rat*(1+np.exp(del_e))*(np.log(k_rat*enz_star+np.exp(del_e+del_mu)
                                               * (1+np.exp(del_e)+k_rat*enz_star)))/(1+np.exp(-del_e-del_mu))
    f -= v_rat*ink_rat*(Dse+Dpe*np.exp(del_e))*(np.log(Dse*k_rat*enz_star+np.exp(del_e+del_mu)
                                                     * (Dse+Dpe*np.exp(del_e)+Dpe*k_rat*enz_star)))/(Dpe+Dse*np.exp(-del_e-del_mu))
    f += enz_star*mu_eff(del_mu, del_e, enz_star, k_spo,
                         k_cat, Dpe, Dse, v_rat)

    return f


def get_steady_te(del_mu=1, del_e=1, enz_star=0.1, k_spo=1, k_cat=1, Dpe=1, Dse=1, Dps=10, v_e=10, v_s=1, v_w=0.5, alpha=0.1, Dew=10, Dsw=10, Dpw=10):

    wat_star = 1-enz_star-alpha

    RP = Rp(del_mu, del_e, enz_star, k_spo, k_cat)
    RS = Rs(del_e, enz_star, k_spo, k_cat)
    RE = Re(del_mu, del_e, enz_star, alpha, k_spo, k_cat)

    sub_star = alpha*RP/(RS+RP)
    prod_star = alpha*RS/(RS+RP)

    Mpe = -Dpe*enz_star*prod_star
    Mse = -Dse*enz_star*sub_star
    Mps = -Dps*prod_star*sub_star

    Mes = Mse*v_e/v_s
    Mep = Mpe*v_e/v_s
    Msp = Mps

    Mew = -Dew*enz_star*wat_star
    Msw = -Dsw*sub_star*wat_star
    Mpw = -Dpw*prod_star*wat_star

    Mwe = Mew*v_w/v_e
    Mws = Msw*v_w/v_s
    Mwp = Mpw*v_w/v_s

    Mss = -(Mes+Mps+Mws)
    Mpp = -(Mep+Msp+Mwp)
    Mee = -(Mse+Mpe+Mwe)
    Mww = -(Mew+Msw+Mpw)

    return sub_star, prod_star, wat_star, RS, RP, RE, Mee, Mes, Mep, Mse, Mss, Msp, Mpe, Mps, Mpp, Mew, Msw, Mpw, Mwe, Mws, Mwp, Mww


def tern_test(sub_star, prod_star, wat_star, RS, RP, RE, Mse, Msp, Mpe, Mwe, Mws, Mwp, enz_star=0.1, v_e=10, v_s=1, v_w=0.5):
    q2 = 0.001

    Cprime = np.zeros([3, 3])
    Cprime[0, 0] = (Mse+Mpe+Mwe)*q2/enz_star + Mwe*v_e*q2/(v_w*wat_star)
    Cprime[0, 1] = -Mse*v_e*q2/(v_s*sub_star) + Mwe*v_e*q2/(v_w*wat_star)
    Cprime[0, 2] = -Mpe*v_e*q2/(v_s*prod_star) + Mwe*v_e*q2/(v_w*wat_star)

    Cprime[1, 0] = -Mse*q2/enz_star - RE + Mws*v_s*q2/(v_w*wat_star)
    Cprime[1, 1] = (Mse*v_e/v_s+Msp+Mws)*q2/sub_star - \
        RS + Mws*v_s*q2/(v_w*wat_star)
    Cprime[1, 2] = -Msp*q2/prod_star + RP + Mws*v_s*q2/(v_w*wat_star)

    Cprime[2, 0] = -Mpe*q2/enz_star + RE + Mwp*v_s*q2/(v_w*wat_star)
    Cprime[2, 1] = -Msp*q2/sub_star + RS + Mwp*v_s*q2/(v_w*wat_star)
    Cprime[2, 2] = (Mpe*v_e/v_s+Msp+Mwp)*q2/prod_star - \
        RP + Mwp*v_s*q2/(v_w*wat_star)

    return np.linalg.det(Cprime)/(q2*q2)


def tern_inst(del_mu=1, del_e=1, enz_star=0.1, k_spo=1, k_cat=1, Dpe=1, Dse=1, Dps=10, v_e=10, v_s=1, v_w=0.5, alpha=0.1, Dew=10, Dsw=10, Dpw=10):

    if alpha + enz_star > 1:
        return np.nan
    sub_star, prod_star, wat_star, RS, RP, RE, Mee, Mes, Mep, Mse, Mss, Msp, Mpe, Mps, Mpp, Mew, Msw, Mpw, Mwe, Mws, Mwp, Mww = get_steady_te(
        del_mu, del_e, enz_star, k_spo, k_cat, Dpe, Dse, Dps, v_e, v_s, v_w, alpha, Dew, Dsw, Dpw)
    val = tern_test(sub_star, prod_star, wat_star, RS, RP, RE,
                    Mse, Msp, Mpe, Mwe, Mws, Mwp, enz_star, v_e, v_s, v_w)

    return val


def tern_inst_scale(del_mu=1, del_e=1, enz_star=0.1, k_spo=1, k_cat=1, Dpe=1, Dse=1, Dps=10, v_e=10, v_s=1, v_w=0.5, alpha=0.1, Dew=10, Dsw=10, Dpw=10):

    Dwe, Dws, Dwp = Dew*v_w/v_e, Dsw*v_w/v_s, Dpw*v_w/v_s

    sub_star, prod_star, wat_star, RS, RP, RE, _, _, _,  _, _, _, _, _, _, _, _, _, _, _, _, _ = get_steady_te(
        del_mu, del_e, enz_star, k_spo, k_cat, Dpe, Dse, Dps, v_e, v_s, v_w, alpha, Dew, Dsw, Dpw)

    detC = -Dpe*(Dwe*enz_star*v_e + Dwp*prod_star*v_s + Dws*sub_star *
                 v_s)*(enz_star*(RE - RS)*v_e - prod_star*(RP + RS)*v_s)
    detC += Dse*(Dwe*enz_star*v_e + Dwp*prod_star*v_s + Dws*sub_star *
                 v_s)*(enz_star*(RE + RP)*v_e + sub_star*(RP + RS)*v_s)
    detC += Dpe*wat_star*(-Dws*enz_star*RE*v_e + Dwe*enz_star *
                          RS*v_e + Dws*prod_star*RP*v_s + Dwp*prod_star*RS*v_s)*v_w
    detC += Dse*wat_star*(Dwp*enz_star*RE*v_e + Dwe*enz_star *
                          RP*v_e + Dws*sub_star*RP*v_s + Dwp*sub_star*RS*v_s)*v_w
    detC += Dwe*wat_star*v_s*(Dws*(enz_star*(RE + RP)*v_e + sub_star*(RP + RS)*v_s + wat_star*RP*v_w) +
                              Dwp*(enz_star*(-RE + RS)*v_e + prod_star*(RP + RS)*v_s + wat_star*RS*v_w))

    detC = detC*-1/(v_w*v_s)

    return detC


def inst_cond(del_mu=1, del_e=1, enz_star=0.1, k_spo=1, k_cat=1, Dpe=1, Dse=1, v_rat=10):
    sub_star, prod_star, RS, RP, RE, Mee, Mes, Mep, Mse, Mss, Msp, Mpe, Mps, Mpp = get_steady_bi(
        del_mu, del_e, enz_star, k_spo, k_cat, Dpe, Dse, v_rat)
    # refers to eq 34 in report (times v_P)
    alpha = 1-enz_star
    LHS = 1/v_rat
    gammaP = Mpe/(Mse+Mpe)
    gammaS = Mse/(Mse+Mpe)
    RHS = enz_star/alpha*(RE*((gammaP/RS)-(gammaS/RP))-1)
    return LHS - RHS  # negative when unstable  

def inst_cond_chi(del_mu=1, del_e=1, enz_star=0.1, k_spo=1, k_cat=1, Dpe=1, Dse=1, v_rat=10, chi=1):
    sub_star, prod_star, RS, RP, RE, Mee, Mes, Mep, Mse, Mss, Msp, Mpe, Mps, Mpp = get_steady_bi(
        del_mu, del_e, enz_star, k_spo, k_cat, Dpe, Dse, v_rat)
    # refers to eq 34 in report (times v_P)
    alpha = 1-enz_star
    LHS = 1/v_rat + chi
    gammaP = Mpe/(Mse+Mpe)
    gammaS = Mse/(Mse+Mpe)
    RHS = enz_star/alpha*(RE*((gammaP/RS)-(gammaS/RP))-1)
    return LHS - RHS  # negative when unstable  


def Rp(del_mu=1, del_e=1, enz_star=0.1, k_spo=1, k_cat=1):
    return k_spo + k_cat*enz_star*np.exp(-del_e-del_mu)


def Rs(del_e=1, enz_star=0.1, k_spo=1, k_cat=1):
    return k_spo*np.exp(del_e) + k_cat*enz_star


def Re(del_mu=1, del_e=1, enz_star=0.1, alpha=0.9, k_spo=1, k_cat=1):
    num = alpha*k_cat*k_spo*(1-np.exp(-del_mu))
    denom = Rp(del_mu, del_e, enz_star, k_spo, k_cat)
    denom += Rs(del_e, enz_star, k_spo, k_cat)
    return num/denom


def binodals(enz_start=1e-8, init_step=1e-4, sample_step=1e-7, e_acc=1e-7, del_mu=10, del_e=-5, k_spo=1, k_cat=1, Dpe=10, Dse=1, v_rat=10):

    enz_sample = np.arange(0, 1, sample_step)
    fs = f_eff(del_mu, del_e, enz_sample, k_spo, k_cat, Dpe, Dse, v_rat)

    elow = 0
    ehi = 1

    enz_check = enz_start
    to_right = (enz_sample > enz_check)

    mu = mu_eff(del_mu, del_e, enz_check, k_spo, k_cat, Dpe, Dse, v_rat)
    f = f_eff(del_mu, del_e, enz_check, k_spo, k_cat, Dpe, Dse, v_rat)

    linear = (enz_sample-enz_check)*mu + f

    while not (linear[to_right] > fs[to_right]).any():

        elow = enz_check
        enz_check += init_step
        to_right = (enz_sample > enz_check)

        mu = mu_eff(del_mu, del_e, enz_check, k_spo, k_cat, Dpe, Dse, v_rat)
        f = f_eff(del_mu, del_e, enz_check, k_spo, k_cat, Dpe, Dse, v_rat)
        linear = (enz_sample-enz_check)*mu + f

        if enz_check > 1:
            return (0, 1), (0, 1)

    ehi = enz_check

    hi_right = (enz_sample > ehi)
    enz_right_hi = enz_sample[(enz_sample > ehi)]
    cross_inds = np.where(linear[to_right] > fs[to_right])[0]

    low_cross, hi_cross = enz_right_hi[cross_inds[0]
                                       ], enz_right_hi[cross_inds[-1]]

    while hi_cross-low_cross > e_acc:

        enz_check = (ehi+elow)/2
        to_right = (enz_sample > enz_check)

        mu = mu_eff(del_mu, del_e, enz_check, k_spo, k_cat, Dpe, Dse, v_rat)
        f = f_eff(del_mu, del_e, enz_check, k_spo, k_cat, Dpe, Dse, v_rat)
        linear = (enz_sample-enz_check)*mu + f

        if (linear[to_right] > fs[to_right]).any():
            ehi = enz_check
        else:
            elow = enz_check

        linearhi = (enz_sample-ehi)*mu + f
        hi_right = (enz_sample > ehi)
        enz_right_hi = enz_sample[(enz_sample > ehi)]
        cross_inds = np.where(linearhi[to_right] > fs[to_right])[0]

        low_cross, hi_cross = enz_right_hi[cross_inds[0]
                                           ], enz_right_hi[cross_inds[-1]]

    return (elow, ehi), (low_cross, hi_cross)


def find_binodals(enz_start=1, init_step=1e-4, sample_step=1e-7, e_acc=1e-5, del_mu=10, del_e=-5, k_spo=1, k_cat=1, Dpe=10, Dse=1, v_rat=10, f_thresh=1e-10):

    enz_sample = np.arange(0, 1, sample_step)
    fs = f_eff(del_mu, del_e, enz_sample, k_spo, k_cat, Dpe, Dse, v_rat)

    elow = 0
    ehi = 1

    enz_check = enz_start
    to_left = (enz_sample < enz_check)

    mu = mu_eff(del_mu, del_e, enz_check, k_spo, k_cat, Dpe, Dse, v_rat)
    f = f_eff(del_mu, del_e, enz_check, k_spo, k_cat, Dpe, Dse, v_rat)

    linear = (enz_sample-enz_check)*mu + f
    
    # while not (linear[to_left] > fs[to_left]).any():
    while not np.nanmin(fs[to_left] - linear[to_left]) < -1*f_thresh:

        ehi = enz_check
        enz_check -= init_step

        mu = mu_eff(del_mu, del_e, enz_check, k_spo, k_cat, Dpe, Dse, v_rat)
        f = f_eff(del_mu, del_e, enz_check, k_spo, k_cat, Dpe, Dse, v_rat)
        linear = (enz_sample-enz_check)*mu + f
        to_left = (enz_sample < enz_check)

        if enz_check < 0 or not np.array(to_left).any():
            return (0, 1), (0, 1)

    elow = enz_check

    while ehi-elow > e_acc:
    # want to split these up
    # start with ehi-elow
    # then effectvieyl repeat with ehi,elow = crosshi,crosslow
        enz_check = (ehi+elow)/2
        to_left = (enz_sample < elow)

        mu = mu_eff(del_mu, del_e, enz_check, k_spo, k_cat, Dpe, Dse, v_rat)
        f = f_eff(del_mu, del_e, enz_check, k_spo, k_cat, Dpe, Dse, v_rat)
        linear = (enz_sample-enz_check)*mu + f

        # if (linear[to_left] > fs[to_left]).any():
        if np.nanmin(fs[to_left] - linear[to_left]) < -1*f_thresh:
            elow = enz_check
        else:
            ehi = enz_check

        if enz_check != elow:
            to_left = (enz_sample < elow)
            mu = mu_eff(del_mu, del_e, elow, k_spo, k_cat, Dpe, Dse, v_rat)
            f = f_eff(del_mu, del_e, elow, k_spo, k_cat, Dpe, Dse, v_rat)
            linear = (enz_sample-elow)*mu + f

        enz_left = enz_sample[to_left]
        cross_inds = np.where(linear[to_left] > fs[to_left])[0]
        crosslow, crosshi = enz_left[cross_inds[0]], enz_left[cross_inds[-1]]
    
    uppers = (elow, ehi)
    elow, ehi = crosslow, crosshi

    while ehi-elow > e_acc:
        enz_check = (ehi+elow)/2
        to_right = (enz_sample > elow)

        mu = mu_eff(del_mu, del_e, enz_check, k_spo, k_cat, Dpe, Dse, v_rat)
        f = f_eff(del_mu, del_e, enz_check, k_spo, k_cat, Dpe, Dse, v_rat)
        linear = (enz_sample-enz_check)*mu + f

        # if (linear[to_left] > fs[to_left]).any():
        if np.nanmin(fs[to_right] - linear[to_right]) < -1*f_thresh:
            ehi = enz_check
        else:
            elow = enz_check
        
    lowers = (elow, ehi)
    
    return uppers, lowers


def read_sol(data_file):
    
    whole_file = open(data_file, 'r').read()
    data_array = np.array([row.split(',') for row in whole_file.split('\n')[1:-2]])
    cols = whole_file.split('\n')[0].split(',')
    sol = {}
    for colname, col_no in zip(cols, range(len(cols))):
        sol[colname] = data_array[:, col_no]
        
    system = {}
    for param in whole_file.split('\n')[-1].split(','):
        label, value = param.split(': ')
        system[label.replace(' ', '')] = float(value)
        
    return sol, system