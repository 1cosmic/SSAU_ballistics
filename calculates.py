from math import *


# Earth cosmic consts.
g_const = 9.801
mu = 398600.44
w_Earth = 7.2921 * (10 ** -5)
R = 6371


# Program memory.
m = [0,0]
mt = [0,0]
z = [0,0]
P_midel = [0,0]
n0 = [0,0]
k_hight_engine = [0,0]
dmt = [0,0]
tk = [0,0]
Vx = [0,0]

def ballistic_parameters(rocket):

    m_start = rocket['m'][0] + rocket['m'][1] + rocket['m_pn']
    last_mass = m_start

    for i in range(2):
        m[i] = last_mass
        mt[i] = rocket['mt'][i]
        P0 = rocket['P0'][i]
        P_vacuum = rocket['P_vacuum'][i]
        P_ud = rocket['P_ud'][i]
        P_ud_vacuum = rocket['P_ud_vacuum'][i]
        D_midel = rocket['D_midel'][i]

        z[i] = m[i] / (m[i] - mt[i])
        P_midel[i] = m[i] / (pi * pow(D_midel / 2, 2)) if P_midel[i] != "-" else P_midel[i]
        n0[i] = P0 / (m[i] * g_const) if P0 != "-" else P_vacuum / (m[i] * g_const)
        k_hight_engine[i] = P_ud_vacuum / P_ud if P_ud != "-" else "-"

        dmt[i] = P_vacuum / P_ud_vacuum
        tk[i] = mt[i] / dmt[i]

        Vx[i] = P_ud_vacuum * log(z[i])

        # Decrement of mass the next step of Rocket.
        last_mass = last_mass - rocket['m'][i]


    # Reformat of return values.
    data = {}
    data['m_start'] = m_start
    data['m_step_1'] = m[0]
    data['m_step_2'] = m[1]

    data['z_1'] = z[0]
    data['z_2'] = z[1]

    data['P_midel_1'] = P_midel[0]
    data['P_midel_2'] = P_midel[1]

    data['n0_1'] = n0[0]
    data['n0_2'] = n0[1]

    data['k_hight_engine_1'] = k_hight_engine[0]
    data['k_hight_engine_2'] = k_hight_engine[1]

    data['dmt_1'] = dmt[0]
    data['dmt_2'] = dmt[1]

    data['tk_1'] = tk[0]
    data['tk_2'] = tk[1]

    data['Vx_1'] = Vx[0]
    data['Vx_2'] = Vx[1]
    data['V_sum'] = Vx[0] + Vx[1]

    # Round all values for human)
    for key in data:
        if data[key] != "-":
            data[key] = round(data[key], 3)

    return rocket | data

def calcs_Vx(ballistic_data):
    H0 = ballistic_data['H']
    i = radians(ballistic_data['i'])
    fi_0 = radians(ballistic_data['fi_0'])

    V0 = sqrt(mu / (R + H0))
    A = asin(cos(i) / cos(fi_0))
    Vw = w_Earth * R * cos(fi_0)

    Vk = sqrt(V0 ** 2 + Vw ** 2 - 2 * V0 * Vw * sin(A))

    return Vk