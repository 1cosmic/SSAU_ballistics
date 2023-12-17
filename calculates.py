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

P_vacuum = 0
P_ud = 0
P_ud_vacuum = 0

Vk = 0
H0 = 0
H = 0
i = 0
fi_0 = 0
A = 0
Vw = 0
Vk = 0

m_correct = 0

def ballistic_parameters(rocket):
    global m, mt, z, P_midel, n0, k_hight_engine, dmt, tk, Vx, P_vacuum, P_ud, P_ud_vacuum, Vk, H0, i, fi_0, A, Vw, Vk, A0

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

    P_vacuum = rocket['P_vacuum']
    P_ud = rocket['P_ud']
    P_ud_vacuum = rocket['P_ud_vacuum']

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
    global H, H0, i, fi_0, V0, A, A0, Vw, Vk

    H = ballistic_data['target_H']
    H0 = ballistic_data['H']
    i = radians(ballistic_data['i'])
    fi_0 = radians(ballistic_data['fi_0'])

    V0 = sqrt(mu / (R + H0))
    A = asin(cos(i) / cos(fi_0))
    Vw = w_Earth * R * cos(fi_0)

    Vk = sqrt(V0 ** 2 + Vw ** 2 - 2 * V0 * Vw * sin(A))
    A0 = acos((V0 / Vk) * cos(A))

    return Vk

def correct_mass_fuel():

    dV = 1
    while abs(dV) > 0.1:
        # Vmax = float(input("Расчитайте Vmax и введите его: "))
        Vmax = 7.59

        dV = Vmax - Vk
        z2_correct = z[1] * exp(-(dV / P_ud_vacuum[1]))
        m_correct = ((z2_correct - 1) / z2_correct) * m[1]
        tk_correct = m_correct / dmt[1]

        dV = round(dV, 3)
        z2_correct = round(z2_correct, 3)
        m_correct = round(m_correct, 3)
        tk_correct = round(tk_correct, 3)

        print(f"\ndV: {dV}\nz2: {z2_correct}\nm_correct: {m_correct}\ntk_correct: {tk_correct}")
        print("Ошибка в скорости выведения: ", dV)

    else:
        print("Погрешность удовлетворяет 0.1 м/с")
        return m_correct


def swap_start_to_iner():

    # Xk = float(input("Введите Xk при Vmax: "))
    # Yk = float(input("Введите Yk: "))
    # Uk = float(input("Введите Uk: "))
    # Wk = float(input("Введите Wk: "))
    Xk = 887.994
    Yk = 128.115
    Uk = 7.521
    Wk = -1.028

    x0 = -Xk * cos(A0) * sin(fi_0) + (R + Yk) * cos(fi_0)
    y0 = Xk * sin(A0)
    z0 = Xk * cos(A0) * cos(fi_0) + (R + Yk) * sin(fi_0)

    r0 = sqrt(x0 **2 + y0 **2 + z0 **2)
    dr = r0 - R

    Vx0 = -Uk * cos(A0) * sin(fi_0) + Wk * cos(fi_0) - w_Earth * y0
    Vy0 = Uk * sin(A0) + w_Earth * x0
    Vz0 = Uk * cos(A0) * cos(fi_0) + Wk * sin(fi_0)

    V0_by_inert = sqrt(Vx0 **2 + Vy0 **2 + Vz0 **2)
    dV = V0 - V0_by_inert
    tetta0 = asin( (x0 * Vx0 + y0 * Vy0 + z0 * Vz0) / (r0 * V0_by_inert))

    C1 = y0 * Vz0 - z0 * Vy0
    C2 = z0 * Vx0 - x0 * Vz0
    C3 = x0 * Vy0 - y0 * Vx0
    C = sqrt(C1 **2 + C2 **2 + C3 **2)

    tg_omega0 = C1 / -C2
    omega0 = atan(tg_omega0)

    cos_i = C3 / C
    i_by_inert = acos(cos_i)

    v0_by_inert = (V0_by_inert **2 * r0) / (mu)
    A = r0 / (2 - v0_by_inert)

    e = sqrt(1 + (v0_by_inert - 2) * v0_by_inert * pow(cos(tetta0), 2))
    p = A * (1 - (e **2))

    tg_vu = (v0_by_inert * sin(tetta0) * cos(tetta0)) / (v0_by_inert * pow(cos(tetta0), 2) -1)
    vu = atan(tg_vu)

    tg_u = z0 / (sin(i) * (x0 * cos(omega0) + y0 * sin(omega0)))
    u = atan(tg_u)
    w = u - vu

    print(f"\nX0: {x0} Y0: {y0} Z0: {z0} ")
    print(f"r0: {r0} разница высот: {dr}, сравните с H0!: {H0}")


    # Exrection from NOO to ragrered orbit:
    r1 = r0
    r2 = R + H
    ra = r2

    Vc1 = sqrt(mu / r1)
    Vc2 = sqrt(mu / r2)

    Vp1 = sqrt((2 * ra * mu) / (r1 * (r1 + ra)))
    Vp2 = sqrt((2 * ra * mu) / (r2 * (r2 + ra)))
    Vp1 = sqrt((2 * r1 * mu) / (ra * (r1 + ra)))
    Vp1 = sqrt((2 * r2 * mu) / (ra * (r2 + ra)))


