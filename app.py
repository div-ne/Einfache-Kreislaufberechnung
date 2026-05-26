import numpy as np
import pandas as pd
import streamlit as st
import CoolProp.CoolProp as cp
import io
from scipy.optimize import brentq

st.set_page_config(page_title="Einfache Kreislaufberechnung", layout="wide")

APP_TITLE = "Einfache Kreislaufberechnung"
APP_VERSION = "0.12.0V"

FLUIDS = {
    "R11": "R11",
    "R12": "R12",
    "R13": "R13",
    "R13I1": "R13I1",
    "R14": "R14",
    "R21": "R21",
    "R22": "R22",
    "R23": "R23",
    "R32": "R32",
    "R40": "R40",
    "R41": "R41",
    "R113": "R113",
    "R114": "R114",
    "R115": "R115",
    "R116": "R116",
    "R123": "R123",
    "R124": "R124",
    "R125": "R125",
    "R141b": "R141b",
    "R142b": "R142b",
    "R143a": "R143a",
    "R152A": "R152A",
    "R161": "R161",
    "R218": "R218",
    "R227EA": "R227EA",
    "R236EA": "R236EA",
    "R236FA": "R236FA",
    "R245ca": "R245ca",
    "R245fa": "R245fa",
    "R290": "Propane",
    "RC318": "RC318",
    "R365MFC": "R365MFC",
    "R401A": "R401A",
    "R401B": "R401B",
    "R401C": "R401C",
    "R402A": "R402A",
    "R402B": "R402B",
    "R403A": "R403A",
    "R403B": "R403B",
    "R404A": "R404A",
    "R405A": "R405A.mix",
    "R406A": "R406A.mix",
    "R407A": "R407A.mix",
    "R407B": "R407B.mix",
    "R407C": "R407C",
    "R407D": "R407D.mix",
    "R407E": "R407E.mix",
    "R407F": "R407F.mix",
    "R408A": "R408A.mix",
    "R409A": "R409A.mix",
    "R409B": "R409B.mix",
    "R410A": "R410A",
    "R410B": "R410B.mix",
    "R411A": "R411A.mix",
    "R411B": "R411B.mix",
    "R412A": "R412A.mix",
    "R413A": "R413A.mix",
    "R414A": "R414A.mix",
    "R414B": "R414B.mix",
    "R415A": "R415A.mix",
    "R415B": "R415B.mix",
    "R416A": "R416A.mix",
    "R417A": "R417A.mix",
    "R417B": "R417B.mix",
    "R417C": "R417C.mix",
    "R418A": "R418A.mix",
    "R419A": "R419A.mix",
    "R419B": "R419B.mix",
    "R420A": "R420A.mix",
    "R421A": "R421A.mix",
    "R421B": "R421B.mix",
    "R422A": "R422A.mix",
    "R422B": "R422B.mix",
    "R422C": "R422C.mix",
    "R422D": "R422D.mix",
    "R422E": "R422E.mix",
    "R423A": "R423A.mix",
    "R424A": "R424A.mix",
    "R425A": "R425A.mix",
    "R426A": "R426A.mix",
    "R427A": "R427A.mix",
    "R428A": "R428A.mix",
    "R429A": "R429A.mix",
    "R430A": "R430A.mix",
    "R431A": "R431A.mix",
    "R432A": "R432A.mix",
    "R433A": "R433A.mix",
    "R433B": "R433B.mix",
    "R433C": "R433C.mix",
    "R434A": "R434A.mix",
    "R435A": "R435A.mix",
    "R436A": "R436A.mix",
    "R436B": "R436B.mix",
    "R437A": "R437A.mix",
    "R438A": "R438A.mix",
    "R439A": "R439A.mix",
    "R440A": "R440A.mix",
    "R441A": "R441A.mix",
    "R442A": "R442A.mix",
    "R443A": "R443A.mix",
    "R444A": "R444A.mix",
    "R444B": "R444B.mix",
    "R445A": "R445A.mix",
    "R446A": "R446A.mix",
    "R447A": "R447A.mix",
    "R448A": "R448A.mix",
    "R449A": "R449A.mix",
    "R449B": "R449B.mix",
    "R450A": "R450A.mix",
    "R451A": "R451A.mix",
    "R451B": "R451B.mix",
    "R452A": "R452A.mix",
    "R453A": "R453A.mix",
    "R454A": "R454A.mix",
    "R454B": "R454B.mix",
    "R500": "R500.mix",
    "R501": "R501.mix",
    "R502": "R502.mix",
    "R503": "R503.mix",
    "R504": "R504.mix",
    "R507A": "R507A",
    "R508A": "R508A.mix",
    "R508B": "R508B.mix",
    "R509A": "R509A.mix",
    "R510A": "R510A.mix",
    "R511A": "R511A.mix",
    "R512A": "R512A.mix",
    "R513A": "R513A.mix",
    "R600": "n-Butane",
    "R601": "n-Pentane",
    "R600a": "IsoButane",
    "R601a": "Isopentane",
    "R702": "Hydrogen",
    "R717": "Ammonia",
    "R729": "Air",
    "R764": "SulfurDioxide",
    "R1233zd(E)": "R1233zd(E)",
    "R1234yf": "R1234yf",
    "R1234ze(E)": "R1234ze(E)",
    "R1234ze(Z)": "R1234ze(Z)",
    "R1243zf": "R1243zf",
    "R1336mzz(E)": "R1336mzz(E)",
}

def getVi(di):
    Vi = np.pi * (di) ** 2 * 2.5 * 1e-7
    return Vi

def getAM(da):
    AM = np.pi * da * 1e-3
    return AM

cu_pipes = [
    ("6x1", 6, 1, 4, getVi(4), getAM(6)),
    ("8x1", 8, 1, 6, getVi(6), getAM(8)),
    ("10x1", 10, 1, 8, getVi(8), getAM(10)),
    ("12x1", 12, 1, 10, getVi(10), getAM(12)),
    ("15x1", 15, 1, 13, getVi(13), getAM(15)),
    ("16x1", 16, 1, 14, getVi(14), getAM(16)),
    ("18x1", 18, 1, 16, getVi(16), getAM(18)),
    ("22x1", 22, 1, 20, getVi(20), getAM(22)),
    ("28x1", 28, 1, 26, getVi(26), getAM(28)),
    ("28x1.5", 28, 1.5, 25, getVi(25), getAM(28)),
    ("35x1.5", 35, 1.5, 32, getVi(32), getAM(35)),
    ("42x1.5", 42, 1.5, 39, getVi(39), getAM(42)),
    ("54x1.5", 54, 1.5, 51, getVi(51), getAM(54)),
    ("54x2", 54, 2, 50, getVi(50), getAM(54)),
    ("64x2", 64, 2, 60, getVi(60), getAM(64)),
    ("76.1x2", 76.1, 2, 74.1, getVi(72.1), getAM(76.1)),
    ("88.9x2", 88.9, 2, 84.9, getVi(84.9), getAM(88.9)),
    ("108x2.5", 108, 2.5, 103, getVi(103), getAM(108)),
    ("133x3", 133, 3, 127, getVi(127), getAM(133))
]

def friction_coefficient(di, w, nu):
    Re = di * w / nu
    k = 0.0015 / 1e3
    epsilon_k = k / di
    lambda_hagenpoiseulle = 64 / Re
    lambda_blasius = 0.3164 / Re**0.25
    lambda_nikuradse = (-2 * np.log10(k / 3.71 / di))**-2
    lambda_prandtl = 0.02
    for _ in range(10):
        lambda_prandtl = (2 * np.log10(Re * np.sqrt(lambda_prandtl)))**-2
    lambda_colebrookwhite = 0.02
    for _ in range(10):
        lambda_colebrookwhite = (-2 * np.log10(2.51 / Re / np.sqrt(lambda_colebrookwhite) + k / 3.71 / di))**-2
    check_moody_diagram = Re * np.sqrt(lambda_nikuradse) * k / di
    if Re < 2320:
        return lambda_hagenpoiseulle
    if check_moody_diagram >= 200:
        return lambda_nikuradse
    if epsilon_k < 0.001 and Re < 10000:
        return lambda_blasius
    if epsilon_k < 0.0002 and Re < 100000:
        return lambda_blasius
    if epsilon_k < 0.00002 and Re < 1000000:
        return lambda_prandtl
    if epsilon_k < 0.00001:
        return lambda_prandtl
    return lambda_colebrookwhite

def find_next_bigger_pipe(di):
    di_values = [pipe[3] for pipe in cu_pipes]
    bigger = [x for x in di_values if x > di * 1e3]
    return min(bigger) * 1e-3 if bigger else None

def find_next_smaller_pipe(di):
    di_values = [pipe[3] for pipe in cu_pipes]
    smaller = [x for x in di_values if x < di * 1e3]
    return max(smaller) * 1e-3 if smaller else None

def pipe_geometry(di_m, length_m):
    key = round(di_m * 1e3, 1)
    row = next(x for x in cu_pipes if abs(x[3] - key) < 1e-6)
    return row[0], round(row[5] * length_m, 3), round(row[4] * length_m * 1e3, 2)

def project_pipe(pre_di, length, V, nu, rho, dp_max, strategy="pressure"):
    di_values = sorted(pipe[3] * 1e-3 for pipe in cu_pipes)

    def calc(di):
        w = 4 * V / (np.pi * di**2)
        lam = float(friction_coefficient(di, w, nu))
        dp = lam * length * rho * w**2 * 0.5 / di
        if rho > 900:
            jg = np.sqrt(0.85 * 9.80665 * di * abs(rho * 1.1 - rho) / rho)
        else:
            jg = np.sqrt(0.85 * 9.80665 * di * abs(900 - rho) / rho)
        return w, jg, dp

    valid = []
    for di in di_values:
        w, jg, dp = calc(di)
        if w > jg and dp < dp_max:
            valid.append((di, w, jg, dp))

    if valid:
        return valid[-1] if strategy == "pressure" else valid[0]

    di = min(di_values, key=lambda x: abs(x - pre_di))
    w, jg, dp = calc(di)
    return di, w, jg, dp

def T_brentq_s(p, s, T_min, T_max, fluid):
    def f(T):
        return cp.PropsSI("S", "T", T, "P", p, fluid) - s
    return brentq(f, T_min, T_max)

def T_brentq_h(p, h, T_min, T_max, fluid):
    def f(T):
        return cp.PropsSI("H", "T", T, "P", p, fluid) - h
    return brentq(f, T_min, T_max)

def run_calculation(inputs):
    project = inputs["project"]
    mode = inputs["mode"]
    fluid = FLUIDS.get(inputs["fluid"], inputs["fluid"])
    errors = []

    if mode == 0:
        Q_0 = inputs["Q_0"] * 1e3
        Q_c = 0.0
        V_com = 0.0
    elif mode == 1:
        Q_c = inputs["Q_c"] * 1e3
        Q_0 = 0.0
        V_com = 0.0
        if Q_c > 0:
            Q_c = -Q_c
    else:
        V_com = inputs["V_com"] / 3.6e3
        Q_0 = 0.0
        Q_c = 0.0

    T_c = inputs["T_c"] + 273.15
    T_0 = inputs["T_0"] + 273.15
    dt_cu = inputs["dt_cu"]
    dt_0h = inputs["dt_0h"]
    dt_sh = inputs["dt_sh"]

    if inputs["if_pipes"] == 0:
        l_hg = l_fl = l_sl = 1.0
    else:
        l_hg = inputs["l_hg"]
        l_fl = inputs["l_fl"]
        l_sl = inputs["l_sl"]

    if T_c < T_0:
        raise ValueError("Die Verflüssigungstemperatur darf nicht unter der Verdampfungstemperatur liegen.")
    if (T_c - T_0) < 20:
        errors.append("Die Differenz zwischen der Verflüssigungstemperatur und Verdampfungstemperatur ist sehr gering.")
    if Q_0 == 0 and mode == 0:
        Q_0 = 1
    if Q_c == 0 and mode == 1:
        Q_c = 1
    if V_com == 0 and mode == 2:
        V_com = 1

    p_0_dew = cp.PropsSI("P", "T", T_0, "Q", 1, fluid)
    p_0_bubble = cp.PropsSI("P", "T", T_0, "Q", 0, fluid)
    p_0 = (p_0_dew + p_0_bubble) / 2 if fluid[1] == "4" else p_0_dew
    T_0_dew = cp.PropsSI("T", "P", p_0, "Q", 1, fluid)
    T_0_bubble = cp.PropsSI("T", "P", p_0, "Q", 0, fluid)
    T_4 = T_0_dew
    h_4 = cp.PropsSI("H", "P", p_0, "Q", 0.5, fluid)

    T_0h = T_0 + dt_0h
    T_sh = T_0 + dt_0h + dt_sh
    p_c = cp.PropsSI("P", "T", T_c, "Q", 1, fluid)
    T_c_dew = cp.PropsSI("T", "P", p_c, "Q", 1, fluid)
    T_c_bubble = cp.PropsSI("T", "P", p_c, "Q", 0, fluid)
    T_cu = T_c_bubble - dt_cu
    h_3 = cp.PropsSI("H", "T", T_c_bubble, "Q", 0, fluid)
    h_4 = h_3
    T_2min = T_c_dew + 1
    T_2max = T_c_dew + 100

    h_1 = cp.PropsSI("H", "P", p_0, "T", T_sh, fluid)
    s_1 = cp.PropsSI("S", "P", p_0, "T", T_sh, fluid)
    T_2s = T_brentq_s(p_c, s_1, T_2min, T_2max, fluid)
    h_2s = cp.PropsSI("H", "P", p_c, "T", T_2s, fluid)
    pi = p_c / p_0
    eta_is = -0.020644445 + 0.68403852 * pi - 0.22147167 * pi**2 + 0.032145926 * pi**3 - 0.00178 * pi**4 if pi <= 5 else 0.821 - 0.0105 * pi
    h_2 = h_1 + (h_2s - h_1) / eta_is
    h_5 = cp.PropsSI("H", "P", p_0, "T", T_0h, fluid)
    rho_com_in = cp.PropsSI("D", "P", p_0, "T", T_sh + 10, fluid)

    q_c = h_3 - h_2
    q_0 = h_5 - h_4
    if mode == 0:
        m_R = Q_0 / q_0
        Q_c = m_R * q_c
        V_com = m_R / rho_com_in
    elif mode == 1:
        m_R = Q_c / q_c
        Q_0 = m_R * q_0
        V_com = m_R / rho_com_in
    else:
        m_R = V_com * rho_com_in
        Q_c = m_R * q_c
        Q_0 = m_R * q_0

    P_com = m_R * (h_2 - h_1)
    COP = -Q_c / P_com
    EER = Q_0 / P_com

    rho_1 = cp.PropsSI("D", "P", p_0, "T", T_sh, fluid)
    T_end = T_brentq_h(p_c, h_2, T_2min, T_2max, fluid)
    rho_2 = cp.PropsSI("D", "P", p_c, "T", T_end, fluid)
    rho_3 = cp.PropsSI("D", "P", p_c, "T", T_cu, fluid)

    V_1 = m_R / rho_1
    V_2 = m_R / rho_2
    V_3 = m_R / rho_3

    mu_1 = cp.PropsSI("V", "P", p_0, "T", T_sh, fluid)
    mu_2 = cp.PropsSI("V", "P", p_c, "T", T_end, fluid)
    mu_3 = cp.PropsSI("V", "P", p_c, "T", T_cu, fluid)

    nu_1 = mu_1 / rho_1
    nu_2 = mu_2 / rho_2
    nu_3 = mu_3 / rho_3

    dp_max_hg = float(p_c - cp.PropsSI("P", "T", T_c_dew - 1, "Q", 1, fluid))
    dp_max_fl = float(p_c - cp.PropsSI("P", "T", T_c_bubble - 1, "Q", 1, fluid))
    dp_max_sl = float(cp.PropsSI("P", "T", T_0_dew + 1, "Q", 1, fluid) - p_0)

    pre_hg = float(find_next_bigger_pipe(np.sqrt(4 * V_2 / np.pi / 15)))
    pre_fl = float(find_next_bigger_pipe(np.sqrt(4 * V_3 / np.pi / 1)))
    pre_sl = float(find_next_bigger_pipe(np.sqrt(4 * V_1 / np.pi / 12)))

    di_hg_dp, w_hg_dp, jg_hg_dp, dp_hg_dp = project_pipe(pre_hg, l_hg * 1.4, V_2, nu_2, rho_2, dp_max_hg, strategy="pressure")
    di_hg_dm, w_hg_dm, jg_hg_dm, dp_hg_dm = project_pipe(pre_hg, l_hg * 1.4, V_2, nu_2, rho_2, 1e5, strategy="diameter")

    di_fl_dp, w_fl_dp, jg_fl_dp, dp_fl_dp = project_pipe(pre_fl, l_fl * 1.4, V_3, nu_3, rho_3, dp_max_fl, strategy="pressure")
    di_fl_dm, w_fl_dm, jg_fl_dm, dp_fl_dm = project_pipe(pre_fl, l_fl * 1.4, V_3, nu_3, rho_3, 1e5, strategy="diameter")

    di_sl_dp, w_sl_dp, jg_sl_dp, dp_sl_dp = project_pipe(pre_sl, l_sl * 1.4, V_1, nu_1, rho_1, dp_max_sl, strategy="pressure")
    di_sl_dm, w_sl_dm, jg_sl_dm, dp_sl_dm = project_pipe(pre_sl, l_sl * 1.4, V_1, nu_1, rho_1, 1e5, strategy="diameter")

    d_hg_dp, A_M_hg_dp, V_i_hg_dp = pipe_geometry(di_hg_dp, l_hg)
    d_hg_dm, A_M_hg_dm, V_i_hg_dm = pipe_geometry(di_hg_dm, l_hg)

    d_fl_dp, A_M_fl_dp, V_i_fl_dp = pipe_geometry(di_fl_dp, l_fl)
    d_fl_dm, A_M_fl_dm, V_i_fl_dm = pipe_geometry(di_fl_dm, l_fl)

    d_sl_dp, A_M_sl_dp, V_i_sl_dp = pipe_geometry(di_sl_dp, l_sl)
    d_sl_dm, A_M_sl_dm, V_i_sl_dm = pipe_geometry(di_sl_dm, l_sl)

    df = pd.DataFrame({
        "P": ["1", "2", "3", "4", "5", "c''", "c'", "0''", "0'"],
        "Temperatur [°C]": [
            round(T_sh - 273.15, 2), round(T_end - 273.15, 2), round(T_cu - 273.15, 2),
            round(T_4 - 273.15, 2), round(T_0h - 273.15, 2),
            round(T_c_dew - 273.15, 2), round(T_c_bubble - 273.15, 2),
            round(T_0_dew - 273.15, 2), round(T_0_bubble - 273.15, 2)
        ],
        "Druck [bar]": [
            round(p_0 / 1e5, 2), round(p_c / 1e5, 2), round(p_c / 1e5, 2),
            round(p_0 / 1e5, 2), round(p_0 / 1e5, 2),
            round(p_c / 1e5, 2), round(p_c / 1e5, 2),
            round(p_0 / 1e5, 2), round(p_0 / 1e5, 2)
        ],
        "Spezifische Enthalpie [kJ/kg]": [
            round(h_1 / 1e3, 2), round(h_2 / 1e3, 2), round(h_3 / 1e3, 2),
            round(h_4 / 1e3, 2), round(h_5 / 1e3, 2),
            round(cp.PropsSI("H", "P", p_c, "Q", 1, fluid) / 1e3, 2),
            round(cp.PropsSI("H", "P", p_c, "Q", 0, fluid) / 1e3, 2),
            round(cp.PropsSI("H", "P", p_0, "Q", 1, fluid) / 1e3, 2),
            round(cp.PropsSI("H", "P", p_0, "Q", 0, fluid) / 1e3, 2),
        ],
        "Dichte [kg/m3]": [
            round(rho_1, 2), round(rho_2, 2), round(rho_3, 2),
            round(cp.PropsSI("D", "P", p_0, "Q", 0.5, fluid), 2),
            round(cp.PropsSI("D", "P", p_0, "T", T_0h, fluid), 2),
            round(cp.PropsSI("D", "P", p_c, "Q", 1, fluid), 2),
            round(cp.PropsSI("D", "P", p_c, "Q", 0, fluid), 2),
            round(cp.PropsSI("D", "P", p_0, "Q", 1, fluid), 2),
            round(cp.PropsSI("D", "P", p_0, "Q", 0, fluid), 2),
        ],
        "Spezifische Entropie [kJ/kg/K]": [
            round(s_1 / 1e3, 4),
            round(cp.PropsSI("S", "P", p_c, "T", T_end, fluid) / 1e3, 4),
            round(cp.PropsSI("S", "P", p_c, "T", T_cu, fluid) / 1e3, 4),
            round(cp.PropsSI("S", "P", p_0, "Q", 0.5, fluid) / 1e3, 4),
            round(cp.PropsSI("S", "P", p_0, "T", T_0h, fluid) / 1e3, 4),
            round(cp.PropsSI("S", "P", p_c, "Q", 1, fluid) / 1e3, 2),
            round(cp.PropsSI("S", "P", p_c, "Q", 0, fluid) / 1e3, 2),
            round(cp.PropsSI("S", "P", p_0, "Q", 1, fluid) / 1e3, 2),
            round(cp.PropsSI("S", "P", p_0, "Q", 0, fluid) / 1e3, 2),
        ],
        "Dampfqualität [%]": ["", "", "", 50, "", 100, 0, 100, 0],
    })

    df_results = pd.DataFrame({
        "Parameter": [
            "Projekt", "Leistungsaufnahme [kW]", "Wärmeleistung [kW]",
            "Kälteleistung [kW]", "Kältemittelmassenstrom [kg/s]",
            "Verdichtervolumenstrom [m3/h]", "COP", "EER"
        ],
        "Wert": [project, round(P_com / 1e3, 2), round(Q_c / 1e3, 2), round(Q_0 / 1e3, 2), round(m_R, 6), round(V_com * 3.6e3, 2), round(COP, 2), round(EER, 2)]
    })

    df_pipes = pd.DataFrame([
        ["Heissgasleitung", "minimaler Druckverlust", d_hg_dp, l_hg, round(w_hg_dp, 2), round(jg_hg_dp, 2), round(dp_hg_dp / 1e5, 2), A_M_hg_dp, V_i_hg_dp],
        ["", "minimaler Durchmesser", d_hg_dm, l_hg, round(w_hg_dm, 2), round(jg_hg_dm, 2), round(dp_hg_dm / 1e5, 2), A_M_hg_dm, V_i_hg_dm],
        ["Flüssigkeitsleitung", "minimaler Druckverlust", d_fl_dp, l_fl, round(w_fl_dp, 2), round(jg_fl_dp, 2), round(dp_fl_dp / 1e5, 2), A_M_fl_dp, V_i_fl_dp],
        ["", "minimaler Durchmesser", d_fl_dm, l_fl, round(w_fl_dm, 2), round(jg_fl_dm, 2), round(dp_fl_dm / 1e5, 2), A_M_fl_dm, V_i_fl_dm],
        ["Saugleitung", "minimaler Druckverlust", d_sl_dp, l_sl, round(w_sl_dp, 2), round(jg_sl_dp, 2), round(dp_sl_dp / 1e5, 2), A_M_sl_dp, V_i_sl_dp],
        ["", "minimaler Durchmesser", d_sl_dm, l_sl, round(w_sl_dm, 2), round(jg_sl_dm, 2), round(dp_sl_dm / 1e5, 2), A_M_sl_dm, V_i_sl_dm],
    ], columns=[
        "Leitung", "Berechnungsmethode", "Durchmesser [mm]", "Länge [m]",
        "Strömungsgeschwindigkeit [m/s]", "Jacobs-Geschwindigkeit [m/s]",
        "Druckverlust [bar]", "Aussenmantelfläche [m2]", "Innenvolumen [dm3]"
    ])

    return df, df_results, df_pipes, errors

st.markdown(
    f"""
    <div style="display:flex; align-items:baseline; gap:0.6rem; flex-wrap:nowrap; white-space:nowrap; margin:0 0 1rem 0; overflow-x:auto;">
        <h1 style="margin:0; display:inline-block;">{APP_TITLE}</h1>
        <span style="font-size:1rem; color:#6b7280; display:inline-block;">{APP_VERSION}</span>
    </div>
    """,
    unsafe_allow_html=True,
)

col1, col2 = st.columns([1.05, 1.0])

with col1:
    row1_col1, row1_col2 = st.columns(2)
    with row1_col1:
        project = st.text_input("Projekt", "Projekt")
    with row1_col2:
        fluid = st.selectbox("Kältemittel", list(FLUIDS.keys()), index=list(FLUIDS.keys()).index("R290"))

    row2_col1, row2_col2 = st.columns(2)
    with row2_col1:
        mode = st.selectbox("Eingabemodus", ["Kälteleistung", "Wärmeleistung", "Verdichtervolumenstrom"])
    with row2_col2:
        show_q0 = mode == "Kälteleistung"
        show_qc = mode == "Wärmeleistung"
        show_vcom = mode == "Verdichtervolumenstrom"
        if show_q0:
            q0 = st.number_input("Kälteleistung [kW]", value=1.0)
            qc = 0.0
            vcom = 0.0
        elif show_qc:
            qc = st.number_input("Wärmeleistung [kW]", value=1.0)
            q0 = 0.0
            vcom = 0.0
        else:
            vcom = st.number_input("Verdichtervolumenstrom [m3/h]", value=1.0)
            q0 = 0.0
            qc = 0.0

    row3_col1, row3_col2 = st.columns(2)
    with row3_col1:
        tc = st.number_input("Verflüssigungstemperatur [°C]", value=45.0)
    with row3_col2:
        dtcu = st.number_input("Unterkühlung [K]", value=2.0)

    row4_col1, row4_col2 = st.columns(2)
    with row4_col1:
        t0 = st.number_input("Verdampfungstemperatur [°C]", value=5.0)
    with row4_col2:
        dt0h = st.number_input("Verdampferüberhitzung [K]", value=6.0)

    dtsh = st.number_input("Saugleitungsüberhitzung [K]", value=9.0)
    ifpipes = st.selectbox("Rohrleitungsdimensionierung", ["Nein", "Ja"])

    show_pipe_inputs = ifpipes == "Ja"
    lhg = st.number_input("Heissgasleitungslänge [m]", value=5.0) if show_pipe_inputs else 5.0
    lfl = st.number_input("Flüssigkeitsleitungslänge [m]", value=2.5) if show_pipe_inputs else 2.5
    lsl = st.number_input("Saugleitungslänge [m]", value=3.0) if show_pipe_inputs else 3.0

    run = st.button("Berechnen")

with col2:
    pass

if run:
    try:
        df, df_results, df_pipes, errors = run_calculation({
            "project": project,
            "fluid": fluid,
            "mode": ["Kälteleistung", "Wärmeleistung", "Verdichtervolumenstrom"].index(mode),
            "Q_0": q0,
            "Q_c": qc,
            "V_com": vcom,
            "T_c": tc,
            "T_0": t0,
            "dt_cu": dtcu,
            "dt_0h": dt0h,
            "dt_sh": dtsh,
            "if_pipes": 0 if ifpipes == "Nein" else 1,
            "l_hg": lhg,
            "l_fl": lfl,
            "l_sl": lsl,
        })

        csv_buffer = io.StringIO()
        csv_buffer.write("Systemdaten\n")
        df_results.to_csv(csv_buffer, index=False, sep=';')
        csv_buffer.write("\nKreislaufpunkte\n")
        df.to_csv(csv_buffer, index=False, sep=';')

        if ifpipes == "Ja":
            csv_buffer.write("\nRohrleitungsdimensionierung\n")
            df_pipes.to_csv(csv_buffer, index=False, sep=';')

        with col2:
            st.subheader("Systemdaten")
            st.dataframe(df_results, use_container_width=True, hide_index=True)
            st.download_button(
                "CSV-Datei erstellen",
                data=csv_buffer.getvalue().encode("utf-8-sig"),
                file_name=f"{project.replace(' ', '_')}_auswertung.csv",
                mime="text/csv"
            )

        st.subheader("Kreislaufpunkte")
        st.dataframe(df, use_container_width=True, hide_index=True)

        if ifpipes == "Ja":
            st.subheader("Rohrleitungsdimensionierung")
            st.dataframe(df_pipes, use_container_width=True, hide_index=True)

        with st.expander("Anleitung"):
            st.markdown(
                """
                Dieses Tool dient zur **einfachen Kreislaufberechnung** von Kälteanlagen und unterstützt dich bei der Bestimmung von **Kreislaufpunkten**, **Systemdaten**, **Stoffdaten** sowie der **Rohrleitungsdimensionierung**. Damit kannst du thermodynamische Zustände schnell überschlagen, Ergebnisse tabellarisch auswerten und bei Bedarf als CSV-Datei exportieren.

                **Bedeutung der Kreislaufpunkte:**

                - **1:** Verdichtereingang nach Überhitzung von Verdampfer und Saugleitung
                - **2:** Verdichterausgang und Verflüssigereintritt
                - **3:** Verflüssigeraustritt und Expansionsventileingang
                - **4:** Expansionsventilausgang und Verdampfereintritt
                - **5:** Verdampferausgang und Eingang Saugleitung
                - **c''**: Gesättigtes Gas bei Verflüssigungsdruck
                - **c'**: Flüssigkeit bei Verflüssigungsdruck
                - **0''**: Gesättigtes Gas bei Verdampfungsdruck
                - **0'**: Flüssigkeit bei Verdampfungsdruck

                Repository: https://github.com/div-ne/Einfache-Kreislaufberechnung
                """
            )

        if errors:
            st.warning("\n".join(errors))
    except Exception as e:
        st.error(str(e))