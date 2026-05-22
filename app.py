import numpy as np
import pandas as pd
import streamlit as st
import CoolProp.CoolProp as cp
from scipy.optimize import brentq

st.set_page_config(page_title="VCCSimple Streamlit", layout="wide")

APP_TITLE = "VCCSimple Streamlit"
APP_VERSION = "0.5.0V"

FLUIDS = {
    "R290": "Propane",
    "R1270": "Propylene",
    "R170": "Ethane",
    "R1150": "Ethylene",
    "R744": "CarbonDioxide",
    "CO2": "CarbonDioxide",
    "R717": "Ammonia",
    "R718": "Water",
    "R600": "n-Butane",
    "R600a": "IsoButane",
    "R601": "n-Pentane",
    "R601a": "Isopentane",
    "R702": "Hydrogen",
    "R729": "Air",
    "R764": "SulfurDioxide",
}

def getVi(di):
    Vi = np.pi * (di)**2 * 2.5 * 1e-7
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

def project_pipe(pre_di, length, V, nu, rho, dp_max):
    di = pre_di
    for _ in range(40):
        w = 4 * V / (np.pi * di**2)
        lam = float(friction_coefficient(di, w, nu))
        dp = lam * length * rho * w**2 * 0.5 / di
        if rho > 900:
            jg = np.sqrt(0.85 * 9.80665 * di * abs(rho * 1.1 - rho) / rho)
        else:
            jg = np.sqrt(0.85 * 9.80665 * di * abs(900 - rho) / rho)
        if w > jg and dp < dp_max:
            return di, w, jg, dp
        nd = find_next_bigger_pipe(di) if dp > dp_max else find_next_smaller_pipe(di)
        if nd is None or nd == di:
            return di, w, jg, dp
        di = nd
    return di, w, jg, dp

def T_brentq_s(p, s, T_min, T_max, fluid):
    def f(T):
        return cp.PropsSI('S', 'T', T, 'P', p, fluid) - s
    return brentq(f, T_min, T_max)

def T_brentq_h(p, h, T_min, T_max, fluid):
    def f(T):
        return cp.PropsSI('H', 'T', T, 'P', p, fluid) - h
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

    di_hg, w_hg, jg_hg, dp_hg = project_pipe(pre_hg, l_hg * 1.4, V_2, nu_2, rho_2, dp_max_hg)
    di_fl, w_fl, jg_fl, dp_fl = project_pipe(pre_fl, l_fl * 1.4, V_3, nu_3, rho_3, dp_max_fl)
    di_sl, w_sl, jg_sl, dp_sl = project_pipe(pre_sl, l_sl * 1.4, V_1, nu_1, rho_1, dp_max_sl)

    d_hg, A_M_hg, V_i_hg = pipe_geometry(di_hg, l_hg)
    d_fl, A_M_fl, V_i_fl = pipe_geometry(di_fl, l_fl)
    d_sl, A_M_sl, V_i_sl = pipe_geometry(di_sl, l_sl)

    df = pd.DataFrame({
        "P": ["1", "2", "3", "4", "5", "c''", "c'", "0''", "0'"],
        "Temperatur [°C]": [round(T_sh - 273.15, 2), round(T_end - 273.15, 2), round(T_cu - 273.15, 2), round(T_4 - 273.15, 2), round(T_0h - 273.15, 2), round(T_c_dew - 273.15, 2), round(T_c_bubble - 273.15, 2), round(T_0_dew - 273.15, 2), round(T_0_bubble - 273.15, 2)],
        "Druck [bar]": [round(p_0 / 1e5, 2), round(p_c / 1e5, 2), round(p_c / 1e5, 2), round(p_0 / 1e5, 2), round(p_0 / 1e5, 2), round(p_c / 1e5, 2), round(p_c / 1e5, 2), round(p_0 / 1e5, 2), round(p_0 / 1e5, 2)],
        "Spezifische Enthalpie [kJ/kg]": [round(h_1 / 1e3, 2), round(h_2 / 1e3, 2), round(h_3 / 1e3, 2), round(h_4 / 1e3, 2), round(h_5 / 1e3, 2), round(cp.PropsSI("H", "P", p_c, "Q", 1, fluid) / 1e3, 2), round(cp.PropsSI("H", "P", p_c, "Q", 0, fluid) / 1e3, 2), round(cp.PropsSI("H", "P", p_0, "Q", 1, fluid) / 1e3, 2), round(cp.PropsSI("H", "P", p_0, "Q", 0, fluid) / 1e3, 2)],
        "Dichte [kg/m3]": [round(rho_1, 2), round(rho_2, 2), round(rho_3, 2), round(cp.PropsSI("D", "P", p_0, "Q", 0.5, fluid), 2), round(cp.PropsSI("D", "P", p_0, "T", T_0h, fluid), 2), round(cp.PropsSI("D", "P", p_c, "Q", 1, fluid), 2), round(cp.PropsSI("D", "P", p_c, "Q", 0, fluid), 2), round(cp.PropsSI("D", "P", p_0, "Q", 1, fluid), 2), round(cp.PropsSI("D", "P", p_0, "Q", 0, fluid), 2)],
        "Spezifische Entropie [kJ/kg/K]": [round(s_1 / 1e3, 4), round(cp.PropsSI("S", "P", p_c, "T", T_end, fluid) / 1e3, 4), round(cp.PropsSI("S", "P", p_c, "T", T_cu, fluid) / 1e3, 4), round(cp.PropsSI("S", "P", p_0, "Q", 0.5, fluid) / 1e3, 4), round(cp.PropsSI("S", "P", p_0, "T", T_0h, fluid) / 1e3, 4), round(cp.PropsSI("S", "P", p_c, "Q", 1, fluid) / 1e3, 2), round(cp.PropsSI("S", "P", p_c, "Q", 0, fluid) / 1e3, 2), round(cp.PropsSI("S", "P", p_0, "Q", 1, fluid) / 1e3, 2), round(cp.PropsSI("S", "P", p_0, "Q", 0, fluid) / 1e3, 2)],
        "Dampfqualität [%]": ["", "", "", 50, "", 100, 0, 100, 0],
    })

    df_results = pd.DataFrame({
        "Parameter": ["Projekt", "Leistungsaufnahme [kW]", "Wärmeleistung [kW]", "Kälteleistung [kW]", "Kältemittelmassenstrom [kg/s]", "Verdichtervolumenstrom [m3/h]", "COP", "EER"],
        "Wert": [project, round(P_com / 1e3, 2), round(Q_c / 1e3, 2), round(Q_0 / 1e3, 2), round(m_R, 6), round(V_com * 3.6e3, 2), round(COP, 2), round(EER, 2)]
    })

    df_pipes = pd.DataFrame({
        "Leitung": ["Heissgasleitung", "Flüssigkeitsleitung", "Saugleitung"],
        "Durchmesser [mm]": [d_hg, d_fl, d_sl],
        "Länge [m]": [l_hg, l_fl, l_sl],
        "Strömungsgeschwindigkeit [m/s]": [round(w_hg, 2), round(w_fl, 2), round(w_sl, 2)],
        "Jacobs-Geschwindigkeit [m/s]": [round(jg_hg, 2), round(jg_fl, 2), round(jg_sl, 2)],
        "Druckverlust [bar]": [round(dp_hg / 1e5, 2), round(dp_fl / 1e5, 2), round(dp_sl / 1e5, 2)],
        "Aussenmantelfläche [m2]": [A_M_hg, A_M_fl, A_M_sl],
        "Innenvolumen [dm3]": [V_i_hg, V_i_fl, V_i_sl],
    })

    return df, df_results, df_pipes, errors

title_col, ver_col = st.columns([6, 1])
with title_col:
    st.title(APP_TITLE)
with ver_col:
    st.markdown(f'<div style="padding-top:1.1rem;text-align:left;color:#666;">{APP_VERSION}</div>', unsafe_allow_html=True)

col1, col2 = st.columns([1.05, 1.0])

with col1:
    project = st.text_input("Projekt", "Projekt 1")
    fluid = st.selectbox("Kältemittel", list(FLUIDS.keys()))
    mode = st.selectbox("Eingabemodus", ["Kälteleistung", "Wärmeleistung", "Verdichtervolumenstrom"])

    show_q0 = mode == "Kälteleistung"
    show_qc = mode == "Wärmeleistung"
    show_vcom = mode == "Verdichtervolumenstrom"

    q0 = st.number_input("Kälteleistung [kW]", value=1.0, disabled=not show_q0) if show_q0 else 0.0
    qc = st.number_input("Wärmeleistung [kW]", value=1.0, disabled=not show_qc) if show_qc else 0.0
    vcom = st.number_input("Verdichtervolumenstrom [m3/h]", value=1.0, disabled=not show_vcom) if show_vcom else 0.0

    tc = st.number_input("Verflüssigungstemperatur [°C]", value=45.0)
    t0 = st.number_input("Verdampfungstemperatur [°C]", value=-10.0)
    dtcu = st.number_input("Unterkühlung [K]", value=5.0)
    dt0h = st.number_input("Verdampferüberhitzung [K]", value=10.0)
    dtsh = st.number_input("Saugleitungsüberhitzung [K]", value=10.0)
    ifpipes = st.selectbox("Rohrleitungsdimensionierung", ["Nein", "Ja"])

    show_pipe_inputs = ifpipes == "Ja"
    lhg = st.number_input("Heissgasleitungslänge [m]", value=25.0) if show_pipe_inputs else 25.0
    lfl = st.number_input("Flüssigkeitsleitungslänge [m]", value=17.8) if show_pipe_inputs else 17.8
    lsl = st.number_input("Saugleitungslänge [m]", value=11.0) if show_pipe_inputs else 11.0

    run = st.button("Berechnen")

with col2:
    st.subheader("Systemdaten")
    st.write("Hier stehen die wichtigsten Ergebniswerte.")
    st.subheader("Hinweise")
    st.write("Bei aktivierter Rohrleitungsdimensionierung wird die Rohrtabelle angezeigt.")

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
        st.subheader("Kreislaufpunkte")
        st.dataframe(df, use_container_width=True, hide_index=True)
        st.subheader("Systemdaten")
        st.dataframe(df_results, use_container_width=True, hide_index=True)
        if ifpipes == "Ja":
            st.subheader("Rohrleitungsdimensionierung")
            st.dataframe(df_pipes, use_container_width=True, hide_index=True)
        if errors:
            st.warning("\n".join(errors))
    except Exception as e:
        st.error(str(e))