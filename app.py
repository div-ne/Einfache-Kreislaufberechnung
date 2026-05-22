import numpy as np
import pandas as pd
import streamlit as st
import CoolProp.CoolProp as cp
from scipy.optimize import brentq

from Cu_pipes import cu_pipes

st.set_page_config(page_title="VCCSimple Streamlit", layout="wide")

APP_TITLE = "VCCSimple Streamlit"
APP_VERSION = "1.0"

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
    if not bigger:
        return None
    return min(bigger) * 1e-3

def find_next_smaller_pipe(di):
    di_values = [pipe[3] for pipe in cu_pipes]
    smaller = [x for x in di_values if x < di * 1e3]
    if not smaller:
        return None
    return max(smaller) * 1e-3

def Q_brentq_h(p, h, x_min, x_max, fluid):
    def f(x):
        return cp.PropsSI('H', 'P', p, 'Q', x, fluid) - h
    return brentq(f, x_min, x_max)

def T_brentq_s(p, s, T_min, T_max, fluid):
    def f(T):
        return cp.PropsSI('S', 'T', T, 'P', p, fluid) - s
    return brentq(f, T_min, T_max)

def T_brentq_h(p, h, T_min, T_max, fluid):
    def f(T):
        return cp.PropsSI('H', 'T', T, 'P', p, fluid) - h
    return brentq(f, T_min, T_max)

def run_calculation(inputs):
    project = inputs['project']
    mode = inputs['mode']
    errors = ''
    fluid = inputs['fluid']
    if fluid == 'CO2':
        fluid = 'CarbonDioxide'
    if mode == 0:
        Q0 = inputs['Q0'] * 1e3
        Qc = 0
        Vcom = 0
    elif mode == 1:
        Qc = inputs['Qc'] * 1e3
        Q0 = 0
        Vcom = 0
        if Qc > 0:
            Qc = -Qc
    else:
        Vcom = inputs['Vcom'] / 3.6e3
        Q0 = 0
        Qc = 0
    Tgco = inputs['Tgco'] + 273.15
    T0 = inputs['T0'] + 273.15
    dt0h = inputs['dt0h']
    dtsh = inputs['dtsh']
    if Tgco < T0:
        raise ValueError('Die Verflüssigungstemperatur darf nicht unter der Verdampfungstemperatur liegen.')
    if abs(Tgco - T0) < 20:
        errors += 'Die Differenz zwischen Verflüssigung und Verdampfung ist sehr gering.
'
    if Q0 == 0 and mode == 0:
        Q0 = 1
    if Qc == 0 and mode == 1:
        Qc = 1
    if Vcom == 0 and mode == 2:
        Vcom = 1
    if_pipes = inputs['ifpipes']
    if if_pipes == 0:
        lhg = lfl = lsl = 1
    else:
        lhg, lfl, lsl = inputs['lhg'], inputs['lfl'], inputs['lsl']
    p0dew = cp.PropsSI('P', 'T', T0, 'Q', 1, fluid)
    p0bubble = cp.PropsSI('P', 'T', T0, 'Q', 0, fluid)
    p0 = (p0dew + p0bubble) / 2
    Tgcodew = cp.PropsSI('T', 'P', p0, 'Q', 1, fluid)
    Tgcobubble = cp.PropsSI('T', 'P', p0, 'Q', 0, fluid)
    h0dew = cp.PropsSI('H', 'P', p0, 'Q', 1, fluid)
    h0bubble = cp.PropsSI('H', 'P', p0, 'Q', 0, fluid)
    T4 = cp.PropsSI('T', 'P', p0, 'Q', 0.5, fluid)
    h4 = cp.PropsSI('H', 'P', p0, 'Q', 0.5, fluid)
    T0h = T0 + dt0h
    Tsh = T0 + dt0h + dtsh
    pc = cp.PropsSI('P', 'T', Tgco, 'Q', 1, fluid)
    h3 = cp.PropsSI('H', 'T', Tgco, 'Q', 0, fluid)
    h1 = cp.PropsSI('H', 'P', p0, 'T', Tsh, fluid)
    s1 = cp.PropsSI('S', 'P', p0, 'T', Tsh, fluid)
    T2s = T_brentq_s(pc, s1, Tgco + 1, Tgco + 100, fluid)
    h2s = cp.PropsSI('H', 'P', pc, 'T', T2s, fluid)
    eta_is = 0.75
    h2 = h1 + (h2s - h1) / eta_is
    h5 = cp.PropsSI('H', 'P', p0, 'T', T0h, fluid)
    rho1 = cp.PropsSI('D', 'P', p0, 'T', Tsh, fluid)
    rho2 = cp.PropsSI('D', 'P', pc, 'T', Tgco, fluid)
    rho3 = cp.PropsSI('D', 'P', pc, 'T', Tgco - 1, fluid)
    mu1 = cp.PropsSI('V', 'P', p0, 'T', Tsh, fluid)
    mu2 = cp.PropsSI('V', 'P', pc, 'T', Tgco, fluid)
    mu3 = cp.PropsSI('V', 'P', pc, 'T', Tgco - 1, fluid)
    nu1, nu2, nu3 = mu1 / rho1, mu2 / rho2, mu3 / rho3
    lhg_aq, lfl_aq, lsl_aq = lhg * 1.4, lfl * 1.4, lsl * 1.4
    q_c = h3 - h2
    q_0 = h5 - h4
    if mode == 0:
        mR = Q0 / q_0
        Qc = mR * q_c
        Vcom = mR / cp.PropsSI('D', 'P', p0, 'T', Tsh + 10, fluid)
    elif mode == 1:
        mR = Qc / q_c
        Q0 = mR * q_0
        Vcom = mR / cp.PropsSI('D', 'P', p0, 'T', Tsh + 10, fluid)
    else:
        mR = Vcom * cp.PropsSI('D', 'P', p0, 'T', Tsh + 10, fluid)
        Qc = mR * q_c
        Q0 = mR * q_0
    Pcom = mR * (h2 - h1)
    COP = -Qc / Pcom
    EER = Q0 / Pcom
    dp_max_hg = float(pc - cp.PropsSI('P', 'T', Tgco - 1, 'Q', 1, fluid))
    dp_max_fl = dp_max_hg
    dp_max_sl = float(cp.PropsSI('P', 'T', T0 + 1, 'Q', 1, fluid) - p0)
    def project_pipe(name, pre_di, length, V, nu, rho, dp_max):
        return _project_pipe(name, pre_di, length, V, nu, rho, dp_max)
    def _project_pipe(name, pre_di, length, V, nu, rho, dp_max):
        di = pre_di
        for _ in range(30):
            w = 4 * V / (np.pi * di**2)
            lam = float(friction_coefficient(di, w, nu))
            dp = lam * length * rho * w**2 * 0.5 / di
            jg = np.sqrt(0.85 * 9.80665 * di * abs(rho * 1.1 - rho) / rho) if rho > 900 else np.sqrt(0.85 * 9.80665 * di * abs(900 - rho) / rho)
            if w > jg and dp < dp_max:
                return di, w, jg, dp
            if dp > dp_max:
                nd = find_next_bigger_pipe(di) or di
            else:
                nd = find_next_smaller_pipe(di) or di
            if nd == di:
                return di, w, jg, dp
            di = nd
        return di, w, jg, dp
    pre_hg = find_next_bigger_pipe(np.sqrt(4 * Vcom / np.pi / 15)) or 0.01
    pre_fl = find_next_bigger_pipe(np.sqrt(4 * Vcom / np.pi / 1)) or 0.01
    pre_sl = find_next_bigger_pipe(np.sqrt(4 * Vcom / np.pi / 12)) or 0.01
    di_hg_dp, w_hg_dp, jg_hg_dp, dp_hg_dp = _project_pipe('Heissgasleitung', pre_hg, lhg_aq, Vcom, nu2, rho2, dp_max_hg)
    di_hg_dm, w_hg_dm, jg_hg_dm, dp_hg_dm = _project_pipe('Heissgasleitung', pre_hg, lhg_aq, Vcom, nu2, rho2, 1e5)
    di_fl_dp, w_fl_dp, jg_fl_dp, dp_fl_dp = _project_pipe('Flüssigkeitsleitung', pre_fl, lfl_aq, Vcom, nu3, rho3, dp_max_fl)
    di_fl_dm, w_fl_dm, jg_fl_dm, dp_fl_dm = _project_pipe('Flüssigkeitsleitung', pre_fl, lfl_aq, Vcom, nu3, rho3, 1e5)
    di_sl_dp, w_sl_dp, jg_sl_dp, dp_sl_dp = _project_pipe('Saugleitung', pre_sl, lsl_aq, Vcom, nu1, rho1, dp_max_sl)
    di_sl_dm, w_sl_dm, jg_sl_dm, dp_sl_dm = _project_pipe('Saugleitung', pre_sl, lsl_aq, Vcom, nu1, rho1, 1e5)
    def pipe_row(leitung, variante, di, length, w, jg, dp):
        d = next(x[0] for x in cu_pipes if x[3] == round(di * 1e3, 1) or x[3] == round(di * 1e3))
        vi = next(x[4] for x in cu_pipes if x[3] == round(di * 1e3, 1) or x[3] == round(di * 1e3)) * length * 1e3
        am = next(x[5] for x in cu_pipes if x[3] == round(di * 1e3, 1) or x[3] == round(di * 1e3)) * length
        return [leitung, variante, d, round(length, 2), round(w, 2), round(jg, 2), round(dp / 1e5, 2), round(am, 3), round(vi, 2)]
    t1 = []
    t1.append(['1', round(Tsh - 273.15, 2), round(p0 / 1e5, 2), round(h1 / 1e3, 2), round(rho1, 2), round(cp.PropsSI('S', 'P', p0, 'T', Tsh, fluid) / 1e3, 4), 0])
    t1.append(['2', round(Tgco - 273.15, 2), round(pc / 1e5, 2), round(h2 / 1e3, 2), round(rho2, 2), round(cp.PropsSI('S', 'P', pc, 'T', Tgco, fluid) / 1e3, 4), 0])
    t1.append(['3', round(Tgco - 273.15, 2), round(pc / 1e5, 2), round(h3 / 1e3, 2), round(rho3, 2), round(cp.PropsSI('S', 'P', pc, 'T', Tgco - 1, fluid) / 1e3, 4), 0])
    t1.append(['4', round(T4 - 273.15, 2), round(p0 / 1e5, 2), round(h4 / 1e3, 2), round(cp.PropsSI('D', 'P', p0, 'Q', 0.5, fluid), 2), round(cp.PropsSI('S', 'P', p0, 'Q', 0.5, fluid) / 1e3, 4), 50])
    t1.append(['5', round(T0h - 273.15, 2), round(p0 / 1e5, 2), round(h5 / 1e3, 2), round(cp.PropsSI('D', 'P', p0, 'T', T0h, fluid), 2), round(cp.PropsSI('S', 'P', p0, 'T', T0h, fluid) / 1e3, 4), 0])
    df1 = pd.DataFrame(t1, columns=['P', 'Temperatur [°C]', 'Druck [bar]', 'Spezifische Enthalpie [kJ/kg]', 'Dichte [kg/m3]', 'Spezifische Entropie [kJ/kg/K]', 'Dampfqualität [%]'])
    df2 = pd.DataFrame([
        ['Leistung [kW]', round(Pcom / 1e3, 2)],
        ['Kälteleistung [kW]', round(Q0 / 1e3, 2)],
        ['Wärmeleistung [kW]', round(Qc / 1e3, 2)],
        ['Massenstrom [kg/s]', round(mR, 6)],
        ['Volumenstrom [m3/h]', round(Vcom * 3.6e3, 2)],
        ['COP', round(COP, 2)],
        ['EER', round(EER, 2)],
    ], columns=['Parameter', 'Wert'])
    df3 = pd.DataFrame([
        pipe_row('Heissgasleitung', 'Auf Druckverlust optimiert', di_hg_dp, lhg, w_hg_dp, jg_hg_dp, dp_hg_dp),
        pipe_row('Heissgasleitung', 'Auf Durchmesser optimiert', di_hg_dm, lhg, w_hg_dm, jg_hg_dm, dp_hg_dm),
        pipe_row('Flüssigkeitsleitung', 'Auf Druckverlust optimiert', di_fl_dp, lfl, w_fl_dp, jg_fl_dp, dp_fl_dp),
        pipe_row('Flüssigkeitsleitung', 'Auf Durchmesser optimiert', di_fl_dm, lfl, w_fl_dm, jg_fl_dm, dp_fl_dm),
        pipe_row('Saugleitung', 'Auf Druckverlust optimiert', di_sl_dp, lsl, w_sl_dp, jg_sl_dp, dp_sl_dp),
        pipe_row('Saugleitung', 'Auf Durchmesser optimiert', di_sl_dm, lsl, w_sl_dm, jg_sl_dm, dp_sl_dm),
    ], columns=['Leitung', 'Variante', 'Durchmesser [mm]', 'Länge [m]', 'Strömungsgeschwindigkeit [m/s]', 'Jacobs-Geschwindigkeit [m/s]', 'Druckverlust [bar]', 'Aussenmantelfläche [m2]', 'Innenvolumen [dm3]'])
    return df1, df2, df3, errors

st.title(APP_TITLE)
col1, col2 = st.columns([1, 1])
with col1:
    project = st.text_input('Projekt', 'Projekt 1')
    fluid = st.selectbox('Fluid', ['R290', 'R1270', 'R170', 'R1150', 'R744', 'R717', 'R718', 'R600', 'R600a', 'R601', 'R601a', 'R702', 'R729', 'R764', 'CO2'])
    mode = st.selectbox('Eingabemodus', ['Kälteleistung', 'Wärmeleistung', 'Verdichtervolumenstrom'])
    q0 = st.number_input('Kälteleistung [kW]', value=1.0)
    qc = st.number_input('Wärmeleistung [kW]', value=1.0)
    vcom = st.number_input('Verdichtervolumenstrom [m3/h]', value=1.0)
    tgc = st.number_input('Gaskühleraustrittstemperatur [°C]', value=45.0)
    t0 = st.number_input('Verdampfungstemperatur [°C]', value=-10.0)
    dt0h = st.number_input('Verdampferüberhitzung [K]', value=10.0)
    dtsh = st.number_input('Saugleitungsüberhitzung [K]', value=10.0)
    ifpipes = st.selectbox('Rohrleitungsdimensionierung', ['Nein', 'Ja'])
    lhg = st.number_input('Heissgasleitungslänge [m]', value=25.0)
    lfl = st.number_input('Flüssigkeitsleitungslänge [m]', value=17.8)
    lsl = st.number_input('Saugleitungslänge [m]', value=11.0)
    run = st.button('Berechnen')
with col2:
    st.markdown('### Hinweise')
    st.write('Bei aktivierter Rohrleitungsdimensionierung werden beide Varianten angezeigt: Druckverlust-optimiert und Durchmesser-optimiert.')
    st.write('Die Rohrdaten stammen aus `Cu_pipes.py`.')

if run:
    inputs = dict(project=project, fluid=fluid, mode=['Kälteleistung','Wärmeleistung','Verdichtervolumenstrom'].index(mode), Q0=q0, Qc=qc, Vcom=vcom, Tgco=tgc, T0=t0, dt0h=dt0h, dtsh=dtsh, ifpipes=0 if ifpipes=='Nein' else 1, lhg=lhg, lfl=lfl, lsl=lsl)
    df1, df2, df3, errors = run_calculation(inputs)
    st.subheader('Kreislaufpunkte')
    st.dataframe(df1, use_container_width=True)
    st.subheader('Systemdaten')
    st.dataframe(df2, use_container_width=True)
    st.subheader('Rohrleitungsdimensionierung')
    st.dataframe(df3, use_container_width=True)
    if errors:
        st.warning(errors)

st.caption(f'{APP_TITLE} – Version {APP_VERSION}')
