import numpy as np
import pandas as pd
import streamlit as st
import CoolProp.CoolProp as cp
from scipy.optimize import brentq

APPTITLE = "Kälteträger-Rechner"
APPVERSION = "0.9.1V"

ROUGHNESSROWS = [
    ["Gezogene und gepresste Rohre aus Kupfer, Messing, Bronze, Aluminium, Glas oder Kunststoff, neu, technisch glatt", "0.001"],
    ["Gezogene und gepresste Rohre aus Kupfer, Messing, Bronze, Aluminium, Glas oder Kunststoff, gebraucht", "0.010 - 0.0300"],
    ["Gummischlauch, neu, handelsüblich", "0.0016"],
    ["Rohre aus Gusseisen, neu, handelsüblich", "0.25 - 0.5"],
    ["Rohre aus Gusseisen, angerostet", "1.00 - 1.5"],
    ["Rohre aus Gusseisen, verkrustet", "1.50 - 3.0"],
    ["Rohre aus Gusseisen, nach mehrjährigem Betrieb gereinigt", "0.30 - 1.5"],
    ["Rohre aus Gusseisen, städtliche Kanalisation", "1.20"],
    ["Neue nahtlose Stahlrohre, gewalzt oder gezogen, mit Walzhaut", "0.02 - 0.06"],
    ["Neue nahtlose Stahlrohre, gewalzt oder gezogen, gebeizt", "0.03 - 0.04"],
    ["Neue nahtlose Stahlrohre, gewalzt oder gezogen, bei engen Rohren", "0.10"],
    ["Neue längsgeschweißte Stahlrohre, mit Walzhaut", "0.04 - 0.1"],
    ["Neue längsgeschweißte Stahlrohre, leicht verkrustet", "1.00 - 1.5"],
    ["Neue längsgeschweißte Stahlrohre, stark verkrustet", "2.00 - 4.0"],
    ["Neue längsgeschweißte Stahlrohre, gebraucht und gereinigt", "0.15 - 0.2"],
    ["Neue Stahlrohre mit Überzug, Metallspritzung", "0.08 - 0.09"],
    ["Neue Stahlrohre mit Überzug, tauchverzinkt", "0.07 - 0.10"],
    ["Neue Stahlrohre mit Überzug, handelsüblich verzinkt", "0.10 - 0.16"],
    ["Neue Stahlrohre mit Überzug, bituminiert", "0.050"],
    ["Neue Stahlrohre mit Überzug, zementiert", "0.180"],
    ["Neue Stahlrohre mit Überzug, galvanisiert", "0.008"],
    ["Gebrauchte Stahlrohre, gleichmäßige Rostnarben", "0.15"],
    ["Gebrauchte Stahlrohre, leichte Verkrustung", "0.15 - 0.4"],
    ["Gebrauchte Stahlrohre, mittlere Verkrustung", "1.50"],
    ["Gebrauchte Stahlrohre, starke Verkrustung", "2.00 - 4.0"],
    ["Asbest-Zementrohre, neu, handelsüblich", "0.03 - 0.1"],
    ["Betonrohre, Druckstollen, handelsüblich glattstrich", "0.3 - 0.8"],
    ["Betonrohre, Druckstollen, handelsüblich mittelglatt", "1.0 - 2.0"],
    ["Betonrohre, Druckstollen, handelsüblich rau", "2.0 - 3.0"],
    ["Betonrohre, Druckstollen, mehrjähriger Betrieb mit Wasser", "0.2 - 0.3"],
    ["Neues Tonrohr, Drainagerohr, gebrannt", "0.6 - 0.8"],
    ["Neues Tonrohr, aus rohen Tonziegeln", "9.0"],
    ["Medizinisches, Kälte- oder Heizungsgewinderohr, neu, handelsüblich", "0.045"],
    ["Medizinisches, Kälte- oder Heizungsstahlrohr nahtlos, neu, handelsüblich", "0.045"],
    ["Medizinisches, Kälte- oder Heizungskupferrohr, neu, handelsüblich", "0.0005 - 0.0015"],
    ["Medizinisches, Kälte- oder Heizungspräzisionsstahlrohr, neu, handelsüblich", "0.001 - 0.0015"],
    ["Medizinisches, Kälte- oder Heizungskunststoffrohr, neu, handelsüblich", "0.001 - 0.0015"],
]

LIMITSROWS = [
    ["Wasser", "100", "ohne Konzentrationsauswahl"],
    ["Antifrogen N", "10 bis 60", "zulässiger Eingabebereich"],
    ["Antifrogen L", "10 bis 60", "zulässiger Eingabebereich"],
    ["Kaliumformiat", "40 bis 100", "zulässiger Eingabebereich"],
]

CUPIPES = [
    (6, 1, 4), (8, 1, 6), (10, 1, 8), (12, 1, 10), (15, 1, 13), (16, 1, 14),
    (18, 1, 16), (22, 1, 20), (28, 1, 26), (28, 1.5, 25), (35, 1.5, 32),
    (42, 1.5, 39), (54, 1.5, 51), (54, 2, 50), (64, 2, 60), (76.1, 2, 74.1),
    (88.9, 2, 84.9), (108, 2.5, 103), (133, 3, 127)
]


def buildfluid(name, concentration):
    if name == "Wasser":
        return "Water", "Wasser", 100
    if name == "Antifrogen N":
        return f"INCOMP::AN{concentration/100:.4f}", f"Antifrogen N {concentration:.0f} %", concentration
    if name == "Antifrogen L":
        return f"INCOMP::AL{concentration/100:.4f}", f"Antifrogen L {concentration:.0f} %", concentration
    if name == "Kaliumformiat":
        return f"INCOMP::KF{concentration/100:.4f}", f"Kaliumformiat {concentration:.0f} %", concentration
    raise ValueError(f"Unbekanntes Fluid: {name}")


def safefreezetemp(coolantname, fluid, atmosphericpressure):
    if coolantname == "Wasser":
        return 273.15
    try:
        return cp.PropsSI("Tfreeze", "T", 0, "P", atmosphericpressure, fluid)
    except Exception:
        return np.nan


def frictioncoefficient(di, w, nu, k):
    Re = w * di / nu
    epsilonk = k / di
    if Re <= 2320:
        return 64 / Re
    lambdablasius = 0.3164 * Re ** -0.25
    lambdanikuradse = (-2 * np.log10(k / (3.71 * di))) ** -2
    lambdaprandtl = 0.02
    for _ in range(10):
        lambdaprandtl = (-2 * np.log10(2.51 / (Re * np.sqrt(lambdaprandtl)))) ** -2
    lambdacolebrookwhite = 0.02
    for _ in range(10):
        lambdacolebrookwhite = (-2 * np.log10((2.51 / (Re * np.sqrt(lambdacolebrookwhite))) + (k / (3.71 * di)))) ** -2
    check = Re * np.sqrt(lambdanikuradse) * k / di
    if check < 200:
        return lambdanikuradse
    if epsilonk <= 0.001 and Re <= 10000:
        return lambdablasius
    if epsilonk <= 0.0002 and Re <= 100000:
        return lambdablasius
    if epsilonk <= 0.00002 and Re <= 1000000:
        return lambdaprandtl
    if epsilonk <= 0.00001:
        return lambdaprandtl
    return lambdacolebrookwhite


def next_pipe(di, bigger=True):
    vals = [p[0] for p in CUPIPES]
    if bigger:
        cands = [x for x in vals if x >= di * 1e3]
        return (min(cands) / 1e3) if cands else None
    cands = [x for x in vals if x <= di * 1e3]
    return (max(cands) / 1e3) if cands else None


def pipe_area(da):
    return np.pi * (da * 1e-3) ** 2 / 4


def pipe_vi(di):
    return np.pi * (di * 1e-3) ** 2 / 4


def pipe_am(da):
    return np.pi * da * 1e-3


def calc_pipe_name(di):
    for da, t, dii in CUPIPES:
        if abs(dii - di * 1e3) < 1e-6:
            return f"{da}x{t}"
    return f"{di*1e3:.1f}"


def findoptpipe(name, wpre, dipre, laq, dpmax, V, nu, rho, errors, pipeopt):
    di = dipre
    nsmallcheck = 0
    nbigcheck = 0
    nmax = 15
    while True:
        w = 4 * V / (np.pi * di ** 2)
        lambdarohr = float(frictioncoefficient(di, w, nu, 0.0015))
        dptat = lambdarohr * laq * rho * w ** 2 / (2 * di)
        if rho < 900:
            jg = np.sqrt(0.85 * 9.80665 * di * (rho * 1.1 - rho) / rho)
        else:
            jg = np.sqrt(0.85 * 9.80665 * di * (900 - rho) / rho)
        if nsmallcheck >= nmax or nbigcheck >= nmax:
            di = dipre
            w = 4 * V / (np.pi * di ** 2)
            errors = f"{name} konnte nicht eindeutig projektiert werden. Es wurde eine vorläufige Auswahl für {wpre} m/s getroffen."
            return di, w, jg, dptat, errors
        if w <= jg and dptat <= dpmax:
            return di, w, jg, dptat, errors
        if pipeopt == 0:
            if w > jg or dptat > dpmax:
                newdi = next_pipe(di, bigger=True)
                if newdi is None:
                    return di, w, jg, dptat, errors
                di = newdi
                nbigcheck += 1
            else:
                newdi = next_pipe(di, bigger=False)
                if newdi is None:
                    return di, w, jg, dptat, errors
                di = newdi
                nsmallcheck += 1
        else:
            if w > jg or dptat > dpmax:
                newdi = next_pipe(di, bigger=False)
                if newdi is None:
                    return di, w, jg, dptat, errors
                di = newdi
                nsmallcheck += 1
            else:
                newdi = next_pipe(di, bigger=True)
                if newdi is None:
                    return di, w, jg, dptat, errors
                di = newdi
                nbigcheck += 1


def runcalculationinputs(inputs):
    project = inputs.get("project", "")
    mode = inputs["mode"]
    fluidinput = inputs["fluidinput"]
    concentration = inputs["concentration"]
    Tc = float(inputs["Tc"]) + 273.15
    T0 = float(inputs["T0"]) + 273.15
    dtcu = float(inputs["dtcu"])
    dt0h = float(inputs["dt0h"])
    dtsh = float(inputs["dtsh"])
    ifpipes = inputs["ifpipes"]
    pipeopt = inputs["pipeopt"]
    lhg = float(inputs.get("lhg", 1))
    lfl = float(inputs.get("lfl", 1))
    lsl = float(inputs.get("lsl", 1))

    errors = "Keine Fehlermeldungen."
    atmosphericpressure = 1.01325e5

    fluid, name, concentration = buildfluid(fluidinput, concentration)
    if fluidinput == "Wasser":
        concentration = 100
    if Tc < T0:
        raise ValueError("Die Verflüssigungstemperatur darf nicht unter der Verdampfungstemperatur liegen.")
    if Tc - T0 < 20:
        st.warning("Die Differenz zwischen Verflüssigungs- und Verdampfungstemperatur ist sehr gering. Dies kann zu fehlerhaften Berechnungen führen.")

    if mode == "Kälteleistung":
        Q0 = float(inputs.get("Q0", 0)) * 1e3
        Qc = 0.0
        Vcom = 0.0
    elif mode == "Wärmeleistung":
        Qc = float(inputs.get("Qc", 0)) * 1e3
        Q0 = 0.0
        Vcom = 0.0
    else:
        Vcom = float(inputs.get("Vcom", 0)) / 3.6e3
        Q0 = 0.0
        Qc = 0.0

    if Q0 == 0 and mode == "Kälteleistung":
        Q0 = 1
    if Qc == 0 and mode == "Wärmeleistung":
        Qc = 1
    if Vcom == 0 and mode == "Verdichtervolumenstrom":
        Vcom = 1

    pc = cp.PropsSI("P", "T", Tc, "Q", 1, fluid)
    p0 = cp.PropsSI("P", "T", T0, "Q", 1, fluid)
    Tcbubble = cp.PropsSI("T", "P", pc, "Q", 0, fluid)
    Tcdew = cp.PropsSI("T", "P", pc, "Q", 1, fluid)
    T0bubble = cp.PropsSI("T", "P", p0, "Q", 0, fluid)
    T0dew = cp.PropsSI("T", "P", p0, "Q", 1, fluid)
    Tcu = Tcbubble - dtcu
    Tsh = T0dew + dtsh
    T4 = T0dew - dt0h
    T0h = T0dew + dt0h

    h1 = cp.PropsSI("H", "P", p0, "T", Tsh, fluid)
    s1 = cp.PropsSI("S", "P", p0, "T", Tsh, fluid)
    h2s = cp.PropsSI("H", "P", pc, "T", cp.PropsSI("T", "P", pc, "S", s1, fluid), fluid)
    pi = pc / p0
    etais = 0.821 - 0.0105 * pi if pi <= 5 else -0.020644445 + 0.68403852 * pi - 0.22147167 * pi**2 + 0.032145926 * pi**3 - 0.00178 * pi**4
    ws = h2s - h1
    wr = ws / etais
    h2 = h1 + wr
    h3 = cp.PropsSI("H", "T", Tcbubble, "Q", 0, fluid)
    h4 = cp.PropsSI("H", "P", p0, "Q", 0, fluid)
    h5 = cp.PropsSI("H", "P", p0, "T", T0h, fluid)
    rho1 = cp.PropsSI("D", "P", p0, "T", Tsh, fluid)
    rho2 = cp.PropsSI("D", "P", pc, "T", cp.PropsSI("T", "P", pc, "S", s1, fluid), fluid)
    rho3 = cp.PropsSI("D", "P", pc, "T", Tcu, fluid)
    rho0f = cp.PropsSI("D", "P", p0, "Q", 0, fluid)
    mu1 = cp.PropsSI("V", "P", p0, "T", Tsh, fluid)
    mu2 = cp.PropsSI("V", "P", pc, "T", cp.PropsSI("T", "P", pc, "S", s1, fluid), fluid)
    mu3 = cp.PropsSI("V", "P", pc, "T", Tcu, fluid)
    nu1 = mu1 / rho1
    nu2 = mu2 / rho2
    nu3 = mu3 / rho3
    q0 = h1 - h5
    qc = h3 - h2
    if mode == "Kälteleistung":
        mR = Q0 / q0
        Pcom = mR * wr
        Qc = mR * qc
        Vcom = mR / rho1
    elif mode == "Wärmeleistung":
        mR = Qc / qc
        Pcom = mR * wr
        Q0 = mR * q0
        Vcom = mR / rho1
    else:
        mR = Vcom * rho1
        Q0 = mR * q0
        Qc = mR * qc
        Pcom = mR * wr

    COP = Q0 / Pcom
    EER = Qc / Pcom
    COPcarnot = Tc / (Tc - T0)
    EERcarnot = T0 / (Tc - T0)
    etaCOP = COP / COPcarnot
    etaEER = EER / EERcarnot

    laqhg = lhg * 1.4
    laqfl = lfl * 1.4
    laqsl = lsl * 1.4

    rho4 = cp.PropsSI("D", "P", p0, "Q", 1, fluid)
    rho5 = cp.PropsSI("D", "P", p0, "T", T0h, fluid)
    mu4 = cp.PropsSI("V", "P", p0, "Q", 1, fluid)
    mu5 = cp.PropsSI("V", "P", p0, "T", T0h, fluid)
    nu4 = mu4 / rho4
    nu5 = mu5 / rho5

    dpmaxhg = float(pc - cp.PropsSI("P", "T", Tcdew - 1, "Q", 1, fluid))
    dpmaxfl = float(pc - cp.PropsSI("P", "T", Tcbubble - 1, "Q", 1, fluid))
    dpmaxsl = float(cp.PropsSI("P", "T", T0dew - 1, "Q", 1, fluid) - p0)

    whgpre = 15
    wflpre = 1
    wslpre = 12

    dihgpre = np.sqrt(4 * Vcom / (np.pi * whgpre))
    diflpre = np.sqrt(4 * Vcom / (np.pi * wflpre))
    dislpre = np.sqrt(4 * Vcom / (np.pi * wslpre))

    dihg, whg, jghg, dphg, errors = findoptpipe("Heissgasleitung", whgpre, dihgpre, laqhg, dpmaxhg, Vcom, nu2, rho2, errors, pipeopt)
    difl, wfl, jgfl, dpfl, errors = findoptpipe("Flüssigkeitsleitung", wflpre, diflpre, laqfl, dpmaxfl, Vcom, nu3, rho3, errors, pipeopt)
    disl, wsl, jgsl, dpsl, errors = findoptpipe("Saugleitung", wslpre, dislpre, laqsl, dpmaxsl, Vcom, nu1, rho1, errors, pipeopt)

    pipe_rows = []
    if ifpipes == 1:
        for label, di, l, w, jg, dp in [
            ("Heissgasleitung", dihg, lhg, whg, jghg, dphg),
            ("Flüssigkeitsleitung", difl, lfl, wfl, jgfl, dpfl),
            ("Saugleitung", disl, lsl, wsl, jgsl, dpsl),
        ]:
            da = next((x[0] for x in CUPIPES if abs(x[2] - di * 1e3) < 1e-6), di * 1e3)
            pipe_rows.append([label, da, l, w, jg, dp / 1e5, pipe_vi(di) * 1e3, pipe_am(da)])

    df = pd.DataFrame([
        ["Bezeichnung", name, name],
        ["Übertragene Leistung [kW]", round(Q0 / 1000, 2), round(Qc / 1000, 2)],
        ["Volumenstrom [m3/h]", round(Vcom * 3.6e3, 2), round(Vcom * 3.6e3 * rho1 / rho3, 2)],
        ["Strömungsgeschwindigkeit [m/s]", round(whg, 2), round(wfl, 2)],
        ["Druckverlust [bar]", round(dphg / 1e5, 2), round(dpfl / 1e5, 2)],
        ["Dichte [kg/m3]", round(rho1, 2), round(rho3, 2)],
        ["Spezifische Wärmekapazität [kJ/kg/K]", round(cp.PropsSI("C", "P", p0, "T", Tsh, fluid) / 1000, 2), round(cp.PropsSI("C", "P", pc, "T", Tcu, fluid) / 1000, 2)],
        ["Mittlere Temperatur [°C]", round(inputs["Tc"], 2), round(inputs["T0"], 2)],
        ["Gefrierpunkt [°C]", round(safefreezetemp(fluidinput, fluid, atmosphericpressure) - 273.15, 2), round(safefreezetemp(fluidinput, fluid, atmosphericpressure) - 273.15, 2)],
    ], columns=["Parameter", "Kälteträger 1", "Kälteträger 2"])

    dfresults = pd.DataFrame([
        ["Leistungsaufnahme [kW]", round(Pcom / 1000, 2)],
        ["Wärmeleistung [kW]", round(Qc / 1000, 2)],
        ["Kälteleistung [kW]", round(Q0 / 1000, 2)],
        ["Kältemittelmassenstrom [kg/s]", round(mR, 6)],
        ["Verdichtervolumenstrom [m3/h]", round(Vcom * 3.6e3, 2)],
        ["COP [-]", round(COP, 2)],
        ["EER [-]", round(EER, 2)],
        ["Gütegrad COP [-]", round(etaCOP * 100, 2)],
        ["Gütegrad EER [-]", round(etaEER * 100, 2)],
        ["Fehlermeldungen", errors],
    ], columns=["Parameter", "Wert"])

    dfpipes = pd.DataFrame(pipe_rows, columns=["Leitung", "Durchmesser [mm]", "Länge [m]", "Strömungsgeschwindigkeit [m/s]", "Jacobs-Geschwindigkeit [m/s]", "Druckverlust [bar]", "Innenvolumen [dm3]", "Aussenmantelfläche [m2]"])
    return df, dfresults, dfpipes

st.set_page_config(page_title=APPTITLE, layout="wide")
st.markdown(f"""
<div style='display:flex; align-items:baseline; gap:14px; flex-wrap:wrap; margin-bottom:0.2rem'>
  <div style='font-size:3rem; font-weight:700; line-height:1.1'>{APPTITLE}</div>
  <div style='color:#9ca3af; font-size:1rem; line-height:1.1'>{APPVERSION}</div>
</div>
""", unsafe_allow_html=True)
st.caption("Rechnet thermohydraulische Kennwerte von einem Kälteträger auf einen anderen um.")

left, right = st.columns([1, 1.25])
with left:
    project = st.text_input("Projekt", value="Projekt")
    c1col, c1pctcol = st.columns(2)
    with c1col:
        coolant1 = st.selectbox("Fluid 1", ["Wasser", "Antifrogen N", "Antifrogen L", "Kaliumformiat"], index=0)
    with c1pctcol:
        concentration1 = st.number_input("Konzentration Fluid 1 [%]", value=100.0, step=1.0, disabled=coolant1 == "Wasser")
        if coolant1 == "Wasser":
            concentration1 = 100.0
    vqcol, vqvalcol = st.columns(2)
    with vqcol:
        VorQ = st.selectbox("Basis", ["Volumenstrom Fluid 1 [m3/h]", "Übertragene Leistung Fluid 1 [kW]", "Verdichtervolumenstrom [m3/h]"], index=0)
    with vqvalcol:
        vqvalue = st.number_input("Wert", value=10.0, step=0.1)
    pressuredropcoolant1 = st.number_input("Druckverlust Fluid 1 [bar]", value=2.5, step=0.1)
    tcol, dtcol = st.columns(2)
    with tcol:
        meantemperature = st.number_input("Mittlere Temperatur [°C]", value=7.5, step=0.1)
    with dtcol:
        temperaturedifference = st.number_input("Temperaturdifferenz [K]", value=5.0, step=0.1)
    dicol, roughcol = st.columns(2)
    with dicol:
        meaninnerpipediameter = st.number_input("Rohrinnendurchmesser [mm]", value=25.0, step=0.1)
    with roughcol:
        roughness = st.number_input("Rohrrauheit [mm]", value=0.0015, step=0.0001, format="%.4f")
    c2col, c2pctcol = st.columns(2)
    with c2col:
        coolant2 = st.selectbox("Fluid 2", ["Wasser", "Antifrogen N", "Antifrogen L", "Kaliumformiat"], index=1)
    with c2pctcol:
        concentration2 = st.number_input("Konzentration Fluid 2 [%]", value=34.0, step=1.0, disabled=coolant2 == "Wasser")
        if coolant2 == "Wasser":
            concentration2 = 100.0
    run = st.button("Berechnen", use_container_width=True)

with right:
    st.subheader("Ergebnis")
    if run:
        try:
            mode = "Kälteleistung" if VorQ == "Übertragene Leistung Fluid 1 [kW]" else ("Verdichtervolumenstrom" if VorQ == "Verdichtervolumenstrom [m3/h]" else "Volumenstrom")
            inputs = {
                "project": project,
                "mode": mode,
                "fluidinput": coolant1,
                "concentration": float(concentration1),
                "Q0": float(vqvalue) if mode == "Kälteleistung" else 0,
                "Qc": float(vqvalue) if mode == "Volumenstrom" else 0,
                "Vcom": float(vqvalue) if mode == "Verdichtervolumenstrom" else 0,
                "Tc": float(meantemperature),
                "T0": float(meantemperature - temperaturedifference),
                "dtcu": float(temperaturedifference),
                "dt0h": 10.0,
                "dtsh": 10.0,
                "ifpipes": 1,
                "pipeopt": 0,
                "lhg": 25.0,
                "lfl": 17.8,
                "lsl": 11.0,
            }
            df, dfresults, dfpipes = runcalculationinputs(inputs)
            st.dataframe(df, use_container_width=True, hide_index=True)
            st.download_button("CSV herunterladen", data=df.to_csv(index=False).encode("utf-8"), file_name="kaeltetraeger-rechner-ergebnis.csv", mime="text/csv", use_container_width=True)
        except Exception as e:
            st.error(f"Fehler bei der Berechnung: {e}")
    else:
        st.info("Eingaben setzen und auf Berechnen klicken.")

st.markdown("---")
with st.expander("Anleitung"):
    st.markdown("""
Mit diesem Tool kannst du von einem Kälteträger 1 auf einen Kälteträger 2 umrechnen.

Für die Berechnung wählst du zuerst beim Kälteträger 1 entweder den Volumenstrom oder die übertragene Leistung als Ausgangsbasis.
Danach gibst du den Druckverlust mit Kälteträger 1 an. Zusätzlich werden allgemeine Angaben benötigt: mittlere Temperatur, Temperaturdifferenz, Rohrinnendurchmesser, Rohrrauheit sowie Kälteträger 2 mit seiner Konzentration.

Ein wesentlicher Nutzen des Rechners ist, dass neben den hydraulischen Kenngrößen auch thermodynamische Eigenschaften der gewählten Kälteträger berechnet und gegenübergestellt werden.
Dazu gehören insbesondere Dichte, spezifische Wärmekapazität und Gefrierpunkt.
""")
with st.expander("Rohrrauheitswerte"):
    st.caption("Orientierungswerte für die Eingabe von k in mm.")
    st.dataframe(pd.DataFrame(ROUGHNESSROWS, columns=["Werkstoff und Rohrart", "k [mm]"],), use_container_width=True, hide_index=True)
with st.expander("Grenzen der Kälteträger"):
    st.caption("Zulässige Konzentrationsbereiche für die Eingabe in diesem Tool.")
    st.dataframe(pd.DataFrame(LIMITSROWS, columns=["Kälteträger", "Mögliche Konzentration", "Hinweis"]), use_container_width=True, hide_index=True)
