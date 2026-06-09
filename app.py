import os
import numpy as np
import pandas as pd
import streamlit as st
import CoolProp.CoolProp as cp

APP_TITLE = "Rohrreibungszahl-Bestimmung"
APP_VERSION = "v1.0.0"
MOODY_SOURCE = "https://en.wikipedia.org/wiki/Moody_chart"
MOODY_IMAGE = os.path.join(os.path.dirname(__file__), "assets", "moody-diagram.jpg")

ROUGHNESS_ROWS = [
    ("Gezogene und gepresste Rohre aus Kupfer, Messing, Bronze, Aluminium, Glas oder Kunststoff", "neu, technisch glatt", "0.001 … 0.0015"),
    ("", "gebraucht", "0.010 … 0.0300"),
    ("Gummischlauch", "neu, handelsüblich", "0.0016"),
    ("Rohre aus Gusseisen", "neu, handelsüblich", "0.25 … 0.5"),
    ("", "angerostet", "1.00 … 1.5"),
    ("", "verkrustet", "1.50 … 3.0"),
    ("", "nach mehrjährigem Betrieb gereinigt", "0.30 … 1.5"),
    ("", "städtliche Kanalisation", "1.20"),
    ("Neue nahtlose Stahlrohre, gewalzt oder gezogen", "mit Walzhaut", "0.02 … 0.06"),
    ("", "gebeizt", "0.03 … 0.04"),
    ("", "bei engen Rohren", "… 0.10"),
    ("Neue längsgeschweisste Stahlrohre", "mit Walzhaut", "0.04 … 0.1"),
    ("", "leicht verkrustet", "1.00 … 1.5"),
    ("", "stark verkrustet", "2.00 … 4.0"),
    ("", "gebraucht und gereinigt", "0.15 … 0.2"),
    ("Neue Stahlrohre mit Überzug", "Metallspritzung", "0.08 … 0.09"),
    ("", "tauchverzinkt", "0.07 … 0.10"),
    ("", "handelsüblich verzinkt", "0.10 … 0.16"),
    ("", "bituminiert", "0.050"),
    ("", "zementiert", "0.180"),
    ("", "galvanisiert", "0.008"),
    ("Gebrauchte Stahlrohre", "gleichmässige Rostnarben", "0.15"),
    ("", "leichte Verkrustung", "0.15 … 0.4"),
    ("", "mittlere Verkrustung", "1.50"),
    ("", "starke Verkrustung", "2.00 … 4.0"),
    ("Asbest-Zementrohre", "neu, handelsüblich", "0.03 … 0.1"),
    ("Betonrohre, Druckstollen", "handelsüblich Glattstrich", "0.3 … 0.8"),
    ("", "handelsüblich mittelglatt", "1.0 … 2.0"),
    ("", "handelsüblich rau", "2.0 … 3.0"),
    ("", "mehrjähriger Betrieb mit Wasser", "0.2 … 0.3"),
    ("Neues Tonrohr", "Drainagerohr, gebrannt", "0.6 … 0.8"),
    ("", "aus rohen Tonziegeln", "9.0"),
    ("Medizinisches, Kälte- oder Heizungsgewinderohr", "neu, handelsüblich", "0.045"),
    ("Medizinisches, Kälte- oder Heizungsstahlrohr nahtlos", "neu, handelsüblich", "0.045"),
    ("Medizinisches, Kälte- oder Heizungskupferrohr", "neu, handelsüblich", "0.0005 … 0.0015"),
    ("Medizinisches, Kälte- oder Heizungspräzisionsstahlrohr", "neu, handelsüblich", "0.001 … 0.0015"),
    ("Medizinisches, Kälte- oder Heizungskunststoffrohr", "neu, handelsüblich", "0.001 … 0.0015"),
]

FLUIDS = {
    "Wasser": "H2O",
    "R11": "R11",
    "R12": "R12",
    "R22": "R22",
    "R32": "R32",
    "R123": "R123",
    "R124": "R124",
    "R125": "R125",
    "R142b": "R142b",
    "R143a": "R143a",
    "R152A": "R152A",
    "R218": "R218",
    "R227EA": "R227EA",
    "R236EA": "R236EA",
    "R236FA": "R236FA",
    "R245fa": "R245fa",
    "R290": "Propane",
    "R404A": "R404A",
    "R407A": "R407A.mix",
    "R407B": "R407B.mix",
    "R407C": "R407C",
    "R410A": "R410A",
    "R411A": "R411A.mix",
    "R411B": "R411B.mix",
    "R415A": "R415A.mix",
    "R415B": "R415B.mix",
    "R417A": "R417A.mix",
    "R417B": "R417B.mix",
    "R417C": "R417C.mix",
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
    "R425A": "R425A.mix",
    "R427A": "R427A.mix",
    "R428A": "R428A.mix",
    "R430A": "R430A.mix",
    "R431A": "R431A.mix",
    "R432A": "R432A.mix",
    "R433A": "R433A.mix",
    "R433B": "R433B.mix",
    "R433C": "R433C.mix",
    "R434A": "R434A.mix",
    "R436A": "R436A.mix",
    "R436B": "R436B.mix",
    "R439A": "R439A.mix",
    "R440A": "R440A.mix",
    "R441A": "R441A.mix",
    "R442A": "R442A.mix",
    "R443A": "R443A.mix",
    "R444A": "R444A.mix",
    "R444B": "R444B.mix",
    "R449A": "R449A.mix",
    "R449B": "R449B.mix",
    "R451A": "R451A.mix",
    "R451B": "R451B.mix",
    "R452A": "R452A.mix",
    "R454A": "R454A.mix",
    "R454B": "R454B.mix",
    "R500": "R500.mix",
    "R501": "R501.mix",
    "R507A": "R507A",
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
    "R744": "CO2",
    "R1233zd(E)": "R1233zd(E)",
    "R1234yf": "R1234yf",
    "R1234ze(E)": "R1234ze(E)",
}

st.set_page_config(page_title=APP_TITLE, layout="wide", page_icon="logo.png")

if "run_calculation" not in st.session_state:
    st.session_state["run_calculation"] = False


def calculate_lambda(fluid, temperature_c, pressure_bar, diameter_mm, velocity_ms, roughness_mm):
    k = roughness_mm * 1e-3
    di = diameter_mm * 1e-3
    p = pressure_bar * 1e5
    t = 273.15 + temperature_c
    nu = cp.PropsSI('V', 'T', t, 'P', p, fluid) / cp.PropsSI('D', 'T', t, 'P', p, fluid)
    re = di * velocity_ms / nu
    epsilon_k = k / di
    lambda_hagenpoiseulle = 64 / re
    lambda_blasius = 0.3164 / re**0.25
    lambda_nikuradse = (-2 * np.log10(k / 3.71 / di))**-2
    lambda_prandtl = 0.02
    for _ in range(10):
        lambda_prandtl = (2 * np.log10(re * np.sqrt(lambda_prandtl)))**-2
    lambda_colebrookwhite = 0.02
    for _ in range(10):
        lambda_colebrookwhite = (-2 * np.log10(2.51 / re / lambda_colebrookwhite + k / 3.71 / di))**-2
    check_moody_diagram = re * np.sqrt(lambda_nikuradse) * k / di
    if re < 2320:
        lambda_result = lambda_hagenpoiseulle
        selected = "Hagen-Poiseuille"
        flow_type = "laminare Strömung"
    elif check_moody_diagram >= 200:
        lambda_result = lambda_nikuradse
        selected = "Nikuradse"
        flow_type = "hydraulisch raues Rohr"
    elif epsilon_k < 0.001 and re < 10000:
        lambda_result = lambda_blasius
        selected = "Blasius"
        flow_type = "hydraulisch glattes Rohr, Re < 100000"
    elif epsilon_k < 0.0002 and re < 100000:
        lambda_result = lambda_blasius
        selected = "Blasius"
        flow_type = "hydraulisch glattes Rohr, Re < 100000"
    elif epsilon_k < 0.00002 and re < 1000000:
        lambda_result = lambda_prandtl
        selected = "Prandtl"
        flow_type = "hydraulisch glattes Rohr, Re > 100000"
    elif epsilon_k < 0.00001:
        lambda_result = lambda_prandtl
        selected = "Prandtl"
        flow_type = "hydraulisch glattes Rohr, Re > 100000"
    else:
        lambda_result = lambda_colebrookwhite
        selected = "Colebrook-White"
        flow_type = "Übergangsbereich zwischen hydraulisch glatt und hydraulisch rau"
    overview = pd.DataFrame(
        [
            ("Fluid", fluid),
            ("Temperatur [°C]", round(temperature_c, 3)),
            ("Druck [bar]", round(pressure_bar, 5)),
            ("Innendurchmesser [mm]", round(diameter_mm, 4)),
            ("Strömungsgeschwindigkeit [m/s]", round(velocity_ms, 4)),
            ("Rohrauheit [mm]", round(roughness_mm, 6)),
            ("Kinematische Viskosität ν [m²/s]", nu),
            ("Reynolds-Zahl Re [-]", f"{re:.0f}"),
            ("Relative Rauheit k/d [-]", f"{epsilon_k:.2e}"),
            ("Berechnungsgrundlage", selected),
            ("Strömungsbedingung", flow_type),
            ("Rohrreibungszahl λ [-]", lambda_result),
        ],
        columns=["Parameter", "Wert"],
    )
    formulas = pd.DataFrame(
        [
            ("Hagen-Poiseuille", lambda_hagenpoiseulle, "laminare Strömung"),
            ("Blasius", lambda_blasius, "hydraulisch glattes Rohr, Re < 100000"),
            ("Prandtl", lambda_prandtl, "hydraulisch glattes Rohr, Re > 100000"),
            ("Colebrook-White", lambda_colebrookwhite, "Übergangsbereich zwischen hydraulisch glatt und hydraulisch rau"),
            ("Nikuradse", lambda_nikuradse, "hydraulisch raues Rohr"),
        ],
        columns=["Formel", "Rohrreibungszahl λ [-]", "Zuordnung"],
    )
    return overview, formulas


st.markdown(
    f"""
    <div style='display:flex; align-items:baseline; gap:14px; flex-wrap:wrap; margin-bottom:0.2rem;'>
        <div style='font-size:3rem; font-weight:700; line-height:1.1;'>{APP_TITLE}</div>
        <div style='color:#9ca3af; font-size:1rem; line-height:1.1;'>{APP_VERSION}</div>
    </div>
    """,
    unsafe_allow_html=True,
)

st.caption("Berechnung der Rohrreibungszahl auf Basis von Fluid, Temperatur, Druck, Geometrie und Strömungsgeschwindigkeit.")

left, right = st.columns([1, 1.2])

with left:
    st.subheader("Eingabe")
    fluid_label = st.selectbox("Fluid", list(FLUIDS.keys()), index=list(FLUIDS.keys()).index("Wasser"))
    fluid = FLUIDS[fluid_label]
    temperature_c = st.number_input("Temperatur [°C]", value=10.0, step=0.1)
    pressure_bar = st.number_input("Druck [bar]", value=1.01325, step=0.00001, format="%.5f")
    diameter_mm = st.number_input("Innendurchmesser [mm]", value=50.0, step=0.1)
    velocity_ms = st.number_input("Strömungsgeschwindigkeit [m/s]", value=2.0, step=0.1)
    roughness_mm = st.number_input("Rohrrauheit [mm]", value=0.0015, step=0.0001, format="%.4f")
    run = st.button("Berechnen", use_container_width=True)

with right:
    st.subheader("Ergebnis")
    if run or st.session_state.get("run_calculation", False):
        if run:
            st.session_state["run_calculation"] = True
        try:
            overview_df, formulas_df = calculate_lambda(fluid, float(temperature_c), float(pressure_bar), float(diameter_mm), float(velocity_ms), float(roughness_mm))
            st.dataframe(overview_df, use_container_width=True, hide_index=True, height=(len(overview_df) + 1) * 35 + 3)
            with st.expander("Formeln im Vergleich"):
                st.dataframe(formulas_df, use_container_width=True, hide_index=True, height=(len(formulas_df) + 1) * 35 + 3)
            csv_export = pd.concat(
                [
                    overview_df.assign(Bereich="Ergebnis"),
                    formulas_df.rename(columns={"Formel": "Parameter", "Rohrreibungszahl λ [-]": "Wert"}).assign(Bereich="Formelvergleich"),
                ],
                ignore_index=True,
            )
            csv_content = f"{APP_TITLE};{APP_VERSION}\n" + csv_export.to_csv(index=False, sep=";")
            st.download_button(
                label="CSV herunterladen",
                data=csv_content.encode("utf-8-sig"),
                file_name="Rohrreibungszahl-Bestimmung.csv",
                mime="text/csv",
                use_container_width=True,
            )
        except Exception as e:
            st.error(f"Fehler bei der Berechnung: {e}")
    else:
        st.info("Eingaben setzen und auf Berechnen klicken.")

st.markdown("---")
with st.expander("Anleitung"):
    st.markdown(
        f"""
Mit diesem Tool wird die **Rohrreibungszahl λ** auf Grundlage der übergebenen Stoff- und Strömungsdaten berechnet.

Berücksichtigt werden:
- Gesetz von Hagen-Poiseuille,
- Formel nach Blasiu,
- Formel nach Prandtl,
- sowie Colebrook-White,
- Formel nach Nikuradse.

Repository: https://github.com/div-ne/Bestimmung-Rohrreibungszahl/
"""
    )
with st.expander("Rohrrauheitswerte"):
    st.caption("Orientierungswerte für die Eingabe von k [mm].")
    roughness_df = pd.DataFrame(ROUGHNESS_ROWS, columns=["Werkstoff und Rohrart", "Zustand der Rohre", "k [mm]"])
    st.dataframe(roughness_df, use_container_width=True, hide_index=True)
with st.expander("Moody-Diagramm"):
    st.caption("Zur Orientierung für Strömungsbereiche und relative Rauheit.")
    if os.path.exists(MOODY_IMAGE):
        st.image(MOODY_IMAGE, use_container_width=True)
    st.markdown(f"Quelle: [Wikipedia – Moody chart]({MOODY_SOURCE})")
