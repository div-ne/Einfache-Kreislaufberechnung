# Rohrreibungszahl Bestimmung

[**Zur Anwendung**](https://rohrreibungszahl.streamlit.app/)

## Funktionen

- Berechnung der Rohrreibungszahl λ.
- Vergleich der Rohrreibungszahlen nach Hagen-Poiseuille, Blasius, Prandtl, Colebrook-White und Nikuradse.
- Export der Ergebnisse als CSV-Datei.

## Eingabeparameter

Für die Berechnung werden folgende Eingaben benötigt:

- Fluid
- Temperatur
- Druck
- Rohrinnendurchmesser
- Strömungsgeschwindigkeit
- Rohrrauhheit

## Rohrrauheitswerte und Moody-Diagramm

Zusätzlich enthält die App einen aufklappbaren Bereich mit Rohrrauheitswerten für verschiedene Rohre und einem Moody-Diagramm zur Orientierung.

## CSV-Export

Über den Button **„CSV-Datei erstellen“** können die Daten exportiert werden.

## Technische Basis

Die Anwendung basiert auf **Streamlit** für die Oberfläche, **CoolProp** für Stoffdaten und thermodynamische Zustandsgrößen und **NumPy** sowie **Pandas** für Berechnung.
