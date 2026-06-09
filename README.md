# Einfache Kreislaufberechnung

[**Zur Anwendung**](https://kreislaufberechnung.streamlit.app/)

## Funktionen

- Berechnung zentraler Systemdaten wie Kältemittelmassenstrom und EER.
- Berechnung von neun Kreislaufpunkten mit spezifischer Enthalpie, Entropie und Dichte.
- Bei Kältemitteln mit Temperaturgleit wird die eingegebene Verflüssigungs- bzw. Verdampfungstemperatur als mittlere Temperatur interpretiert. Für die Verflüssigung wird die mittlere Temperatur zwischen linker und rechter Grenzkurve verwendet, für die Verdampfung die mittlere Temperatur zwischen Verdampfereintritt und rechter Grenzkurve.
- Auslegung von jeweils zwei Rohrleitungsvarianten für Heißgasleitung, Flüssigkeitsleitung und Saugleitung inklusive Innenvolumen und Mantelfläche.
- Export aller berechneten Daten als CSV-Datei.

## Eingabeparameter

Die Berechnung wird mit folgenden Eingaben durchgeführt

- Projektname
- Kältemittel
- Eingabemodus
- Kälteleistung, Wärmeleistung oder Verdichtervolumenstrom, abhängig vom gewählten Modus
- Verflüssigungstemperatur
- Unterkühlung
- Verdampfungstemperatur
- Verdampferüberhitzung
- Saugleitungsüberhitzung
- Aktivierung der Rohrleitungsdimensionierung
- Leitungslängen für Heißgas-, Flüssigkeits- und Saugleitung

## Kreislaufpunkte

Folgende Kreislaufpunkte werden berechnet:
- **1:** Verdichtereingang nach Überhitzung von Verdampfer und Saugleitung
- **2:** Verdichterausgang und Verflüssigereintritt
- **3:** Verflüssigeraustritt und Expansionsventileingang
- **4:** Expansionsventilausgang und Verdampfereintritt
- **5:** Verdampferausgang und Eingang Saugleitung
- **c''**: Gesättigtes Gas bei Verflüssigungsdruck
- **c'**: Flüssigkeit bei Verflüssigungsdruck
- **0''**: Gesättigtes Gas bei Verdampfungsdruck
- **0'**: Flüssigkeit bei Verdampfungsdruck

## Rohrleitungsdimensionierung

Die Rohrleitungsberechnung verwendet hinterlegte Kupferrohrabmessungen und bewertet je Leitung Strömungsgeschwindigkeit, Jacobs-Geschwindigkeit, Druckverlust, Außenmantelfläche und Innenvolumen.

Für jede Leitung werden zwei Varianten bestimmt:

- eine Variante mit minimalem Druckverlust
- eine Variante mit minimalem Durchmesser

## CSV-Export

Über den Button **„CSV-Datei erstellen“** können die Daten exportiert werden

## Technische Basis

Die Anwendung basiert auf **Streamlit** für die Oberfläche, **CoolProp** für Stoffdaten und thermodynamische Zustandsgrößen, **NumPy** und **Pandas** für Berechnung und Datenaufbereitung sowie **SciPy** für numerische Lösungsverfahren mit `brentq`.
