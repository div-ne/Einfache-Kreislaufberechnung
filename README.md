# Einfache Kreislaufberechnung

[**Zur Live-Anwendung**](https://kreislaufberechnung.streamlit.app/)

## Funktionen

- Berechnung auf Basis von **Kälteleistung**, **Wärmeleistung** oder **Verdichtervolumenstrom**.
- Auswahl des Kältemittels aus einer umfangreichen Stoffdatenliste.
- Ausgabe der **Kreislaufpunkte** mit Temperatur, Druck, spezifischer Enthalpie, Dichte, spezifischer Entropie und Dampfqualität; die Punkte werden in der integrierten Anleitung erklärt.
- Optionale **Rohrleitungsdimensionierung** für Heißgas-, Flüssigkeits- und Saugleitung mit den Varianten „minimaler Druckverlust“ und „minimaler Durchmesser“.
- Export aller berechneten Ergebnisse als **CSV-Datei**, bestehend aus Systemdaten, Kreislaufpunkten und optional Rohrleitungsdimensionierung.
- Integrierte **Anleitung** mit Kurzbeschreibung der Anwendung und Erläuterung aller Kreislaufpunkte.

## Oberfläche

Die Oberfläche ist auf ein kompaktes Arbeiten ausgelegt: Eingabefelder sind in zwei Spalten gruppiert, Systemdaten werden rechts ausgegeben und die weiteren Ergebnistabellen darunter dargestellt.

Die App zeigt Titel und Versionsnummer im Kopfbereich an; in der aktuellen Vorlage ist die Version auf **0.9.1V** gesetzt.

## Eingabeparameter

Folgende Eingaben stehen in der App zur Verfügung:

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

Die aktuell gesetzten Standardwerte sind unter anderem `Projekt`, `R290`, `Unterkühlung = 2 K`, `Verdampfungstemperatur = 5 °C`, `Verdampferüberhitzung = 6 K`, `Saugleitungsüberhitzung = 9 K`, `Heissgasleitungslänge = 5 m`, `Flüssigkeitsleitungslänge = 2.5 m` und `Saugleitungslänge = 3 m`.

## Rohrleitungsdimensionierung

Die Rohrleitungsberechnung verwendet hinterlegte Kupferrohrabmessungen und bewertet je Leitung Strömungsgeschwindigkeit, Jacobs-Geschwindigkeit, Druckverlust, Außenmantelfläche und Innenvolumen.

Die Variante **minimaler Druckverlust** wählt die größte zulässige Rohrgröße innerhalb des jeweiligen Druckverlustkriteriums, während **minimaler Durchmesser** die kleinste zulässige Rohrgröße unter demselben Grenzwert auswählt.

## CSV-Export

Über den Button **„CSV-Datei erstellen“** wird eine Auswertungsdatei mit Semikolon als Trennzeichen und UTF-8-BOM erzeugt, damit sich die Datei in Tabellenprogrammen wie Excel in der Regel direkt sauber öffnen lässt.

Die CSV enthält folgende Abschnitte:

- Systemdaten
- Kreislaufpunkte
- Rohrleitungsdimensionierung, sofern aktiviert

## Technische Basis

Die Anwendung basiert auf **Streamlit** für die Oberfläche, **CoolProp** für Stoffdaten und thermodynamische Zustandsgrößen, **NumPy** und **Pandas** für Berechnung und Datenaufbereitung sowie **SciPy** für numerische Lösungsverfahren mit `brentq`.

## Hinweise

Die Anwendung ist für schnelle Überschlagsrechnungen und die kompakte Darstellung thermodynamischer Zustände gedacht.

Für die Ergebnisqualität sind plausible Eingabedaten entscheidend; insbesondere bei der Rohrleitungsdimensionierung sollten Druckverlustgrenzen, Leitungslängen und gewähltes Kältemittel fachlich passend zur betrachteten Anlage sein.
