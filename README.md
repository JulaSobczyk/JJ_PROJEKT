Projekt 1

**Transformacje - dokumentacja**

Ten program został stworzony do przeprowadzania transformacji współrzędnych pomiędzy różnymi układami. 
Możliwe jest przeprowadzenie transformacji na elipsoidach WGS84, GRS80 oraz Krasowskiego:

**1.	XYZ (geocentryczne) -> PLH (elipsoidalne fi, lam, h)**
-  Funkcja ta konwertuje współrzędne kartezjańskie X, Y i Z na szerokość geograficzną P, długość geograficzną L oraz wysokość H nad poziomem elipsoidy,
- Argumenty XYZ (typ float) to wartości odpowiednio współrzędnej x, y, z w metrach w układzie kartezjańskim. 
- Funkcja zwraca trzy wartości typu float w postaci krotki (P, L, H), gdzie P to szerokość geograficzna w radianach, L to długość geograficzna w radianach, a H to wysokość nad poziomem elipsoidy w metrach.

**2.	PLH -> XYZ**
- Funkcja ta konwertuje współrzędne geograficzne P, L i H na współrzędne kartezjańskie X, Y i Z. 
- Argumentami metody są wartości P, L i H typu float wyrażone w radianach dla odpowiedniej elipsoidy.
- Funkcja zwraca trzy wartości typu float w postaci krotki (X, Y, Z), gdzie X, Y, Z to wartości współrzędnych w układzie kartezjańskim.

**3.	XYZ -> NEUp (topocentryczne northing, easting, up)**
- Funkcja przelicza zmiany we współrzędnych XYZ na zmiany wzdłuż północnego, wschodniego i pionowego kierunku. Wartości X, Y i Z odpowiadają współrzędnym geocentrycznym, a dX to wektor zmian we współrzędnych XYZ.
- Argumentami są zmienne XYZ typu float wyrażone w metrach w układzie kartezjańskim oraz dx, dy, dz typu float.

**4.	BL(GRS80, WGS84, Krasowski) -> PL2000**
- Funkcja przelicza współrzędne geograficzne na płaskie współrzędne PL-2000 dla wybranej elipsoidy.
- Argumentami metody są wartości B, L typu float wyrażone w stopniach dziesiętnych dla odpowiedniej elipsoidy.
- Funkcja zwraca wartości x i y typu float gdzie X, Y to wartości współrzędnych płaskich w układzie PL-2000 w metrach.

**5.	BL(GRS80, WGS84, Krasowski) -> PL1992**
- Funkcja przelicza współrzędne geograficzne na płaskie współrzędne PL-1992 dla wybranej elipsoidy.
- Argumentami metody są wartości B, L typu float wyrażone w stopniach dziesiętnych dla odpowiedniej elipsoidy.
- Funkcja zwraca wartości x i y typu float gdzie X, Y to wartości współrzędnych płaskich w układzie PL-1992 w metrach.
  
Program został napisany dla systemu operacyjnego Windows 10, w programie Spyder w języku programowania Python. Wymagania, które należy spełnić aby program działał prawidłowo:
- posiadać zainstalowany program, w którym zaimplementowana jest obsługa języka Python (wersja 3.11.2 lub 3.11.8) preferowanie Spyder (v 5.4.3 lub nowsze)
 - posiadać zainstalowany program git 
- posiadać konto na stronie GitHub
- mieć zainstalowane biblioteki argparse, numpy, math, sys

**Instrukcja użycia programu:**
1.	Pobierz kod źródłowy programu z repozytorium z GitHub
2.	Otwórz program git i pobierz folder z gałezi main
3.	W tym samym folderze, w którym masz zapisany program utwórz plik z danymi do transformacji
4.	Otwórz Wiersz poleceń
5.	W wierszu poleceń otwórz odpowiednią ścieżkę (tam gdzie masz zapisany program)
6.	Do wiersza poleceń wpisz wybraną nazwę transformacji wg schematu poniżej. (TUTAJ)
7.	Program wykona transformację i zapisze plik z wynikami do folderu, w którym pracujesz.

Problemy napotkane przy sprawdzaniu programu:
- Występował problem podczas pisania kodu z wykorzystaniem i użytkowaniem bloków <try:, except:> oraz z samym zapisem do plików.
- Pojawił się problem przy napisaniu kodu programu, w którym użytkownik sam wybiera elipsoidę w jakiej chce liczyć. Zamiast tego elipsoidy i transformacje do nich są odgórnie narzucone w postaci osobnych flag w terminalu.

