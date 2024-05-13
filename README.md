# Projekt 1

## Transformacje - dokumentacja

Ten program został stworzony do przeprowadzania transformacji współrzędnych pomiędzy różnymi układami. 
Możliwe jest przeprowadzenie transformacji na elipsoidach WGS84, GRS80 oraz Krasowskiego:

### 1. XYZ (geocentryczne) -> PLH (elipsoidalne fi, lam, h)
-  Funkcja ta konwertuje współrzędne kartezjańskie X, Y i Z na szerokość geograficzną P, długość geograficzną L oraz wysokość H nad poziomem elipsoidy,
- Argumenty XYZ (typ float) to wartości odpowiednio współrzędnej x, y, z w metrach w układzie kartezjańskim. 
- Funkcja zwraca trzy wartości typu float w postaci krotki (P, L, H), gdzie P to szerokość geograficzna w radianach, L to długość geograficzna w radianach, a H to wysokość nad poziomem elipsoidy w metrach.

### 2. PLH -> XYZ
- Funkcja ta konwertuje współrzędne geograficzne P, L i H na współrzędne kartezjańskie X, Y i Z. 
- Argumentami metody są wartości P, L i H typu float wyrażone w radianach dla odpowiedniej elipsoidy.
- Funkcja zwraca trzy wartości typu float w postaci krotki (X, Y, Z), gdzie X, Y, Z to wartości współrzędnych w układzie kartezjańskim.

### 3. XYZ -> NEUp (topocentryczne -> northing, easting, up)
- Funkcja przelicza zmiany we współrzędnych XYZ na zmiany wzdłuż północnego, wschodniego i pionowego kierunku. Wartości X, Y i Z odpowiadają współrzędnym geocentrycznym, a dX to wektor zmian we współrzędnych XYZ.
- Argumentami są zmienne XYZ typu float wyrażone w metrach w układzie kartezjańskim oraz x_0, y_0, z_0 typu float.

### 4. BL(GRS80, WGS84, Krasowski) -> PL2000
- Funkcja przelicza współrzędne geograficzne na współrzędne płaskie PL-2000 dla wybranej elipsoidy.
- Argumentami metody są wartości B, L typu float wyrażone w stopniach dziesiętnych dla odpowiedniej elipsoidy.
- Funkcja zwraca wartości x i y typu float gdzie X, Y to wartości współrzędnych płaskich w układzie PL-2000 w metrach.

### 5. BL(GRS80, WGS84, Krasowski) -> PL1992
- Funkcja przelicza współrzędne geograficzne na współrzędne płaskie PL-1992 dla wybranej elipsoidy.
- Argumentami metody są wartości B, L typu float wyrażone w stopniach dziesiętnych dla odpowiedniej elipsoidy.
- Funkcja zwraca wartości x i y typu float gdzie X, Y to wartości współrzędnych płaskich w układzie PL-1992 w metrach.

  
## Wymagania uruchomienia programu  
Program został napisany dla systemu operacyjnego Windows 10, w programie Spyder w języku programowania Python. Wymagania, które należy spełnić aby program działał prawidłowo:
- posiadać zainstalowany program, w którym zaimplementowana jest obsługa języka Python (wersja 3.11.2 lub 3.11.8) preferowanie Spyder (v 5.4.3 lub nowsze)
 - posiadać zainstalowany program git 
- posiadać konto na stronie GitHub
- mieć zainstalowane biblioteki argparse, numpy, math, sys

  

## Instrukcja użycia programu:
1.	Pobierz kod źródłowy programu z repozytorium z GitHub
2.	Otwórz program git i pobierz folder z gałezi main
3.	W tym samym folderze, w którym masz zapisany program utwórz plik z danymi do transformacji
4.	Otwórz Wiersz poleceń
5.	W wierszu poleceń otwórz odpowiednią ścieżkę (tam gdzie masz zapisany program)
6.	Do wiersza poleceń wpisz wybraną nazwę transformacji wg schematu poniżej.
7.	Program wykona transformację i zapisze plik z wynikami do folderu, w którym pracujesz.



## Schemat obsługi programu:

Program przelicza współrzędne w trzech różnych elipsoidach (WGS84, GRS80 oraz Krasowskiego) oraz zmienia je między pięcioma różnym transformacjami (XYZ -> PLH, PLH -> XYZ, XYZ -> NEU, BL -> PL1992, BL -> PL2000). Wszystkie tansformacje wykonuje się z poziomu Wiersza poleceń, dalej nazywane skrótem 'cmd'. Aby uruchomić program, w pierwszej kolejnośc należy włączyć cmd i przejść w nim do folderu w któym znajduje się pobrany z repozytorium skrypt kodu. Następnie kolejność wpisywania transformacji w konsoli wygląda następująco:


**python JJ_komendy.py --header_lines <cyfra oznaczająca od którego miejsca w pliku z danymi, te dane się rozpoczynają po nagłówku, licząc linijki kolejno od 0>  <flaga transformacji> <nazwa pliku z danymi, włączcie z jego rozszerzeniem>**


*W terminalu wygląda to dla przykładu nastepująco*

![image](https://github.com/JulaSobczyk/JJ_projekt_repository/assets/166397896/70362744-6c41-4e94-8cd6-279a7ab87deb)

Użytkownik wybiera spośród zaprogramowanych piętnastu flag, które odpowiadają różnym transformacjom, dla różnych elipsoid. I tak, kolejno wszystkie flagi obslugiwane przez program to:

**Flagi dla elipsoidy WGS84**
- *--xyz2plh_geo* -> Transformacja współrzędnych XYZ do PLH dla elipsoidy WGS84
- *--plh2xyz_geo* -> Transformacja współrzędnych PLH do XYZ dla elipsoidy WGS84
- *--xyz2neu_geo* -> Transformacja współrzędnych XYZ do NEU dla elipsoidy WGS84
- *--blto92_geo* -> Transformacja współrzędnych PLH (BL) do PL1992 dla elipsoidy WGS84
- *--blto2000_geo* -> Transformacja współrzędnych PLH (BL) do PL2000 dla elipsoidy WGS84

**Flagi dla elipsoidy GRS80** 
- *--xyz2plh_g80* -> Transformacja współrzędnych XYZ do PLH dla elipsoidy GRS80
- *--plh2xyz_g80* -> Transformacja współrzędnych PLH do XYZ dla elipsoidy GRS80
- *--xyz2neu_g80* -> Transformacja współrzędnych XYZ do NEU dla elipsoidy GRS80
- *--blto92_g80* -> Transformacja współrzędnych PLH (BL) do PL1992 dla elipsoidy GRS80
- *--blto2000_g80* -> Transformacja współrzędnych PLH (BL) do PL2000 dla elipsoidy GRS80

**Flagi dla elipsoidy Krasowskiego**
- *--xyz2plh_kra* -> Transformacja współrzędnych XYZ do PLH dla elipsoidy Krasowskiego
- *--plh2xyz_kra* -> Transformacja współrzędnych PLH do XYZ dla elipsoidy Krasowskiego
- *--xyz2neu_kra* -> Transformacja współrzędnych XYZ do NEU dla elipsoidy Krasowskiego
- *--blto92_kra* -> Transformacja współrzędnych PLH (BL) do PL1992 dla elipsoidy Krasowskiego
- *--blto2000_kra* -> Transformacja współrzędnych PLH (BL) do PL2000 dla elipsoidy Krasowskiego

Użytkownik nie ma możliwości dodawać własnej elipsoidy, dla której program transformował by współrzędne, więc jedyne czym może się podługiwać to wybieraniem odpowiednich flag. 

Należy pilnować, aby flaga jak i nazwa pliku z danymi były poprawnie wpisane. W przeciwnym wypadku program albo nie pokaże nic, albo pokaże, że wystapił błąd i co jest jego przyczyną.

Dla przykładu, jak zachowuje się program, jeśli transformacja została wpisana poprawnie:

![image](https://github.com/JulaSobczyk/JJ_projekt_repository/assets/166397896/4937dd3b-9462-439b-aa8a-922c1add3e8c)

A tak zachowuje się, gdy pojawił się jakiś błąd w zapisie:
- gdy błędnie zapiszemy flagę -> program nie pokazuje nic
- 
  ![image](https://github.com/JulaSobczyk/JJ_projekt_repository/assets/166397896/31640da8-af72-4b81-a2e1-1798b8928e64)

- jeśli np. jest błędna nazwa pliku z danymi
- 
  ![image](https://github.com/JulaSobczyk/JJ_projekt_repository/assets/166397896/633ab701-4df8-4d6c-8993-1f879689ed72)

Uwagi końcowe przy użytkowaniu programu:
- Plik ze współrzędnymi do transformacji musi się znajdować w tym samym folderze co skryp kodu i w którym otwarty mamy 'cmd'. Jednoczeście w samym pliku współrzędne muszą być zapisane poprzez oddzielenie kolejnej współrzędnej przecinkiem oraz, aby np. współrzędne X były w jednej kolumnie, współrzędne Y w drugiej kolumnie obok i współrzędne Z w trzeciej kolumnie. Dla przykładu tak wygląda dobry zapis pliku z danymi:
  
![image](https://github.com/JulaSobczyk/JJ_projekt_repository/assets/166397896/efeb2dbd-f2a3-4e7a-9ea9-25b44a1eac31)

- Przy transformacji XYZ -> NEU w każdej z trzech elipsoid, po zapisaniu w cmd flagi, należy równiez podać współrzędne środka układu, które oddzielone muszą być spacją. Dla przykładu, wygląda to następująco w konsoli:

![image](https://github.com/JulaSobczyk/JJ_projekt_repository/assets/166397896/6fb0d8b7-1a2e-4856-a80c-ecb69e95f588)

- Plik z przeliczonymi współrzędnymi zapisuje się w tym samym folderze, w którym znajduje się skryp kodu oraz plik z danymi współrzędnymi
  
## Problemy napotkane przy sprawdzaniu programu:
- Pojawił się problem przy napisaniu kodu programu, w którym użytkownik sam wybiera elipsoidę w jakiej chce liczyć. Zamiast tego elipsoidy i transformacje do nich są odgórnie narzucone w postaci osobnych flag w terminalu.
- Jeżeli użytkownik wpisze błędnie nazwę flagi, to program nie pokazuje nic, gdzie w przypadku poprawnie napisanego wykonanego polecenia w terminalu, program pokazuje komunikat informujący o poprawnym zapisie do pliku. Trudności sprawiało napisanie warunku, który sprawdzałby, czy podana flaga jest jedną z tych dozwolonych, a jeżeli nią nie jest to zwracałoby wtedy adekwatną informację.

