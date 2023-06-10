# Aplikacja do przeliczania współrzędnych pomiędzy układami

### Założenia projektu
Głównym celem projektu jest szybsze i wygodniejsze przeliczanie współrzędnych pomiędzy danymi układami.

# Co robi aplikacja:
-Przeliczenie współrzednyvh geocentrycznych (X,Y,Z) na współrzedne geodezyjne (φ, λ, h)

-Przeliczenie współrzednych geodezyjnych (φ, λ, h) na współrzedne geocentryczne (X,Y,Z)

-Wyznaczenie współrzednych topocentrycznych (N,E,U)

-Wyznaczenie współrzednych w układzie PL-2000

-Wyznaczenie współrzednych w układzie PL-1992

## Jakie elipsy obsługuje:
-wrs84
-grs80
-ewKrasowskiego
## Definicje wykorzystane w skrypcie 
### Definicja __init__
Definicja służy nam do zdefiniowania odpowiedniego systemu odneisienia. Zawarliśmy tutaj trzy podstawowe systemy:  wgs84, grs80 oraz ewKrasowskiego. Oraz daliśmy funkcję która sami możemy zdefiniować.
### Definicja dms
Służy zamianie stopni dziesiętnych na stopnie, minutty,sekundy.
### Definicja fromdms
Funkcja zamieniająca stopnie w układzie (d m s) na radiany i stopnie dziesiętne.
### Definicja sigma
Funkcja licząca sigmę wykorzystywaną w funkcjach fl2pl1992 oraz BL2200.
### Definicja fl2pl1992
Funkcja przeliczająca współrzędne geodezyjne na współrzędne w układzie 1992.Początkiem układu jest punkt przecięcia południka 19°E z obrazem równika. Południk środkowy odwzorowuje się na linię prostą w skali m0 = 0,9993, na południku środkowym zniekształcenie wynosi –70 cm/km i rośnie do +90 cm/km na skrajnych wschodnich obszarach Polski. Tą definicje zaleca się używać miedzy wartościa 15°E a 24°E ze względu że ten układ jest dostsowany do powierzchni Polski. Stosowanie jej na warościach większych niż wymienonych powyżej jest obarczone dużym ryzykiem błędu.
### Definicja BL22000
Funkcja przeliczająca współrzędne geodezyjne na współrzędne w układzie 2000. Charakteryzuje się tym że powierzchnia Polski jest podzielona na cztery trzystopniowe strefy o południkach 15°E, 18°E, 21°E i 24°E, oznaczone odpowiednio numerami – 5, 6, 7 i 8. Skala długości odwzorowania na południkach osiowych wynosi m0 = 0,999923. Zniekształcenia na południku osiowym wynoszą −7,7 cm/km zaś na styku stref +7 cm/km. W definicji wpisanie innej strefy skutkuje tym że program nie zadziała.
### Definicja hirvonen
Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z) na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.
### Definicja NEU
Funckja obliczająca wektor w układzie neu
### Definicja XYZ2neu
Funkcja przeliczająca współrzędzne kartezjańskie X,Y,Z do układu NEU przy pomocy wektora przestrzennego w układzie topocentrycznym.
# Wprowadzanie danych
## Za pomocą pythona
```
def hirvonen(self,X,Y,Z):
```
Aby wywołać dowolną funkcję należy podać jej argumenty które znajdują się w nawiasie, w tym przypadku współrzędne X,Y,Z. W przypadku innych funkcji mogą to być:
f - szerokość geodezyjna
l - długość geodezyjna
h - wysokość elipsoidalna
s - odległość
alfa - kąt poziomy
z - kąt pionowy
## Za pomocą pliku .txt
Plikiem wejściowym jest plik wsp_inp.txt ze współrzędnymi X,Y,Z
![image](https://user-images.githubusercontent.com/129080884/235367364-26d0a6fb-2402-4695-9468-fefd9bd603a4.png)

Aby otrzymać wyniki należy wprowadzić współrzędne jak na obrazku wyżej, ważne jest aby zgadzała się liczba kolumn(3) i wierszy(12) w przeciwnym wypadku program nie wykona obliczeń. Gdybyśmy chcieli wykonać obliczenia z pliku txt o innej liczbie kolumn i wierszy musilibyśmy zmieniać skrypt programu.

Plik wynikowy prezentuje się następująco:
![image](https://user-images.githubusercontent.com/129080884/235369548-ee8f9fad-45e1-4aa6-bb35-fbdb1739e761.png)


## Za pomocą konsoli GITa
Na początku musimy otworzyć nasz program przy pomocy aplikacji Git Bash w tym celu klikamy prwawym przyciskiem myszy na nasz plik skrypt.py i wybieramy Git Bash Here. Następnie w otwartej już kosoli wpisujemy:
```
python skrypt.py 
```
Następnie dopisujemy interesujące nas argumenty z listy poniżej.
```
Transformacje wspolrzednych

options:
  -h, --help            show this help message and exit
  -X X                  wartosc wspolrzednej pierwszegi punktu X [m]
  -Y Y                  wartosc wspolrzednej pierwszego punktu Y [m]
  -Z Z                  wartosc wspolrzednej pierwszego punktu Z [m]
  -f F                  wartosc wspolrzednej f [° ' '']
  -l L                  wartosc wspolrzednej l [° ' '']
  -H H                  wartosc wspolrzednej H [m]
  -s S                  wartosc dlugosci miedzy dwoma punktami [m]
  -alfa ALFA            wartosc kat poziomego [° ' '']
  -z Z                  wartosc kat zenitalnego [° ' '']
  -ns NS                numer strefy w odwzorowaniu PL-2000 [° ' '']
  -model {grs80,wgs84,ewKrasowski}
                        model elipsoidy
```
Przykładowe wprowadzenie danych:
```
python skrypt.py -model grs80 -X 100.00 -Y 200.00 -Z 300.00
```

# Przykładowe wywołania funkcji
Aby przystapić do pierwszych prac na naszych funkcjach należy zainstalować biblioteki numpy i argparse w tym celu należy wpisać w konsole pythona nastepujące komendy:
```
python - pip install numpy
```
```
python - pip install argparse
```
Następnie tworzymy nowy plik i importujemy wymienione niżej biblioteki:
```
import numpy as np
import argparse
from math import &
from skrypt import *
```
Teraz tworzymy obiekt i wywołujemy klase "Transformacje", w nawiasie podajemy model elipsoidy:
```
test = Transformacje(model = "grs80")
```
Teraz, aby wywołać jakąś funkcje postępujemy w następujący sposób:
```
 phi, lam, h = geo.hirvonen(X, Y, Z)
    print(phi, lam, h)
```
W zależności od interesującej nas transforamcji zmienia się tylko nazwa funkcji i jej argumenty.

# Wymagania 
Program ten można uruchomić na komputerze z systemem operacyjnym Windows.Oraz z odpowiednim programem Phython (wersja 3.9) i musi posiadać bibioteki Numpy oraz argparse.

# Znane błędy:
Transformacja BL-->2000 oraz BL-->1992 dla elipsoidy krasowskiego daje błędne wyniki i nie powinna być używana. 
