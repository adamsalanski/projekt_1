# Aplikacja do przeliczania współrzędnych pomiędzy układami

### Założenia projektu
Głównym celem projektu jest szybsze i wygodniejsze przeliczanie współrzędnych pomiędzy pomiędzy danymi układami. W skrypcie "master" oraz gałęziach łączacych się z nią zostały zawarte  potrzebne tam definicje.

# Co robi aplikacja:
-Przeliczenie współrzednyvh geocentrycznych (X,Y,Z) na współrzedne geodezyjne (φ, λ, h)

-Przeliczenie współrzednych geodezyjnych (φ, λ, h) na współrzedne geocentryczne (X,Y,Z)

-Wyznaczenie współrzednych topocentrycznych (E,N,Up)

-Wyznaczenie współrzednych w układzie 2000

-Wyznaczenie współrzednych w układzie 92

##Jakie elipsy obsługuje:
-wrs80
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


#Przykładowy wygląd pliku

#Wymagania 
Program ten można uruchomić na komputerze z systemem operacyjnym Windows.Oraz z odpowiednim programem Phython (wersja 3.9) i musi posiadać bibioteki Numpy oraz argparse.
