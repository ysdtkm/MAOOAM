set -e
rm -f *.pdf *.dat
make clean
COMPILER=ifort make
./maooam
python call_ft.py
