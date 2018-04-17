set -e
rm -f *.pdf *.dat
make clean
make
./maooam
python call_ft.py
