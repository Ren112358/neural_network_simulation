rm -rf ../data
mkdir ../data
g++ setups.cpp -o a.out
./a.out
g++ main.cpp -o b.out
./b.out
python3 plot.py
exit 0
