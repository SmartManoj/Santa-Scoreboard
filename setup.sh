sudo apt update 
sudo apt-get install libgomp1 g++ tmux && tmux new -d && tmux set -g mouse on && tmux a
wget -q https://raw.githubusercontent.com/SmartManoj/Santa-Scoreboard/main/submission.csv?cb=$(date +%s) -O submission.csv
wget -q https://raw.githubusercontent.com/SmartManoj/Santa-Scoreboard/main/single_group_optimizer.cpp -O single_group_optimizer.cpp
wget -q https://raw.githubusercontent.com/SmartManoj/Santa-Scoreboard/main/run_until_converge_single_group.py -O run_until_converge_single_group.py
wget -q https://raw.githubusercontent.com/SmartManoj/Santa-Scoreboard/main/notify.py -O notify.py
g++ -O3 -march=native -std=c++17 -fopenmp -o single_group_optimizer single_group_optimizer.cpp
