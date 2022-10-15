source ~/.config/.tianff
source ${homedir}/.local/bin/bashrc_plm.sh
plumed sum_hills --hills ../../HILLS --stride 20000 --outfile fes.dist_vp_c. --min 0 --max 14 --bin 400 --negbias
