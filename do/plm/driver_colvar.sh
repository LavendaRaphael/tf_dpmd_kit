source ~/.config/.tianff
source ${homedir}/.local/bin/bashrc_plm.sh
plumed driver --plumed plm_driver.in --noatoms --timestep 0.0005 --trajectory-stride 10