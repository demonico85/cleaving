###%***************************************************************************%
###%    Created by:                                                           *%
###%    Nicodemo Di Pasquale                                                  *%
###%    Department of Mathematics                                             *%
###%    University of Leicester                                               *%
###%    24 March 2018                                                         *%
###%***************************************************************************%




print       "Wall distance for solid ${zwc} and liquid ${zwl}"

fix  fl2 FreeL wallforce ${eps} ${sigma} ${zwl} ${delta} ${rw} ${clwall1} file walls.lmp
fix  fl3 FreeL wallforce ${eps} ${sigma} ${zwl} ${delta} ${rw} ${clwall2} file walls.lmp
fix  fc2 FreeC wallforce ${eps} ${sigma} ${zwc} ${delta} ${rw} ${clwall1} file walls.lmp
fix  fc3 FreeC wallforce ${eps} ${sigma} ${zwc} ${delta} ${rw} ${clwall2} file walls.lmp

variable totWC   equal "f_fc2 + f_fc3"
variable totWL   equal "f_fl2 + f_fl3"

variable totW equal "f_fc2 + f_fc3 + f_fl2 + f_fl3"


run   ${eqnts}

# Density
fix fd all ave/chunk ${nts} 1 ${nts} binchunk density/number file data/dens-liq.${cnt}.out
#
fix  f6 all ave/time  ${Nevery} ${Nrepeat} ${Nfreq}  c_thermo_temp c_thermo_pe v_totW v_zwl  v_totWC v_totWL f_fc2 f_fc3 f_fl2 f_fl3 file  out/ave.F.liq.${cnt}.out

run ${nts}

unfix fl2
unfix fl3
unfix fc2
unfix fc3
unfix f6
unfix fd
 
variable    zwl equal ${zzw}+${dew}
variable    zzw equal ${zwl}


#%*****************************************************************************%
