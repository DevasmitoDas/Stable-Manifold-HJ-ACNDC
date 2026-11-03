VariousFunction  (utilities used by the soft pneumatic problem)

- f_SatTanh.m, f_DSatTanh.m, f_SatTanh_inv.m
    Smooth tanh saturation sigma(u), its slope sigma'(u) and its inverse
    (paper eq. (sat)). sigma is a model of the actuator and enters the drift;
    it replaces the f_Sat / f_SatS input-saturation functions of the reference
    toolbox, which are not needed here (satFlag = 0).
- f_smooth_abs.m   smooth |x| surrogate s_eps of the Bouc-Wen law (eq. (hyst))
- f_halton.m       Halton points (seeds, collocation and validation states)
- f_percentile.m   percentile without the Statistics Toolbox
- f_settling_time.m settling time used in Table 2
- ordblk.m, reig.m real ordered (block) eigen-decomposition of F (from the
                   reference toolbox, used in P1 for Trs)
