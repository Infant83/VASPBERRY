# Coarse standard-WAVEDER integration

These CSV, DAT and NPZ tables contain the same unrounded T=0 sheet response from the actual six-chunk 6×6 VASP 5.4.4 optical calculation. The occupied bundle is 1–123, with 192 stored spinor bands. The result is sigma_xy = 163.1092059559 e²/h throughout the sampled gap.

This value is **strongly unconverged in k sampling**: the grid overweights the narrow Gamma curvature peak. It is retained as a numerical diagnostic, not as a quantized Hall reference. The Chern invariant from direct Fukui overlaps is separately C=-1. The sign convention is sigma_xy=-(e²/h)C; no additional spin factor or factor of one half is applied.

At zero temperature, moving the chemical potential within this gap leaves the occupied bundle unchanged. The resulting flat curve does not establish quantization; the integrated value must also converge to the value implied by the Chern invariant.
