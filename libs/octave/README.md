# Octave functions


We occasionally use specific functions from Octave (in substitution to Matlab functions we do not have access to, i.e. from Matlab toolboxes).

Those functions are listed below.

The functions are modified to adapt to Matlab syntax requirements, e.g.
- endif
- endfor
- endfunction
- ! ismatrix(x) replaced by ~ismatrix(x) etc
- != replaced by ~=


```
postpad.m
```

Downloaded from https://github.com/RoaldFre/NoiseBarrier/blob/master/matlab/postpad.m on 2022-04-14

```
prepad.m 
```

Downloaded from https://github.com/RoaldFre/NoiseBarrier/blob/master/matlab/prepad.m on 2022-04-14

```
upsample.m 
```

Downloaded from [https://github.com/fieldtrip/fieldtrip/tree/master/external/signal](https://github.com/fieldtrip/fieldtrip/tree/master/external/signal) on 2026-05-12.

Function originally from the [Signal Toolkit](https://gnu-octave.github.io/octave-signal/). Direct link: https://github.com/gnu-octave/octave-signal/blob/main/inst/upsample.m


```
xcorr.m
```

Downloaded from Octave Forge (https://sourceforge.net/p/octave/signal/ci/default/tree/inst/xcorr.m) on 2022-04-14. 
 
Signal processing package (https://octave.sourceforge.io/signal/) version 1.4.1 (2019-02-08).



