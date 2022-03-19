# cone_recon

The code is for simple reconstruction from raw data, with some additional scripts to analyze the raw data.

`recon_py` is the python implementation using `torchkbnufft`

- `recon`
    - Load in twix data `.dat`
    - Load gradient files, add dead time, convert to ADC points in kspace
    - Recon the raw data

- `simulation`
    - Load in twix data `.dat`
    - Load gradient files, add dead time, convert to ADC points in kspace
    - Recon the numerical phantom

- Check the raw twix data
    - Plot signal change in each shot, within the phase
    - `plot_k`, `plot_sig`