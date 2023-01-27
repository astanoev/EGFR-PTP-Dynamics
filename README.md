# EGFR-PTP-Dynamics
 Code for simulating dynamics of EGFR-PTP interaction networks and fitting to experimental data
========================

MATLAB library developed to produce the results for the paper 

[Maitreyi S. Joshi, Angel Stanoev, Birga Soetje, Jan Huebinger, Veronika Zorina, Lisaweta Rossmannek, Kirsten Michel, and Philippe IH Bastiaens. "RPTPγ is a redox-regulated suppressor of promigratory EGFR signaling." BioRxiv (2022) [doi: 10.1101/2022.06.01.494340]](https://www.biorxiv.org/content/10.1101/2022.06.01.494340v1.full.pdf).

-------------------------
Requirements
-------------------------

The code has been tested successfully with older versions of MATLAB (at least R2016b (9.1)) and with one of the latest version (R2022b (9.13)) on Windows.

Using the framework
===================

Good starting point is running the figure-generating script:

```matlab
plot_bif_fit
```

This script reproduces the main bifurcation diagrams in the paper, overlaid with the corresponding experimental data. Separate figure is plotted for each experimental condition, where the bifurcation diagrams are estimated and plotted using estimated parameters by fitting to the data points. Dose-response profile, a bifurcation diagram using 'g1' as a bifurcation parameter and a 3D-bifurcation diagram using both, 'g1' and 'EGF_EGFRtt' (liganded fraction), as parameters. Details about the model are presented in the manuscript, and the implementation is given in the 'models.egfr_ptprg_model' file.

Parameter estimation was done using the Metropolis-Hastings algorithm. Many of the model parameters estimated from the data were shared between the experimental conditions, where appropriate, as most of the biochemical constants were presumed not to be affected by the perturbations. Some of the model parameters (typically protein concentrations) were presumed to be affected by the knockout/rescue perturbations, and were thus individually estimated from the respective data sets. Details about the parameter sharing and bounds are presented in the manuscript. Fresh estimation can be performed using the following script, for example:

```matlab
df = data_fitting();
df.fit_data('results_fit');
```

Bifurcation diagrams are generated with a numerical continuation algorithm, using a predefined model. In the following example a continuation object is defined from a model object, and a bifurcation profile is estimated using the 'g1' bifurcation parameter in the [0,1] range:

```matlab
model = models.egfr_ptprg_model;
cn = dynamics.continuation(model);
cn.calc_profile('g1',[0,1],[0,1],true);
```
