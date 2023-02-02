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

Another more graphical way of probing the model capabilities and performing parameter estimation using the experimental data is by running the Bifurcation analysis GUI, which is our bifurcation analysis playground. The GUI is developed using the Matlab App Designer, and can be ran on more recent versions of Matlab (≥ R2018b) If only a model testing is desired, one can simply run the GUI either without parameters (in which case the predefined model will be used), or supplied with a model, and additionally a parameter set can be included as well:

```matlab
apps.app_bif('model',model);
```

The lower panel contains controls for changing the ranges and values for each parameter. Clicking on the controls will select the corresponding parameter as a bifurcation parameter. The bifurcation diagram is shown in the upper panel.

More interesting analysis can be done using the experimental datasets, that contain the dose-response experiments. The datasets and the respective estimated parameter sets can be loaded in a data_fitting object as before, which can then be supplied as an input argument to the GUI:

```matlab
apps.app_bif('data_fit',df);
```

In this case the datasets and parameter sets can be selected from the upper-right drop-down list, and the dose-response data will be overlaid with the bifurcation profile when the EGF_EGFRtt is selected as a bifurcation parameter. The loss function (RMSE in this case) is also calculated, and can be used to optimize each parameter individually.