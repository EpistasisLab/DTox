# *DTox*: Deep learning for Toxicology

**Yun Hao, Joseph D. Romano, and Jason H. Moore**

**University of Pennsylvania, Cedars-Sinai Medical Center**

DTox paper: https://doi.org/10.1016/j.patter.2022.100565

Analysis repository: https://github.com/yhao-compbio/DTox (codes and datasets to reproduce the results of DTox paper)

In drug development, a major reason for attrition is the lack of understanding of cellular mechanisms governing drug toxicity. The black-box nature of conventional classification models has limited their utility in identifying toxicity pathways. Here we developed DTox (Deep learning for Toxicology), an interpretation framework for knowledge-guided neural networks, which can predict compound response to toxicity assays and infer toxicity pathways of individual compounds. We demonstrate that DTox can achieve the same level of predictive performance as conventional models with a significant improvement in interpretability. Using DTox, we were able to rediscover mechanisms of transcription activation by three nuclear receptors, recapitulate cellular activities induced by aromatase inhibitors and PXR agonists, and differentiate distinctive mechanisms leading to HepG2 cytotoxicity. Virtual screening by DTox revealed that compounds with predicted cytotoxicity are at higher risk for clinical hepatic phenotypes. In summary, DTox provides a framework for deciphering cellular mechanisms of toxicity in silico. 

## File Structure
```
.
├── docs               	# DTox documentation
├── data                # DTox data  
├── code                # DTox code 
└── tmp                 # Temporary data ignored by git. See ".gitignore" 
```

## Setup

The conda environment for DTox is specified in `environment.yml`. To build and activate this environment, run:

```shell
# conda version 4.7.5
conda env create --file environment.yml

conda activate DTox
```

Once the conda environment is created, users can implement DTox as instructed in the [tutorial](DTox-implementation.ipynb). To open the tutorial with jupyter notebook, run 

```shell
jupyter notebook
```

