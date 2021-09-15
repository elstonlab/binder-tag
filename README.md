# binder-tag

Initial release (v1.0.0): <br>
[![badge](https://zenodo.org/badge/209335720.svg)]()

Accompanies Liu B\*, Stone OJ\*, Pablo M\*, Herron JC, Nogueira AT, Dagliyan O, Grimm JB, Lavis LD, Elston TC, Hahn KM. <i> Biosensors based on peptide exposure show single molecule conformations in live cells. </i> Cell 2021, Forthcoming. (\*Equal contribution)


__Contents include code for ...__
- analyzing clusters from single particle tracking (SPT) maps (`cluster_analysis`) (Figure 6).
- modeling and fitting cluster kinetics (`cluster_dynamics`) (Figure 6).
- performing diffusional dynamics analysis and linking it to cluster analysis (`diffusional_dynamics`) (Figure 6).
- analyzing regulatory kinetics from tagSrc/Binder SPT co-diffusion events (`kinetics_analysis`) (Figure 7).

Note: The cluster analysis pipeline requires diffusional dynamics to be completed, since we annotate the diffusional states along clustered tracks. However, the diffusional dynamics analysis can be done separately (and is in fact the variational Bayes SPT method from the Elf lab (Persson F et al. Nature Methods 2013).
