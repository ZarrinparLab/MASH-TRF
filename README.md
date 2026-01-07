# TRF-MASH

This is the code for the analysis of phenotypic, microbial (16S), and untargeted metabolomics changes in streptozotocin + high-fat diet (STAM/HFD) mice subjected to time-restricted feeding (TRF). The goal was to determine if this intervention can ameliorate the development of dysfunction-associated steatohepatitis (MASH). The phenotypic results are published and in the `JBR_2025` folder while the microbial and metabolomics is in the `MASHomics` one.

## Instructions on Data Access

### Mouse

16S microbiome data is available in Qiita [https://qiita.ucsd.edu] under the study ID 13785. Untargeted metabolomics data is available in the MassIVE database [https://massive.ucsd.edu/] under the MassIVE ID MSV000088063. Feature based molecular networking data is available on [GNPS](https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=f41b3349635c4db3951083ad68727fc2).

### Human

We used the [Caussy & Tripathy, et al. 2019](https://www.nature.com/articles/s41467-019-09455-9) twin MAFLD cohort 16S microbiome and untargeted metabolomics data for validation. 16S microbiome data is available in Qiita [https://qiita.ucsd.edu] under the study ID 11635. Untargeted metabolomics data is available in the MassIVE database [https://massive.ucsd.edu/] under the MassIVE ID MSV000082374. Feature based molecular networking data is available on [GNPS](https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=778cfcb06dec4af3bf8bc8c50bb388dc). 

## Organization of files

The two main folders are split up by publication. `JBR_2025` holds all the phenotypic related files while `MASHomics` contains the microbiome and metabolomics related files. Within these folder we have `data`, `results`, and `code`.

- `code` holds the scripted used for analysis.

- `data` has the required input files needed to run the code.

- `results` contains the outputs of running the code.

## Citations

Flores Ramos, S.*, Fogelson, K.*, Muti, V.B., Zhong, W., Hu, J., Hosseini, M., Loomba, R., Zarrinpar, A. “Time-Restricted Feeding Is Not Effective in Modulating Fibrosis in a Male MASH Model”. Journal of Biological Rhythms (2025)

Flores Ramos, S., Fogelson, K., Muti, V.B., Aron, A.T., Salido, R.A., Richter, R.A., Mannochio-Russo, H., Loomba, R., Dorrestein, P.C., Knight, R., Zarrinpar, A. “Microbial and Metabolomic early MASH biomarkers from a STAM Mouse Model”. (in prep)
