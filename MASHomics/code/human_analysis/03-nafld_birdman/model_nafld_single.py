from pkg_resources import resource_filename

import biom
from birdman import SingleFeatureModel
import numpy as np
import pandas as pd

MODEL_PATH = "/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/code/human_analysis/03-nafld_birdman/stan/negative_binomial_single.stan"
MD = pd.read_table("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/data/human_analysis/caussy_nafld_16s/NAFLD_metadata_clean_g2p.tsv",
                   sep="\t", index_col='sample_name',keep_default_na=False, na_values=['',])

class nafldModelSingle(SingleFeatureModel):
    def __init__(
        self,
        table: biom.Table,
        feature_id: str,
        beta_prior: float = 5.0,
        inv_disp_sd: float = 1.0,
        num_iter: int = 500,
        num_warmup: int = 500,
        **kwargs
    ):
        super().__init__(
            table=table,
            feature_id=feature_id,
            model_path=MODEL_PATH,
            num_iter=num_iter,
            num_warmup=num_warmup,
            **kwargs
        )


        D = table.shape[0]
        A = np.log(1 / D) 
	# build formula
        self.create_regression(formula="C(groups, Treatment('g1p'))", metadata=MD)

        param_dict = {
            "depth": np.log(table.sum(axis="sample")),
            "B_p": beta_prior,
            "inv_disp_sd": inv_disp_sd,
	    "A": A
        }
        self.add_parameters(param_dict)

        self.specify_model(
            params=["beta_var", "inv_disp"],
            dims={
                "beta_var": ["covariate"],
                "log_lhood": ["tbl_sample"],
                "y_predict": ["tbl_sample"]
            },
            coords={
                "covariate": self.colnames,
                "tbl_sample": self.sample_names,
            },
            include_observed_data=True,
            posterior_predictive="y_predict",
            log_likelihood="log_lhood"

        )
