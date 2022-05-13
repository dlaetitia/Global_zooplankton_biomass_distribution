# Global zooplankton biomass distribution

This repository contains the scripts as well as the inputs and outputs of the habitat models presented in the article Drago et al., 2022.

The script 0.1.master.R has been used to run the scripts from 1.model-fit_habitat_models.R to 2.global-7.Tara_data.R.

To comply with the size limitations of github repositories, the output files for Copepoda (0-200 m and 200-500 m) and Phaeodaria (200-500 m) were separated in three files, e.g.:
- model_full from resample 1 to 50 named Copepoda_world_0_200m_model_full_part1.rds for Copepoda (0-200 m)
- model_full from resample 51 and 100 named Copepoda_world_0_200m_model_full_part2.rds for Copepoda (0-200 m)
- the rest of the models output named Copepoda_world_0_200m_without_model_full.rds for Copepoda (0-200 m)
