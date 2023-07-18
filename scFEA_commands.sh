#!/bin/sh
#day 0
python src/scFEA.py --data_dir data --input_dir CLL_rendeiro_grouped --test_file Data_CLL_D0.csv --moduleGene_file module_gene_m168.csv --stoichiometry_matrix cmMat_c70_m168.csv --output_flux_file output/CLL_grouped_D0_flux.csv --output_balance_file output/CLL_grouped_D0_balance.csv --train_epoch 500

#day 30
python src/scFEA.py --data_dir data --input_dir CLL_rendeiro_grouped --test_file Data_CLL_D30.csv --moduleGene_file module_gene_m168.csv --stoichiometry_matrix cmMat_c70_m168.csv --output_flux_file output/CLL_grouped_D30_flux.csv --output_balance_file output/CLL_grouped_D30_balance.csv --train_epoch 500

#day 120
python src/scFEA.py --data_dir data --input_dir CLL_rendeiro_grouped --test_file Data_CLL_D120.csv --moduleGene_file module_gene_m168.csv --stoichiometry_matrix cmMat_c70_m168.csv --output_flux_file output/CLL_grouped_D120_flux.csv --output_balance_file output/CLL_grouped_D120_balance.csv --train_epoch 500

#day 150
python src/scFEA.py --data_dir data --input_dir CLL_rendeiro_grouped --test_file Data_CLL_D150.csv --moduleGene_file module_gene_m168.csv --stoichiometry_matrix cmMat_c70_m168.csv --output_flux_file output/CLL_grouped_D150_flux.csv --output_balance_file output/CLL_grouped_D150_balance.csv --train_epoch 500

#day 280
python src/scFEA.py --data_dir data --input_dir CLL_rendeiro_grouped --test_file Data_CLL_D280.csv --moduleGene_file module_gene_m168.csv --stoichiometry_matrix cmMat_c70_m168.csv --output_flux_file output/CLL_grouped_D280_flux.csv --output_balance_file output/CLL_grouped_D280_balance.csv --train_epoch 500

