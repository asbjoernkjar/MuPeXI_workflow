import pandas as pd

def merging(df_list, sample_list, split_dict, output, rankEL_thresh, BA_thresh):

    aff = "Mut_MHCaffinity"
    rankEL = "Mut_MHCrank_EL"
    master_DF = None
    print(rankEL_thresh, BA_thresh)
    print(type(rankEL_thresh), type(BA_thresh))

    for mupexi_df,sample_name in zip(df_list,sample_list):
        
        tem_df = pd.read_csv(mupexi_df, sep = "\t", comment="#")
        print(sample_name,mupexi_df ,"whith unfiltered shape:", tem_df.shape)
        tem_df = tem_df[ ((tem_df[aff] < BA_thresh) | (tem_df[rankEL] < rankEL_thresh)) ]
        print("filtered:", tem_df.shape)
        tem_df["Sample_ID"] = split_dict[sample_name] #remapping split samples back to original sample name    x.1, x.2 --> x

        if master_DF is None:
            master_DF = tem_df
        else:
            master_DF = master_DF.append(tem_df)

        print("master_df shape:", master_DF.shape)

    master_DF.to_csv(output, sep = "\t", index = False)

if __name__ == "__main__":

    merging(df_list=snakemake.input.file_list ,
            sample_list=snakemake.params.sample_list,
            split_dict = snakemake.params.split_dict,
            output=snakemake.output.output_file,
            rankEL_thresh=snakemake.params.rankEL_thresh,
            BA_thresh=snakemake.params.BA_thresh )

