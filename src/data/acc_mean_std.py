import pandas as pd

freqs = snakemake.params["freq"]

column_names = ["mean_" + str(freq) + "hz" for freq in freqs] + ["std_" + str(freq) + "hz" for freq in freqs]

compare_freqs = pd.DataFrame(columns=column_names)
for file_path in snakemake.input:
    freq = file_path.split("_")[-1][:-4]

    # mean & std
    data = pd.read_csv(file_path)
    compare_freqs["mean_" + freq] = data.mean(axis=0)
    compare_freqs["std_" + freq] = data.std(axis=0)

compare_freqs.to_csv(snakemake.output[0])

