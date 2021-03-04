import pandas as pd

original_acc = pd.read_csv(snakemake.input[0])
requested_freq = int(snakemake.params["freq"])

# check original frequency
original_freq = int(original_acc.groupby(["local_date_time"])[["timestamp"]].count().mean())



print("original num of rows: ", original_acc.shape[0])

# downsample to 30Hz
if original_freq == 30:
    # p02: 30Hz => 30Hz
    acc_30hz = original_acc
if original_freq == 36:
    # p05 36Hz => 30Hz: select 5 from 6 rows
    acc_30hz = original_acc[original_acc.index % 6 != 0]
if original_freq == 50:
    # p06 & p09 50Hz: select 3 rows from 5 rows
    acc_30hz = original_acc[(original_acc.index % 5 != 0) & (original_acc.index % 5 != 3)]
else:
    raise ValueError("Only support 30Hz/36Hz/50Hz currently")

acc_30hz.reset_index(inplace=True, drop=True)
if requested_freq == 30:
    print("30Hz num of rows: ", acc_30hz.shape[0])
    acc_30hz.to_csv(snakemake.output[0], index=False)

# 30Hz => 15Hz: select 1 row every 2 rows
if requested_freq == 15:
    acc_15hz = acc_30hz[acc_30hz.index % 2 == 0]
    print("15Hz num of rows: ", acc_15hz.shape[0])
    acc_15hz.to_csv(snakemake.output[0], index=False)

# 30Hz => 10Hz: select 1 row every 3 rows
if requested_freq == 10:
    acc_10hz = acc_30hz[acc_30hz.index % 3 == 0]
    print("10Hz num of rows: ", acc_10hz.shape[0])
    acc_10hz.to_csv(snakemake.output[0], index=False)

# 30Hz => 5Hz: select 1 row every 6 rows
if requested_freq == 5:
    acc_5hz = acc_30hz[acc_30hz.index % 6 == 0]
    print("5Hz num of rows: ", acc_5hz.shape[0])
    acc_5hz.to_csv(snakemake.output[0], index=False)

# 30Hz => 1Hz: select 1 row every 30 rows
if requested_freq == 1:
    acc_1hz = acc_30hz[acc_30hz.index % 30 == 0]
    print("1Hz num of rows: ", acc_1hz.shape[0])
    acc_1hz.to_csv(snakemake.output[0], index=False)

