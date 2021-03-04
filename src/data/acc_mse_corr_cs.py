import pandas as pd
from sklearn.metrics import mean_squared_error
from scipy.spatial.distance import cosine
import numpy as np





data_30hz = pd.read_csv(snakemake.input[0])
features = ['phone_accelerometer_rapids_maxmagnitude',
            'phone_accelerometer_rapids_minmagnitude',
            'phone_accelerometer_rapids_avgmagnitude',
            'phone_accelerometer_rapids_medianmagnitude',
            'phone_accelerometer_rapids_stdmagnitude',
            'phone_accelerometer_panda_sumdurationexertionalactivityepisode',
            'phone_accelerometer_panda_maxdurationexertionalactivityepisode',
            'phone_accelerometer_panda_mindurationexertionalactivityepisode',
            'phone_accelerometer_panda_avgdurationexertionalactivityepisode',
            'phone_accelerometer_panda_mediandurationexertionalactivityepisode',
            'phone_accelerometer_panda_stddurationexertionalactivityepisode',
            'phone_accelerometer_panda_sumdurationnonexertionalactivityepisode',
            'phone_accelerometer_panda_maxdurationnonexertionalactivityepisode',
            'phone_accelerometer_panda_mindurationnonexertionalactivityepisode',
            'phone_accelerometer_panda_avgdurationnonexertionalactivityepisode',
            'phone_accelerometer_panda_mediandurationnonexertionalactivityepisode',
            'phone_accelerometer_panda_stddurationnonexertionalactivityepisode']


compare_freqs = pd.DataFrame(columns=["mse_15hz", "mse_10hz", "mse_5hz", "mse_1hz", "corr_15hz", "corr_10hz", "corr_5hz", "corr_1hz", "cs_15hz", "cs_10hz", "cs_5hz", "cs_1hz"])
compare_freqs["feature"] = features
compare_freqs.set_index("feature", inplace=True)
missing_ratio_thresh = 0.75

for feature in features:
    y_30hz = data_30hz[feature]

    for idx in range(1, 5):
        file_path = snakemake.input[idx]
        data = pd.read_csv(file_path)
        freq = file_path.split("_")[-1][:-4]
        y_freq = data[feature]

        # mean squared error
        compare_freqs.loc[feature, "mse_"+freq] = ((y_30hz - y_freq)**2).sum() / y_30hz.count() if y_30hz.count() != 0 else None
        # pearson correlation
        compare_freqs.loc[feature, "corr_"+freq] = y_30hz.corr(y_freq, method="pearson", min_periods=missing_ratio_thresh*len(y_30hz))
        # cosine similarity
        compare_freqs.loc[feature, "cs_"+freq] = np.nansum(y_30hz * y_freq) / (np.sqrt((y_30hz**2).sum()) * np.sqrt((y_freq**2).sum())) if y_freq.count() >= missing_ratio_thresh*len(y_30hz) else None

compare_freqs.to_csv(snakemake.output[0])





