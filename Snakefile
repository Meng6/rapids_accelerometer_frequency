configfile: "config.yaml"
include: "rules/packrat.snakefile"
include: "rules/preprocessing.snakefile"
include: "rules/features.snakefile"
include: "rules/models.snakefile"
include: "rules/reports.snakefile"
include: "rules/mystudy.snakefile" # You can add snakfiles with rules tailored to your project

rule all:
    input:
        # My study (this is an example of a rule created specifically for a study)
        expand("data/interim/{pid}/days_to_analyse_{days_before_surgery}_{days_in_hospital}_{days_after_discharge}.csv",
                            pid = config["PIDS"],
                            days_before_surgery = config["PARAMS_FOR_ANALYSIS"]["DAYS_TO_ANALYSE"]["DAYS_BEFORE_SURGERY"],
                            days_after_discharge = config["PARAMS_FOR_ANALYSIS"]["DAYS_TO_ANALYSE"]["DAYS_AFTER_DISCHARGE"],
                            days_in_hospital = config["PARAMS_FOR_ANALYSIS"]["DAYS_TO_ANALYSE"]["DAYS_IN_HOSPITAL"]),
        expand("data/processed/{pid}/targets_{summarised}.csv", 
                            pid = config["PIDS"],
                            summarised = config["PARAMS_FOR_ANALYSIS"]["SUMMARISED"]),
        expand("data/processed/{pid}/demographic_features.csv", pid=config["PIDS"]),
        # Feature extraction
        expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),
        expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["FITBIT_TABLE"]),
        expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),
        expand("data/processed/{pid}/battery_deltas.csv", pid=config["PIDS"]),
        expand("data/interim/{pid}/applications_foreground_with_datetime_with_genre.csv", pid=config["PIDS"]),
        expand("data/processed/{pid}/screen_deltas.csv", pid=config["PIDS"]),
        expand("data/processed/{pid}/plugin_google_activity_recognition_deltas.csv", pid=config["PIDS"]),
        expand("data/interim/{pid}/phone_valid_sensed_days.csv", pid=config["PIDS"]),
        expand("data/interim/{pid}/phone_sensed_bins.csv", pid=config["PIDS"]),
        expand("data/processed/{pid}/sms_{sms_type}_{day_segment}.csv",
                            pid=config["PIDS"],
                            sms_type = config["SMS"]["TYPES"],
                            day_segment = config["SMS"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/call_{call_type}_{segment}.csv",
                            pid=config["PIDS"], 
                            call_type=config["CALLS"]["TYPES"],
                            segment = config["CALLS"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/location_barnett_{segment}.csv", 
                            pid=config["PIDS"],
                            segment = config["BARNETT_LOCATION"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/bluetooth_{segment}.csv",
                            pid=config["PIDS"], 
                            segment = config["BLUETOOTH"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/activity_recognition_{segment}.csv",pid=config["PIDS"], 
                            segment = config["ACTIVITY_RECOGNITION"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/battery_{day_segment}.csv",
                            pid = config["PIDS"],
                            day_segment = config["BATTERY"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/screen_{day_segment}.csv",
                            pid = config["PIDS"],
                            day_segment = config["SCREEN"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/light_{day_segment}.csv",
                            pid = config["PIDS"],
                            day_segment = config["LIGHT"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/accelerometer_{day_segment}.csv",
                            pid = config["PIDS"],
                            day_segment = config["ACCELEROMETER"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/applications_foreground_{day_segment}.csv",
                            pid = config["PIDS"],
                            day_segment = config["APPLICATIONS_FOREGROUND"]["DAY_SEGMENTS"]),
        expand("data/raw/{pid}/fitbit_{fitbit_sensor}_with_datetime.csv",
                            pid=config["PIDS"],
                            fitbit_sensor=config["FITBIT_SENSORS"]),
        expand("data/processed/{pid}/fitbit_heartrate_{day_segment}.csv",
                            pid = config["PIDS"],
                            day_segment = config["HEARTRATE"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/fitbit_step_{day_segment}.csv",
                            pid = config["PIDS"],
                            day_segment = config["STEP"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/wifi_{segment}.csv",
                            pid=config["PIDS"], 
                            segment = config["WIFI"]["DAY_SEGMENTS"]),
        # Models
        expand("data/processed/{pid}/data_for_individual_model/{source}_{day_segment}_original.csv",
                            pid = config["PIDS"],
                            source = config["PARAMS_FOR_ANALYSIS"]["SOURCES"],
                            day_segment = config["PARAMS_FOR_ANALYSIS"]["DAY_SEGMENTS"]),
        expand("data/processed/data_for_population_model/{source}_{day_segment}_original.csv",
                            source = config["PARAMS_FOR_ANALYSIS"]["SOURCES"],
                            day_segment = config["PARAMS_FOR_ANALYSIS"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/data_for_individual_model/{source}_{day_segment}_clean.csv",
                            pid = config["PIDS"],
                            source = config["PARAMS_FOR_ANALYSIS"]["SOURCES"],
                            day_segment = config["PARAMS_FOR_ANALYSIS"]["DAY_SEGMENTS"]),
        expand("data/processed/data_for_population_model/{source}_{day_segment}_clean.csv",
                            source = config["PARAMS_FOR_ANALYSIS"]["SOURCES"],
                            day_segment = config["PARAMS_FOR_ANALYSIS"]["DAY_SEGMENTS"]),
        expand("data/processed/data_for_population_model/demographic_features.csv"),
        expand("data/processed/data_for_population_model/targets_{summarised}.csv",
                            summarised = config["PARAMS_FOR_ANALYSIS"]["SUMMARISED"]),
        expand("data/processed/output_population_model/{rows_nan_threshold}_{cols_nan_threshold}_{days_before_threshold}|{days_after_threshold}_{cols_var_threshold}/{model}/{cv_method}/{source}_{day_segment}_{summarised}_{scaler}/{result_component}.csv",
                            rows_nan_threshold = config["PARAMS_FOR_ANALYSIS"]["ROWS_NAN_THRESHOLD"],
                            cols_nan_threshold = config["PARAMS_FOR_ANALYSIS"]["COLS_NAN_THRESHOLD"],
                            days_before_threshold = config["PARAMS_FOR_ANALYSIS"]["PARTICIPANT_DAYS_BEFORE_THRESHOLD"],
                            days_after_threshold = config["PARAMS_FOR_ANALYSIS"]["PARTICIPANT_DAYS_AFTER_THRESHOLD"],
                            cols_var_threshold = config["PARAMS_FOR_ANALYSIS"]["COLS_VAR_THRESHOLD"],
                            model = config["PARAMS_FOR_ANALYSIS"]["MODEL_NAMES"],
                            cv_method = config["PARAMS_FOR_ANALYSIS"]["CV_METHODS"],
                            source = config["PARAMS_FOR_ANALYSIS"]["SOURCES"],
                            day_segment = config["PARAMS_FOR_ANALYSIS"]["DAY_SEGMENTS"],
                            summarised = config["PARAMS_FOR_ANALYSIS"]["SUMMARISED"],
                            scaler = config["PARAMS_FOR_ANALYSIS"]["SCALER"],
                            result_component = config["PARAMS_FOR_ANALYSIS"]["RESULT_COMPONENTS"]),
        # Vizualisations
        expand("reports/figures/{pid}/{sensor}_heatmap_rows.html", pid=config["PIDS"], sensor=config["SENSORS"]),
        expand("reports/figures/{pid}/compliance_heatmap.html", pid=config["PIDS"]),
        expand("reports/figures/{pid}/battery_consumption_rates_barchart.html", pid=config["PIDS"]),
        expand("reports/compliance/{pid}/compliance_report.html", pid=config["PIDS"]),
        expand("reports/figures/overall_compliance_heatmap.html"),

rule clean:
    shell:
        "rm -rf data/raw/* && rm -rf data/interim/* && rm -rf data/processed/* && rm -rf reports/figures/* && rm -rf reports/*.zip && rm -rf reports/compliance/*"