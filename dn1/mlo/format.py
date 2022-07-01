import pandas as pd

data = pd.read_csv("co2_data.csv")
data.drop(
    columns=[
        "site_code",
        "hour",
        "minute",
        "second",
        "value_std_dev",
        "year",
        "month",
        "day",
        "nvalue",
        "latitude",
        "altitude",
        "longitude",
        "elevation",
        "intake_height",
        "qcflag",
    ],
    inplace=True,
)
data = data[data["value"] != -999.99]
data.to_csv("data_filtered.csv", index=False)
