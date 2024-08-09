import pandas as pd

meta = pd.read_csv("metadata/metadata.csv")
meta = meta.set_index("seqid")
read_file_path = meta.S3_R1_path.str.split("/").str[-1]
read_basename = read_file_path.str.split(".").str[0]

meta["mouse_type"] = read_basename.str.split("-").str[1]
meta["polysacc"] = read_basename.str.split("-").str[2]
meta["replicate"] = read_basename.str.split("-").str[3].str.split("_").str[0]

meta.to_csv("metadata/metadata.csv")