import pandas as pd
import argparse

def load_otu(zotutab):
    otu = pd.read_table(zotutab)
    otu = otu.set_index("#OTU ID")
    otu.index.name = "OTU"
    otu.columns.name = 'seqid'
    otu = otu.transpose()
    return otu


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-z', '--zotutab',
                    help="usearch zotutab",
                    required=True)
    
    parser.add_argument('-o', '--otutable',
                    help="CSV filepath to write out structured z/otu table",
                    required=True)

    args = parser.parse_args()

    otu = load_otu(args.zotutab)
    otu.to_csv(args.otutable)