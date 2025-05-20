import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.pyplot import rc_context
from matplotlib import rcParams

ad = sc.read(f"ff5a921e-6e6c-49f6-9412-ad9682d23307.h5ad")
ad.obs.to_csv("ff5a921e-6e6c-49f6-9412-ad9682d23307.metadata", sep = '\t')
ad.obsm.to_csv("ff5a921e-6e6c-49f6-9412-ad9682d23307.umap", sep = '\t')